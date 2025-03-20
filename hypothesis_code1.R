## -------------------------------------------------------------------------------
## Stephanie S. Hong
## Capstone, March 19, 2025
## Title: Strengths and weaknesses of in silo replication of epidemiological studies:  Benzodiazepines treatment and patient outcome
## Hypothesis 1:  Patients who were acutely (within past 30 days of infection) exposed to N-unsubstituted benzodiazepines, which activate GPR68, would have worsening of acute lung injury caused by COVID-19.  
## benzodiazepines drugs of interest: alprazolam (Xanax), diazepam (Valium), lorazepam (Ativan), or clonazepam (Klonopin). 
## --------------------------------------------------------------------------------

# R Code for Survival Analysis of Benzodiazepine Exposure in COVID-19 Patients
# Retrospective Cohort Study - Analysis Part 1

# Load required libraries
library(survival)      # For survival analysis
library(survminer)     # For visualization of survival curves
library(dplyr)         # For data manipulation
library(tidyr)         # For data tidying
library(ggplot2)       # For additional plotting
library(cmprsk)        # For competing risks analysis
library(rms)           # For regression modeling
library(tableone)      # For creating descriptive tables

# Set seed for reproducibility
set.seed(123)

# Load data (replace with your actual data loading code)
# This assumes your data is in a CSV file with appropriate columns
covid_benzo_data <- read.csv("path_to_your_data.csv")

# Data preparation
# Ensure dates are in proper date format
covid_benzo_data$index_date <- as.Date(covid_benzo_data$covid_diagnosis_date)
covid_benzo_data$death_date <- as.Date(covid_benzo_data$death_date)
covid_benzo_data$last_followup_date <- as.Date(covid_benzo_data$last_followup_date)

# Create benzo exposure variables 
# (1 = exposed, 0 = not exposed in the 30 days prior to COVID diagnosis)
covid_benzo_data <- covid_benzo_data %>%
  mutate(
    any_benzo_exposure = ifelse(!is.na(benzo_prescription_date) & 
                               (index_date - benzo_prescription_date <= 30) & 
                               (index_date - benzo_prescription_date >= 0), 1, 0),
    alprazolam_exposure = ifelse(benzo_type == "alprazolam" & any_benzo_exposure == 1, 1, 0),
    diazepam_exposure = ifelse(benzo_type == "diazepam" & any_benzo_exposure == 1, 1, 0),
    lorazepam_exposure = ifelse(benzo_type == "lorazepam" & any_benzo_exposure == 1, 1, 0),
    clonazepam_exposure = ifelse(benzo_type == "clonazepam" & any_benzo_exposure == 1, 1, 0)
  )

# Create categorical variable for specific benzo exposure
covid_benzo_data <- covid_benzo_data %>%
  mutate(
    benzo_category = case_when(
      any_benzo_exposure == 0 ~ "No Exposure",
      alprazolam_exposure == 1 ~ "Alprazolam",
      diazepam_exposure == 1 ~ "Diazepam",
      lorazepam_exposure == 1 ~ "Lorazepam",
      clonazepam_exposure == 1 ~ "Clonazepam",
      TRUE ~ "Other/Multiple"
    )
  )

# Calculate follow-up time and event indicator
covid_benzo_data <- covid_benzo_data %>%
  mutate(
    # Event indicator (1 = death, 0 = censored)
    event = ifelse(!is.na(death_date), 1, 0),
    
    # Calculate follow-up time in days
    followup_time = case_when(
      event == 1 ~ as.numeric(death_date - index_date),
      TRUE ~ as.numeric(last_followup_date - index_date)
    )
  )

# Handle missing data if needed
# For this example, we'll use complete case analysis, but multiple imputation could be applied
covid_benzo_data_complete <- covid_benzo_data %>%
  filter(!is.na(followup_time) & !is.na(event) & !is.na(benzo_category))

# Create descriptive table
vars_to_summarize <- c("age", "sex", "race", "charlson_comorbidity_index", 
                      "covid_severity", "benzo_category", "event", "followup_time")

table1 <- CreateTableOne(
  vars = vars_to_summarize,
  strata = "benzo_category",
  data = covid_benzo_data_complete,
  test = TRUE
)
print(table1, showAllLevels = TRUE)

# -----------------------------
# Kaplan-Meier Analysis
# -----------------------------

# Create survival object
surv_obj <- Surv(time = covid_benzo_data_complete$followup_time, 
                event = covid_benzo_data_complete$event)

# Fit Kaplan-Meier curves by benzo category
km_fit <- survfit(surv_obj ~ benzo_category, data = covid_benzo_data_complete)

# Plot Kaplan-Meier curves
km_plot <- ggsurvplot(
  km_fit,
  data = covid_benzo_data_complete,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlim = c(0, 90),  # Limit to 90 days
  xlab = "Days from COVID-19 Diagnosis",
  ylab = "Survival Probability",
  title = "Survival by Benzodiazepine Exposure",
  legend.title = "Exposure Type",
  palette = "jco",
  risk.table.height = 0.25
)
print(km_plot)

# Log-rank test for overall difference
log_rank <- survdiff(surv_obj ~ benzo_category, data = covid_benzo_data_complete)
print(log_rank)

# -----------------------------
# Cox Proportional Hazards Model
# -----------------------------

# Simple model with just benzodiazepine exposure
cox_model1 <- coxph(surv_obj ~ benzo_category, data = covid_benzo_data_complete)
summary(cox_model1)

# Adjusted model with covariates
cox_model2 <- coxph(
  surv_obj ~ benzo_category + age + sex + charlson_comorbidity_index + covid_severity,
  data = covid_benzo_data_complete
)
summary(cox_model2)

# Forest plot of hazard ratios
ggforest(cox_model2, data = covid_benzo_data_complete)

# Test proportional hazards assumption
ph_test <- cox.zph(cox_model2)
print(ph_test)
plot(ph_test)

# If proportional hazards assumption is violated, consider stratified model
# or time-dependent covariates
if (any(ph_test$table[, "p"] < 0.05)) {
  # Example of stratified model by covid severity
  cox_model_stratified <- coxph(
    surv_obj ~ benzo_category + age + sex + charlson_comorbidity_index + strata(covid_severity),
    data = covid_benzo_data_complete
  )
  summary(cox_model_stratified)
}

# -----------------------------
# Mixed Effects Cox Model (for hospital-level clustering)
# -----------------------------

# Assuming hospital_id is available in your data
if ("hospital_id" %in% names(covid_benzo_data_complete)) {
  # Load the coxme package for mixed effects Cox models
  library(coxme)
  
  # Fit mixed effects Cox model with hospital as random effect
  mixed_cox <- coxme(
    surv_obj ~ benzo_category + age + sex + charlson_comorbidity_index + covid_severity + 
      (1 | hospital_id),
    data = covid_benzo_data_complete
  )
  summary(mixed_cox)
}

# -----------------------------
# Competing Risks Analysis (if applicable)
# -----------------------------

# Create competing risk event indicator 
# 1 = death, 2 = hospital discharge, 0 = censored
covid_benzo_data_complete <- covid_benzo_data_complete %>%
  mutate(
    competing_event = case_when(
      event == 1 ~ 1,  # Death
      discharge_status == "Discharged" ~ 2,  # Discharge
      TRUE ~ 0  # Censored
    )
  )

# Fit competing risks model
cr_model <- crr(
  ftime = covid_benzo_data_complete$followup_time,
  fstatus = covid_benzo_data_complete$competing_event,
  cov = as.matrix(covid_benzo_data_complete %>% 
                  select(age, sex, charlson_comorbidity_index, covid_severity)),
  failcode = 1,  # Death is the event of interest
  cencode = 0    # Censoring code
)
summary(cr_model)

# -----------------------------
# Sensitivity Analyses
# -----------------------------

# 1. Analysis with different exposure windows (e.g., 15 days)
covid_benzo_data <- covid_benzo_data %>%
  mutate(
    benzo_15day_exposure = ifelse(!is.na(benzo_prescription_date) & 
                                 (index_date - benzo_prescription_date <= 15) & 
                                 (index_date - benzo_prescription_date >= 0), 1, 0)
  )

# 2. Restricted mean survival time (RMST)
rmst_result <- rmst2(
  time = covid_benzo_data_complete$followup_time,
  status = covid_benzo_data_complete$event,
  arm = covid_benzo_data_complete$any_benzo_exposure,
  tau = 30  # Restriction time (e.g., 30 days)
)
print(rmst_result)

# 3. Dose-response analysis (if dose information is available)
if ("benzo_dose" %in% names(covid_benzo_data_complete)) {
  # Create dose categories (example)
  covid_benzo_data_complete <- covid_benzo_data_complete %>%
    mutate(
      dose_category = case_when(
        any_benzo_exposure == 0 ~ "No Exposure",
        benzo_dose <= low_dose_threshold ~ "Low Dose",
        benzo_dose <= medium_dose_threshold ~ "Medium Dose",
        TRUE ~ "High Dose"
      )
    )
  
  # Cox model with dose categories
  dose_model <- coxph(
    surv_obj ~ dose_category + age + sex + charlson_comorbidity_index + covid_severity,
    data = covid_benzo_data_complete
  )
  summary(dose_model)
}

# Save results
# saveRDS(list(km_fit = km_fit, cox_model = cox_model2), "survival_analysis_results.rds")

# -----------------------------
# Output tables and figures for publication
# -----------------------------

# Generate publication-ready tables
library(knitr)
library(kableExtra)

# Hazard ratios table
hr_table <- broom::tidy(cox_model2, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    "Variable" = term,
    "Hazard Ratio" = estimate,
    "95% CI Lower" = conf.low,
    "95% CI Upper" = conf.high,
    "P-value" = p.value
  )

# Print nicely formatted table
kable(hr_table, digits = 3) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

# The end
