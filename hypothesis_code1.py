## -------------------------------------------------------------------------------
## Stephanie S. Hong
## Capstone, March 19, 2025
## Title: Strengths and weaknesses of in silo replication of epidemiological studies:  Benzodiazepines treatment and patient outcome
## Hypothesis 1:  Patients who were acutely (within past 30 days of infection) exposed to N-unsubstituted benzodiazepines, which activate GPR68, would have worsening of acute lung injury caused by COVID-19.  
## benzodiazepines drugs of interest: alprazolam (Xanax), diazepam (Valium), lorazepam (Ativan), or clonazepam (Klonopin). 
## --------------------------------------------------------------------------------
# Python Code for Survival Analysis of Benzodiazepine Exposure in COVID-19 Patients
# Retrospective Cohort Study - Analysis Part 1

# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, timedelta

# Statistical and survival analysis libraries
from lifelines import KaplanMeierFitter, CoxPHFitter, CoxTimeVaryingFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats

# For handling missing data
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

# Set seed for reproducibility
np.random.seed(123)

# Data loading and preparation
# ------------------------------

# Load the data (replace with your actual data loading code)
# This assumes your data is in a CSV file with appropriate columns
df = pd.read_csv("path_to_your_data.csv")

# Convert date strings to datetime objects
df['index_date'] = pd.to_datetime(df['covid_diagnosis_date'])
df['death_date'] = pd.to_datetime(df['death_date'])
df['last_followup_date'] = pd.to_datetime(df['last_followup_date'])
df['benzo_prescription_date'] = pd.to_datetime(df['benzo_prescription_date'])

# Create benzo exposure variables (exposed within 30 days prior to COVID diagnosis)
df['any_benzo_exposure'] = ((df['index_date'] - df['benzo_prescription_date'] <= pd.Timedelta(days=30)) & 
                           (df['index_date'] - df['benzo_prescription_date'] >= pd.Timedelta(days=0))).astype(int)

# Create specific benzo exposure indicators
df['alprazolam_exposure'] = ((df['benzo_type'] == 'alprazolam') & (df['any_benzo_exposure'] == 1)).astype(int)
df['diazepam_exposure'] = ((df['benzo_type'] == 'diazepam') & (df['any_benzo_exposure'] == 1)).astype(int)
df['lorazepam_exposure'] = ((df['benzo_type'] == 'lorazepam') & (df['any_benzo_exposure'] == 1)).astype(int)
df['clonazepam_exposure'] = ((df['benzo_type'] == 'clonazepam') & (df['any_benzo_exposure'] == 1)).astype(int)

# Create categorical variable for specific benzo exposure
conditions = [
    (df['any_benzo_exposure'] == 0),
    (df['alprazolam_exposure'] == 1),
    (df['diazepam_exposure'] == 1),
    (df['lorazepam_exposure'] == 1),
    (df['clonazepam_exposure'] == 1)
]
choices = ['No Exposure', 'Alprazolam', 'Diazepam', 'Lorazepam', 'Clonazepam']
df['benzo_category'] = np.select(conditions, choices, default='Other/Multiple')

# Calculate follow-up time and event indicator
df['event'] = (~df['death_date'].isna()).astype(int)  # 1 = death, 0 = censored
df['followup_time'] = np.where(
    df['event'] == 1,
    (df['death_date'] - df['index_date']).dt.days,
    (df['last_followup_date'] - df['index_date']).dt.days
)

# Handle missing data
# For this example, we'll use complete case analysis
df_complete = df.dropna(subset=['followup_time', 'event', 'benzo_category'])

# Alternatively, use multiple imputation for missing values
# Only do this for covariates, not for outcome or exposure variables
covariates = ['age', 'sex', 'charlson_comorbidity_index', 'covid_severity']
if df[covariates].isna().sum().sum() > 0:
    print("Handling missing data in covariates...")
    imp = IterativeImputer(max_iter=10, random_state=123)
    df[covariates] = pd.DataFrame(
        imp.fit_transform(df[covariates]), 
        columns=covariates,
        index=df.index
    )

# Data exploration and descriptive statistics
# ------------------------------

# Basic descriptive statistics
desc_stats = df_complete.groupby('benzo_category')[
    ['age', 'followup_time', 'event', 'charlson_comorbidity_index']
].describe()
print(desc_stats)

# Crosstab for categorical variables
for var in ['sex', 'covid_severity', 'race']:
    if var in df_complete.columns:
        crosstab = pd.crosstab(df_complete['benzo_category'], df_complete[var], normalize='index')
        print(f"\nCrosstab for {var} by benzodiazepine category (proportions):")
        print(crosstab)

# Visualize distribution of exposure
exposure_counts = df_complete['benzo_category'].value_counts()
plt.figure(figsize=(10, 6))
sns.barplot(x=exposure_counts.index, y=exposure_counts.values)
plt.title('Distribution of Benzodiazepine Exposure')
plt.ylabel('Count')
plt.xlabel('Benzodiazepine Type')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('benzo_exposure_distribution.png', dpi=300)
plt.close()

# Kaplan-Meier Survival Analysis
# ------------------------------

# Overall survival
print("\nOverall Kaplan-Meier Survival Analysis")
kmf = KaplanMeierFitter()
kmf.fit(durations=df_complete['followup_time'], event_observed=df_complete['event'])

plt.figure(figsize=(10, 6))
kmf.plot_survival_function(ci_show=True)
plt.title('Overall Survival Probability')
plt.xlabel('Days from COVID-19 Diagnosis')
plt.ylabel('Survival Probability')
plt.grid(True)
plt.tight_layout()
plt.savefig('overall_survival.png', dpi=300)
plt.close()

# Survival by benzodiazepine category
print("\nKaplan-Meier Survival Analysis by Benzodiazepine Category")
fig, ax = plt.subplots(figsize=(12, 8))

colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#999999']
benzo_categories = df_complete['benzo_category'].unique()
for i, category in enumerate(benzo_categories):
    category_data = df_complete[df_complete['benzo_category'] == category]
    kmf = KaplanMeierFitter()
    kmf.fit(
        durations=category_data['followup_time'], 
        event_observed=category_data['event'], 
        label=category
    )
    kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[i % len(colors)])

ax.set_title('Survival by Benzodiazepine Exposure')
ax.set_xlabel('Days from COVID-19 Diagnosis')
ax.set_ylabel('Survival Probability')
ax.grid(True)
ax.set_xlim(0, 90)  # Limit to 90 days
plt.legend(loc='lower left')
plt.tight_layout()
plt.savefig('survival_by_benzo_category.png', dpi=300)
plt.close()

# Log-rank test
print("\nLog-rank test for differences in survival between exposure groups:")
groups = []
for category in benzo_categories:
    category_data = df_complete[df_complete['benzo_category'] == category]
    groups.append((
        category_data['followup_time'].values, 
        category_data['event'].values, 
        category
    ))

results = multivariate_logrank_test(
    [group[0] for group in groups],
    [group[1] for group in groups],
    [group[2] for group in groups]
)
print(results.summary)

# Pairwise log-rank tests (comparing each benzodiazepine to no exposure)
print("\nPairwise log-rank tests (each benzodiazepine type vs. no exposure):")
no_exposure = df_complete[df_complete['benzo_category'] == 'No Exposure']
for category in [cat for cat in benzo_categories if cat != 'No Exposure']:
    benzo_data = df_complete[df_complete['benzo_category'] == category]
    results = logrank_test(
        benzo_data['followup_time'], 
        no_exposure['followup_time'],
        benzo_data['event'], 
        no_exposure['event']
    )
    print(f"\n{category} vs. No Exposure:")
    print(results.summary)

# Cox Proportional Hazards Models
# ------------------------------

# Prepare data for Cox model
cox_data = df_complete.copy()

# Transform categorical variables into dummy variables
if 'sex' in cox_data.columns:
    cox_data['sex_male'] = (cox_data['sex'] == 'Male').astype(int)

# Convert benzo_category to dummy variables
benzo_dummies = pd.get_dummies(cox_data['benzo_category'], prefix='benzo', drop_first=False)
# Rename the columns to make them more interpretable
benzo_dummies.columns = [col.replace('benzo_', '') for col in benzo_dummies.columns]
cox_data = pd.concat([cox_data, benzo_dummies], axis=1)

# Drop the original categorical column
cox_data = cox_data.drop('benzo_category', axis=1)

# Standardize continuous variables (optional)
scaler = StandardScaler()
continuous_vars = ['age', 'charlson_comorbidity_index']
cox_data[continuous_vars] = scaler.fit_transform(cox_data[continuous_vars])

# Simple Cox model with just benzodiazepine exposure
print("\nSimple Cox model with just benzodiazepine exposure:")
cph1 = CoxPHFitter()
cols_to_use = [col for col in cox_data.columns if col.startswith(('Alprazolam', 'Diazepam', 'Lorazepam', 'Clonazepam'))]
cph1.fit(
    cox_data[cols_to_use + ['followup_time', 'event']], 
    duration_col='followup_time', 
    event_col='event'
)
print(cph1.summary)

# Adjusted Cox model with covariates
print("\nAdjusted Cox model with covariates:")
cph2 = CoxPHFitter()
cols_to_use = [
    'Alprazolam', 'Diazepam', 'Lorazepam', 'Clonazepam',
    'age', 'sex_male', 'charlson_comorbidity_index'
]

# Add covid_severity dummy variables if they exist
if 'covid_severity' in cox_data.columns:
    severity_dummies = pd.get_dummies(cox_data['covid_severity'], prefix='severity', drop_first=True)
    cox_data = pd.concat([cox_data, severity_dummies], axis=1)
    severity_cols = [col for col in cox_data.columns if col.startswith('severity_')]
    cols_to_use.extend(severity_cols)

cph2.fit(
    cox_data[cols_to_use + ['followup_time', 'event']], 
    duration_col='followup_time', 
    event_col='event'
)
print(cph2.summary)

# Visualize hazard ratios from the adjusted model
plt.figure(figsize=(10, 8))
cph2.plot()
plt.title('Hazard Ratios with 95% Confidence Intervals')
plt.tight_layout()
plt.savefig('hazard_ratios.png', dpi=300)
plt.close()

# Test proportional hazards assumption
print("\nTest of proportional hazards assumption:")
proportional_test = cph2.check_assumptions(
    cox_data[cols_to_use + ['followup_time', 'event']], 
    p_value_threshold=0.05, 
    show_plots=True
)

# Check for non-proportional covariates
non_proportional = proportional_test[proportional_test['p'] < 0.05]['test_statistic'].index
if len(non_proportional) > 0:
    print(f"Non-proportional covariates: {non_proportional}")
    print("Consider stratified model or time-dependent covariates")
    
    # Example of stratified Cox model (if covid_severity violates PH assumption)
    if any('severity' in var for var in non_proportional):
        print("\nFitting stratified Cox model by COVID severity...")
        # Create stratification variable
        cox_data['severity_strata'] = cox_data['covid_severity']
        
        # Fit stratified model
        strat_cols = [col for col in cols_to_use if not col.startswith('severity_')]
        cph_strat = CoxPHFitter()
        cph_strat.fit(
            cox_data[strat_cols + ['followup_time', 'event', 'severity_strata']], 
            duration_col='followup_time', 
            event_col='event',
            strata='severity_strata'
        )
        print(cph_strat.summary)
    
    # Example of time-dependent Cox model
    # Only implemented if we have time-varying data
    if 'time_varying_data.csv' in ['your_file_list.csv']:  # Replace with actual check
        print("\nFitting time-dependent Cox model...")
        # This would require time-varying data in the proper format
        # Example placeholder - replace with actual code if you have such data
        ctv = CoxTimeVaryingFitter()
        # ctv.fit(time_varying_data, id_col='patient_id', event_col='event', 
        #        start_col='start_time', stop_col='stop_time')
        # print(ctv.summary)

# Mixed Effects Cox Model
# For hospital-level clustering, Python currently has limited packages
# Sticking with stratification as a simpler alternative
if 'hospital_id' in cox_data.columns:
    print("\nStrafitied analysis by hospital (alternative to mixed effects):")
    cph_hospital = CoxPHFitter()
    cph_hospital.fit(
        cox_data[cols_to_use + ['followup_time', 'event', 'hospital_id']], 
        duration_col='followup_time', 
        event_col='event',
        strata='hospital_id'
    )
    print(cph_hospital.summary)

# Sensitivity Analyses
# ------------------------------

# 1. Analysis with different exposure windows (e.g., 15 days)
print("\nSensitivity analysis with 15-day exposure window:")
df['benzo_15day_exposure'] = ((df['index_date'] - df['benzo_prescription_date'] <= pd.Timedelta(days=15)) & 
                              (df['index_date'] - df['benzo_prescription_date'] >= pd.Timedelta(days=0))).astype(int)

# 2. Restricted mean survival time (RMST)
print("\nRestricted Mean Survival Time (RMST) Analysis:")
from lifelines import restricted_mean_survival_time

# For overall data
rmst = restricted_mean_survival_time(
    df_complete['followup_time'], 
    df_complete['event'], 
    t=30  # 30-day restriction
)
print(f"Overall RMST for 30 days: {rmst}")

# By benzo exposure
for category in benzo_categories:
    cat_data = df_complete[df_complete['benzo_category'] == category]
    rmst = restricted_mean_survival_time(
        cat_data['followup_time'], 
        cat_data['event'], 
        t=30
    )
    print(f"RMST for {category} at 30 days: {rmst}")

# 3. Dose-response analysis (if dose information is available)
if 'benzo_dose' in df_complete.columns:
    print("\nDose-response analysis:")
    # Create dose categories
    low_threshold = df_complete['benzo_dose'].quantile(0.33)
    medium_threshold = df_complete['benzo_dose'].quantile(0.67)
    
    conditions = [
        (df_complete['any_benzo_exposure'] == 0),
        (df_complete['benzo_dose'] <= low_threshold),
        (df_complete['benzo_dose'] <= medium_threshold),
        (df_complete['benzo_dose'] > medium_threshold)
    ]
    choices = ['No Exposure', 'Low Dose', 'Medium Dose', 'High Dose']
    df_complete['dose_category'] = np.select(conditions, choices, default='Unknown')
    
    # Cox model with dose categories
    dose_dummies = pd.get_dummies(df_complete['dose_category'], prefix='dose', drop_first=False)
    df_complete = pd.concat([df_complete, dose_dummies], axis=1)
    
    dose_cols = [col for col in df_complete.columns if col.startswith('dose_')]
    dose_cols = [col for col in dose_cols if 'No Exposure' not in col]  # Use No Exposure as reference
    
    cph_dose = CoxPHFitter()
    cph_dose.fit(
        df_complete[dose_cols + ['age', 'sex_male', 'charlson_comorbidity_index'] + 
                   ['followup_time', 'event']], 
        duration_col='followup_time', 
        event_col='event'
    )
    print(cph_dose.summary)

# Save models and results if needed
# import pickle
# with open('survival_models.pkl', 'wb') as f:
#     pickle.dump({
#         'cox_model': cph2,
#         'survival_data': df_complete
#     }, f)

print("Survival analysis completed.")
