# Study Definition

## Research Question

Among adults without a prior diabetes diagnosis and with HbA1c <6.5%, can HbA1c, age, BMI, and sex help identify those with fasting plasma glucose (FPG) ≥126 mg/dL and guide selection for additional fasting glucose testing?

## Data

We used CDC NHANES 2013–2014, 2015–2016, and 2017–2018. The shared dataset contains all 29,400 NHANES participants. Our analysis will use variables from five NHANES components:

- **DEMO**: age, sex, pregnancy
- **DIQ**: diabetes history
- **GHB**: HbA1c
- **GLU**: fasting glucose and fasting subsample weights
- **BMX**: body measurements (BMI)

The master table version is data-v0.1 (SHA-256 `7ab87b7f62a3ee56f059429943dc26bce6ee298bb792ae9e3c88688c56a833aa`, recorded in `release_metadata.json`).

## Eligible population

The modeling sample includes participants with `cohort_b = True` and `LBXGH < 6.5` in the shared master table `master_v0.1.csv`.

Cohort B includes adults aged ≥20 years who explicitly reported no prior diabetes diagnosis, were not known to be pregnant, completed the medical examination with a positive examination weight and complete survey design information, and had valid HbA1c and harmonized FPG measurements and a positive fasting subsample weight.

Exclude participants with HbA1c ≥6.5% as they already meet the HbA1c diabetes threshold (American Diabetes Association Professional Practice Committee for Diabetes, 2026, Table 2.1).

Expected counts:

- Cohort B: n = 5,658 (`data/processed/v0.1/qc/sample_flow.csv`)
- Modeling sample: 5,490 participants, of whom 170 have FPG ≥126 mg/dL

## Inputs

HbA1c, age, BMI, and sex.

| Input | Column | Unit / coding |
|---|---|---|
| HbA1c | `LBXGH` | %, continuous |
| Age | `RIDAGEYR` | years, continuous; 80 represents age ≥80 |
| BMI | `BMXBMI` | kg/m², continuous (56 missing) |
| Sex | `RIAGENDR` | 1 = male, 2 = female; treated as categorical |

Within each cross-validation fold, missing BMI values will be imputed using the median BMI from the training portion. The same median will then be applied to the validation portion.

All are routinely available without fasting.

A screening score developed and validated using NHANES data defined undiagnosed diabetes as FPG ≥126 mg/dL in participants without known diabetes. Older age, male sex, and obesity were independently associated with this outcome in the final multivariable model (Bang et al., 2009). These findings support including age, sex, and BMI as candidate predictors in our study.

A study of diabetic male patients in Saudi Arabia found that despite fasting blood glucose (FPG) ≥126 mg/dL, some patients' HbA1c remained below 6.5%, with many having HbA1c values between 6.0% and 6.5% (Khan et al., 2014). Therefore, HbA1c is included in this study as a continuous predictor to evaluate whether variation below the diagnostic threshold helps identify participants with FPG ≥126 mg/dL.

## Outcome

Whether a participant's FPG meets the ADA diabetes threshold at the same NHANES visit (American Diabetes Association Professional Practice Committee for Diabetes, 2026, Table 2.1):

- 1 = FPG ≥126 mg/dL
- 0 = FPG <126 mg/dL

FPG values come from `fpg_harmonized`. This indicates whether FPG reaches the diabetes-range threshold.

Participants with outcome = 1 have FPG ≥126 mg/dL but HbA1c <6.5%, so they are missed when only HbA1c ≥6.5% is used.

Our model will use HbA1c, age, BMI, and sex to estimate how likely each participant is to have outcome = 1.

The outcome will be used only to train and evaluate the model, never as an input.

## Model

Logistic regression will use the four inputs to estimate each participant's probability of FPG ≥126 mg/dL (the probability between 0 and 1).

**Alert cutoff.** A participant will receive an alert if their probability is above a cutoff. Role C will use the training data to select a cutoff that flags a similar proportion of participants as the HbA1c ≥6.0% rule, then apply it to the validation data.

## Comparator rule

Suggest an FPG test for participants with HbA1c ≥6.0% and <6.5%.

**Reason.** In a US cohort without known diabetes, the risk of being subsequently diagnosed with diabetes is higher when the HbA1c level is between 6.0% and 6.5% (Selvin et al., 2010). Among Saudi men with type 2 diabetes, many of those with FPG ≥126 mg/dL but HbA1c <6.5% had HbA1c 6.0–6.4% (Khan et al., 2014).

## Primary comparison

Our model and the comparator rule will be compared at a similar alert rate. Sensitivity and PPV are the main measures. Secondary measures are specificity, alert rate, and model AUC. Model performance will be measured only on participants the model was not trained on.

| Measure | Definition |
|---|---|
| Sensitivity | TP / (TP + FN) |
| Specificity | TN / (TN + FP) |
| PPV | TP / (TP + FP) |
| Alert rate | (TP + FP) / all participants |
| AUC (model only) | probability that a random participant with outcome = 1 receives a higher model score than a random participant with outcome = 0 |

## Validation

We will test the model with five-fold cross-validation. We will split the participants into five groups, train the model on four, and test it on the fifth, repeating this five times. Stratified cross-validation will be used to maintain a similar proportion of participants with FPG ≥126 mg/dL across folds.

Validation is internal only. No independent external dataset will be used.

## Intended use

This prototype is designed to help clinicians decide whether to order an FPG test. It is for adults aged ≥20 years with no known diabetes, who are not known to be pregnant, and whose HbA1c is below 6.5%. The tool estimates the probability that FPG is ≥126 mg/dL. If this probability is above the alert cutoff, it prompts the clinician to consider an FPG test for further diabetes evaluation. It does not diagnose diabetes.

## Hypothesis

The model will identify more participants with FPG ≥126 mg/dL than the comparator rule at a similar alert rate.

## References

American Diabetes Association Professional Practice Committee for Diabetes, Bajaj, M., McCoy, R.G., Balapattabi, K., Bannuru, R.R., Bellini, N.J., Bennett, A.K., Beverly, E.A., Briggs Early, K., ChallaSivaKanaka, S., Echouffo-Tcheugui, J.B., Everett, B.M., Garg, R., Laffel, L.M., Lal, R., Matfin, G., Pandya, N., Pekas, E.J., Peters, A.L., Pilla, S.J., Romeo, G.R., Rosas, S.E., Segal, A.R., Simmons, K.M., Szmuilowicz, E.D. and ElSayed, N.A. (2026) '2. Diagnosis and classification of diabetes: standards of care in diabetes—2026', *Diabetes Care*, 49(Supplement_1), pp. S27–S49. Available at: https://doi.org/10.2337/dc26-S002.

Bang, H., Edwards, A.M., Bomback, A.S., Ballantyne, C.M., Brillon, D., Callahan, M.A., Teutsch, S.M., Mushlin, A.I. and Kern, L.M. (2009) 'Development and validation of a patient self-assessment score for diabetes risk', *Annals of Internal Medicine*, 151(11), pp. 775–783. Available at: https://doi.org/10.7326/0003-4819-151-11-200912010-00005.

Khan, H.A., Ola, M.S., Alhomida, A.S., Sobki, S.H. and Khan, S.A. (2014) 'Evaluation of HbA1c criteria for diagnosis of diabetes mellitus: a retrospective study of 12 785 type 2 Saudi male patients', *Endocrine Research*, 39(2), pp. 61–65. Available at: https://doi.org/10.3109/07435800.2013.828740.

Selvin, E., Steffes, M.W., Zhu, H., Matsushita, K., Wagenknecht, L., Pankow, J., Coresh, J. and Brancati, F.L. (2010) 'Glycated hemoglobin, diabetes, and cardiovascular risk in nondiabetic adults', *New England Journal of Medicine*, 362(9), pp. 800–811. Available at: https://doi.org/10.1056/NEJMoa0908359.
