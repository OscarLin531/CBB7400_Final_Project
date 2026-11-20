# CDS Page Text

Text for the CDS prototype page. This page will collect HbA1c, age, BMI, and sex, then display one of the following messages.

## Inputs

| Field | Unit / options |
|---|---|
| HbA1c (LOINC 4548-4) | % |
| Age | years (≥20) |
| BMI | kg/m² |
| Sex | Male / Female |

## Messages

| Situation | Message |
|---|---|
| 1. Probability above the alert cutoff | "The predicted probability of FPG ≥126 mg/dL exceeds the alert cutoff. Consider additional fasting plasma glucose testing based on the clinical context." |
| 2. Probability at or below the cutoff | "No alert triggered. This result does not rule out abnormal glucose levels or replace routine screening and clinical judgment." |
| 3. A required input is missing | "Insufficient information. Probability cannot be calculated. Please complete all required fields." |
| 4. HbA1c ≥6.5% | "Outside the scope of this model. HbA1c is in the diabetes range. Please follow the appropriate clinical evaluation and confirmation pathway." |
| 5. Age <20, known diabetes, or known pregnancy | "Outside the scope of this model. This model is intended for adults aged ≥20 years without known diabetes who are not known to be pregnant." |

## Footer (if necessary)

"Model probabilities have not been assessed for calibration. For research demonstration only; not for clinical use."

## Implementation Notes

The alert cutoff value comes from Role C's final model.
