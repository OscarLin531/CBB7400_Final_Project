# Data processing rules

Version: v0.1, ready for review. The master keeps all 29,400 participants,
with one row per `cycle` + `SEQN`. The counts below are unweighted QC counts.

## 1. Inputs and joins

We use NHANES 2013-2014, 2015-2016, and 2017-2018. Within each cycle, we
start with DEMO and left join DIQ, GHB, GLU, FASTQX, HIQ, HUQ, and BMX by
SEQN. DEMO includes everyone, so this keeps people who have no lab or
questionnaire record. We then stack the three cycles.

The script checks source SHA-256 values against `download_manifest.csv`,
nonmissing unique keys, source IDs inside DEMO, and unchanged join row counts.
It stops if a check fails. Source-presence flags tell us whether an entire
source row is absent or just an answer is missing. The original XPT files stay
unchanged; only the fields in the variable dictionary enter the master.

## 2. Zero values and survey weights

Some zeros read as about `5.3976e-79` in this environment. We found this in
weights, fasting time, HUQ051, and RIDAGEYR. The 405, 396, and 357 age repairs
restore age 0 for infants in the three cycles.

- Keep weight reader values in `WTMEC2YR_read` and `WTSAF2YR_read` for checking.
- Set nonnegative weights below `1e-10` to zero. Real positive weights in these
  files are in the thousands, far above this cutoff. Negative weights stop the
  script; missing weights stay missing.
- For other selected numeric fields, change only positive values below `1e-70`
  to zero. This restores the tiny zero artifact without rounding ordinary values.
  Counts by source and field are in `zero_repairs.csv`.
- Check each record: `WTMEC2YR=0` should agree with `RIDSTATR=1`. There are
  362, 427, and 550 interview-only participants, with no mismatches.

We keep the two-year weights and add six-year weights by dividing by three,
because we combine three regular two-year cycles. Use the normalized weights
for analysis, not the `_read` columns. Notes on NA weights and R design setup
are in the [README](README.md#survey-weights-in-r).

## 3. Missing responses and categories

| Variable | How we handle it |
| --- | --- |
| DIQ010 | Keep Yes, No, Borderline, Refused, and Don't know separate. Use `not_in_diq` for an absent source row and `item_missing` for a missing answer in an existing row. Only code 2 meets the no-prior-diagnosis rule. |
| RIDEXPRG | Exclude only 1 (known pregnancy). Keep 3 and flag it with `pregnancy_uncertain`. NA becomes `no_public_information`. Public status is restricted to females aged 20-44, so NA does not prove nonpregnancy or explain whether a question was asked. |
| HIQ011 | Keep Yes, No, Refused, Don't know, and missing separate. |
| HUQ030 | Keep the original categories. In `usual_care`, codes 1 and 3 become Yes; 2 becomes No. Keep nonresponse separate. |
| HUQ051 | Keep visit ranges as categories. Code 0 means no visits; 7 means 13-15 visits; 77 means Refused; 99 means Don't know. It is not an exact visit count. |

We handle response codes field by field: 7 or 9 does not mean missing in every
variable. We do not impute missing data or delete lab/BMI outliers. Positive,
finite lab results meet the availability checks; unusual results can still be
reviewed against the original files and documentation.

## 4. Fasting measurements and instrument bridging

We use FASTQX for `PHAFSTHR` and `PHAFSTMN` in every cycle. The GLU_H versions
stay in `PHAFSTHR_glu` and `PHAFSTMN_glu` for comparison; later GLU files do not
have them, so these columns are NA in later cycles.

In 2013-2014, both sources have values for 3,291 people and agree exactly for
each field. Another 5,891 records have a value on only one side for each field.
So the difference here is coverage, rather than different answers from the same
person. The two source-comparison exports show both types of checks.

Fasting time is an extra QC field in v0.1, not another inclusion cutoff. A
positive published fasting weight indicates eligibility for the weighted fasting
analysis. Having a glucose result alone is not enough.

Keep published `LBXGLU` in mg/dL and create `fpg_harmonized`:

- 2013-2014: `1.023 * LBXGLU - 0.5108` (CDC's forward bridge).
- 2015-2016 and 2017-2018: use `LBXGLU` as published.

CDC already forward-adjusted the 2015-2016 samples measured on C501 before
release; see the [GLU_J comparability notes](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/GLU_J.htm).
Applying the bridge again would double-correct these values. We do not use
LBDGLUSI. QC tables show threshold crossings in all available measurements and
in B separately, with raw and harmonized values in the changed-record export.

## 5. Analysis samples and screening rules

We keep everyone in the master and add eligibility flags. A required condition
with missing information gives a False eligibility flag; source fields and
response labels help explain the reason.

A (`cohort_a`) requires:

1. Age >=20 and DIQ010=2.
2. RIDEXPRG is not 1. Uncertain or unavailable public pregnancy status stays in.
3. RIDSTATR=2, positive MEC weight, and nonmissing SDMVSTRA and SDMVPSU.
4. Positive, finite HbA1c.

B (`cohort_b`) also requires positive, finite harmonized FPG and positive fasting
weight. B is therefore a subset of A.

HbA1c positivity is >=6.5%; FPG positivity is >=126 mg/dL after bridging.
Unavailable tests keep NA positivity. We assign the four paired-test categories
only inside B; others receive `not_in_cohort_b`. These are screening categories,
not confirmed clinical diagnoses.

The flow table counts exclusions in this order: DEMO -> adult -> explicit No ->
pregnancy rule -> MEC/design eligibility -> HbA1c (A) -> FPG -> fasting weight
(B). Each person counts at their first failed step. Both flow tables give the
previous, excluded, and remaining counts for each cycle and the total.

## 6. Current unweighted checks

| Check | 2013-2014 | 2015-2016 | 2017-2018 | Total |
| --- | ---: | ---: | ---: | ---: |
| All DEMO records | 10,175 | 9,971 | 9,254 | 29,400 |
| A | 4,486 | 4,234 | 4,004 | 12,724 |
| B | 2,020 | 1,831 | 1,807 | 5,658 |
| FPG present, fasting weight zero, all ages | 245 | 229 | 180 | 654 |
| FPG classification changed, all valid measurements | 20 | 0 | 0 | 20 |
| FPG classification changed within B | 8 | 0 | 0 | 8 |

These are not population estimates. If data or rules change, rerun the script
and use the updated QC exports.

## 7. QC file guide

- `file_inventory.csv`: source dimensions and checksums.
- `join_checks.csv`: row counts, IDs outside DEMO, and absent source rows.
- `zero_repairs.csv` / `weight_checks.csv`: zero repairs, missing weights,
  minimum positive weights, and MEC-status checks.
- `sample_flow.csv` / `sample_flow_long.csv`: sequential exclusions.
- `fpg_checks.csv` / `fpg_classification_changes.csv`: glucose eligibility and
  bridge effects, including records that cross the threshold.
- `fasting_source_comparison.csv` / `fasting_source_differences.csv`: agreement
  denominators and records with unequal values or different missingness.
- `missingness.csv`: missing counts by field and cycle.
- `category_counts.csv`: response reasons and unweighted phenotype counts.
