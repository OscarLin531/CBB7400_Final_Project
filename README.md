# HbA1c and Fasting Glucose Screening

CBB 7400 Final Project

We want to understand how many adults with diabetes-range lab results would be missed by HbA1c screening alone. Using NHANES 2013–2018, we will compare HbA1c and fasting plasma glucose (FPG) among adults aged 20 or older who report no prior diabetes diagnosis. We will also explore whether HbA1c screen positivity is associated with insurance and access to healthcare.

A positive screening result does not confirm diabetes.

## Data

The data come from [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/Default.aspx). We use eight files from each of the 2013–2014, 2015–2016, and 2017–2018 cycles:

- `DEMO`: demographics and survey design
- `DIQ`: diabetes history
- `GHB`: HbA1c
- `GLU`: fasting glucose and fasting subsample weights
- `FASTQX`: fasting status
- `HIQ`: health insurance
- `HUQ`: healthcare use and access
- `BMX`: body measurements

Original XPT files belong in `data/raw/`, with their documentation in `data/codebooks/`. The download manifest records the source URLs and file checksums. We keep the raw data unchanged and link participants using `SEQN`.

## Getting started

Create a Python environment and install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

If the raw data are already included in your copy of the repository, you can skip downloading them. Otherwise, run:

```bash
python download_nhanes.py
```

The script downloads 24 data files and their codebooks. It also creates CSV copies for convenience; these copies are excluded from Git.

## Planned analysis

We will classify participants with both test results into four groups: both negative, HbA1c-only positive, FPG-only positive, and both positive. The screening thresholds are HbA1c ≥6.5% and FPG ≥126 mg/dL.

Our main outcome is the weighted proportion who are FPG-positive but HbA1c-negative. We will also estimate test disagreement and examine healthcare access. The analysis will account for the NHANES survey design and differences in glucose measurement across cycles.

The analysis notebook, figures, two-page report, slides, and a simple follow-up alert prototype will be added as the project progresses. These dependencies cover downloading data and basic notebook work; survey-analysis dependencies will be added when that implementation is ready.
