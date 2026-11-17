# NHANES data cleaning

This folder prepares one shared dataset for the team. It keeps all 29,400 DEMO
participants and adds clear flags for the two analysis samples. It does not
calculate survey-weighted estimates or diagnose diabetes.

## Run

From the repository root:

```bash
python3 -m pip install -r data_cleaning/requirements.txt
python3 data_cleaning/clean_data.py
python3 -m unittest discover -s data_cleaning/tests -v
```

The script also works from another working directory when called by absolute
path. It reads the existing local XPT files and does not download anything.

## Files in this folder

- `clean_data.py`: reading, merging, cleaning, sample flags, and QC exports.
- `variable_dictionary.csv`: definition and processing rule for every master column.
- `PROCESSING.md`: cleaning decisions, sample definitions, and QC notes.
- `SOURCES.md`: official documentation links for each cycle and module.
- `tests/`: boundary, missing-data, zero-value, and actual-data checks.

## Outputs

Outputs are written to `data/processed/v0.1/`:

- `master_v0.1.csv`: full participant table; generated locally and ignored by Git.
- `release_metadata.json`: master checksum, source checksums, code/dictionary
  checksums, and Python package versions.
- `qc/`: small review tables tracked with the project.

CSV missing values are written as `NA`. Boolean columns use `True` and `False`;
measurement flags can also be `NA`. In R, read with `na.strings = "NA"`. Keep
`cycle` and categorical labels as categories for analysis, not continuous values.
There are no separate A/B exports: use the flags on the complete master table.

To check a revision without replacing the default outputs:

```bash
python3 data_cleaning/clean_data.py --output-dir /tmp/nhanes-review
```

The output location does not change the dataset version. Change `VERSION` and
the matching dictionary entry before publishing a new data version. Rerunning
with the same output path replaces generated files in that folder.

## Getting the master table

The master CSV is about 9.5 MiB and is not stored in Git. R users also need to
install Python and the cleaning requirements, then run the commands above once.
Alternatively, one member can run the pipeline and share `master_v0.1.csv` with
its matching `release_metadata.json` through the team's shared drive. No shared
drive copy is uploaded by this script.

Check the CSV against `master_sha256` in that metadata file before using a shared
copy. On macOS or Linux:

```bash
shasum -a 256 data/processed/v0.1/master_v0.1.csv
```

The printed hash must match the metadata. Keep the data version and checksum
with each analysis so everyone uses the same input.

## Survey weights in R

The full master intentionally keeps missing fasting weights. There are 19,844
records with `WTSAF2YR=NA`, 1,175 with zero fasting weight, and 8,381 with a
positive fasting weight. These counts include all ages. Do not pass NA weights
to `svydesign()`.

Use this common preparation rule for the fasting design:

1. Keep records with nonmissing, positive `fasting_weight_6yr` and nonmissing
   `SDMVSTRA` and `SDMVPSU`. This selects the applicable weighted fasting sample.
2. Build `svydesign()` on that entire sample, including participants outside B.
   Use `SDMVSTRA`, `SDMVPSU`, `fasting_weight_6yr`, and `nest=TRUE`.
3. Use `subset()` on the design object to select `cohort_b`.

For the MEC design, follow the same steps with `mec_weight_6yr` and then select
`cohort_a`. Do not filter by adulthood, diabetes history, lab availability, or
A/B membership before constructing the corresponding design. Check that imported
Boolean columns are logical before using them as R subset conditions.

Keep NA weights in the shared master; do not overwrite them with zero. Filtering
nonmissing positive weights for design construction makes the distinction between
an absent weight and a published zero explicit. This follows the
[CDC R guidance](https://wwwn.cdc.gov/nchs/nhanes/tutorials/softwaretips.aspx)
to retain all nonzero-weight records before selecting the target domain.

## Handoff checklist

1. A second member reviews `PROCESSING.md`, the dictionary, and QC tables.
2. The statistics member confirms that the full table has the required weights,
   strata, PSU identifiers, and domain flags.
3. Resolve review issues and rerun the script and tests.
4. Commit the code, dictionary, QC tables, and metadata together. After review,
   create a release tag such as `data-v0.1` on that commit. No release tag has
   been created by this workflow.
5. Each analysis records the data version and `master_sha256`, not only a filename.

Python is used here for preparation. The downstream survey analysis can use R
`survey`; it should use the full master to build designs before domain subsetting.
