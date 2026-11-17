"""Clean and merge the NHANES files, then save the data and QC tables.

Run from any directory: python data_cleaning/clean_data.py
This script does not estimate survey-weighted results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VERSION = '0.1'
CYCLES = {'2013-2014': 'H', '2015-2016': 'I', '2017-2018': 'J'}
FIELDS = {
    'DEMO': ['SEQN', 'RIDSTATR', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH3',
             'RIDEXPRG', 'WTMEC2YR', 'SDMVSTRA', 'SDMVPSU'],
    'DIQ': ['SEQN', 'DIQ010'],
    'GHB': ['SEQN', 'LBXGH'],
    'GLU': ['SEQN', 'WTSAF2YR', 'LBXGLU'],
    'FASTQX': ['SEQN', 'PHAFSTHR', 'PHAFSTMN'],
    'HIQ': ['SEQN', 'HIQ011'],
    'HUQ': ['SEQN', 'HUQ030', 'HUQ051'],
    'BMX': ['SEQN', 'BMXBMI'],
}
WEIGHT_ZERO_TOLERANCE = 1e-10
XPT_ZERO_TOLERANCE = 1e-70


def normalize_weight(values):
    """Fix tiny weight values without changing missing values."""
    # Real positive weights here are in the thousands. The 1e-10 cutoff
    # separates them from the roughly 5e-79 values that represent zero.
    if values.lt(0).any():
        raise ValueError('Negative survey weight found')
    return values.mask(values.ge(0) & values.lt(WEIGHT_ZERO_TOLERANCE), 0.0)


def label_codes(values, mapping):
    """Map response codes to labels and flag any unexpected codes."""
    unknown = values.notna() & ~values.isin(mapping)
    if unknown.any():
        raise ValueError(f'Unexpected codes in {values.name}: {values[unknown].unique()}')
    return values.map(mapping).fillna('missing')


def derive_variables(data):
    """Add analysis columns without dropping any participant."""
    d = data.copy()
    for weight in ['WTMEC2YR', 'WTSAF2YR']:
        d[f'{weight}_read'] = d[weight]
        d[weight] = normalize_weight(d[weight])
    # We combine three regular two-year cycles, so each cycle contributes
    # one third of its two-year weight to the six-year average population.
    d['mec_weight_6yr'] = d.WTMEC2YR / 3
    d['fasting_weight_6yr'] = d.WTSAF2YR / 3
    d['diabetes_history'] = label_codes(d.DIQ010, {
        1: 'yes', 2: 'no', 3: 'borderline', 7: 'refused', 9: 'dont_know'})
    d.loc[~d.diq_record_present, 'diabetes_history'] = 'not_in_diq'
    d.loc[d.diq_record_present & d.DIQ010.isna(), 'diabetes_history'] = 'item_missing'
    d['pregnancy_status'] = label_codes(d.RIDEXPRG, {
        1: 'known_pregnant', 2: 'not_pregnant', 3: 'cannot_determine'})
    d['pregnancy_status'] = d.pregnancy_status.replace({'missing': 'no_public_information'})
    d['known_pregnancy'] = d.RIDEXPRG.eq(1)
    d['pregnancy_uncertain'] = d.RIDEXPRG.eq(3)
    d['insured'] = label_codes(d.HIQ011, {
        1: 'yes', 2: 'no', 7: 'refused', 9: 'dont_know'})
    d['usual_care_category'] = label_codes(d.HUQ030, {
        1: 'one_place', 2: 'no_place', 3: 'multiple_places', 7: 'refused', 9: 'dont_know'})
    d['usual_care'] = d.HUQ030.map({1: 'yes', 2: 'no', 3: 'yes',
                                  7: 'refused', 9: 'dont_know'}).fillna('missing')
    d['visit_frequency'] = label_codes(d.HUQ051, {
        0: 'none', 1: '1', 2: '2_to_3', 3: '4_to_5', 4: '6_to_7',
        5: '8_to_9', 6: '10_to_12', 7: '13_to_15', 8: '16_or_more',
        77: 'refused', 99: 'dont_know'})
    # Keep published mg/dL values. Only 2013-2014 needs forward bridging.
    d['fpg_harmonized'] = d.LBXGLU
    old = d.cycle.eq('2013-2014')
    d.loc[old, 'fpg_harmonized'] = 1.023 * d.loc[old, 'LBXGLU'] - 0.5108
    d['adult_eligible'] = d.RIDAGEYR.ge(20)
    d['no_prior_diabetes'] = d.DIQ010.eq(2)
    # No arbitrary outlier removal: positive finite lab values are available.
    d['hba1c_available'] = np.isfinite(d.LBXGH) & d.LBXGH.gt(0)
    d['fpg_available'] = np.isfinite(d.fpg_harmonized) & d.fpg_harmonized.gt(0)
    d['mec_weight_positive'] = np.isfinite(d.WTMEC2YR) & d.WTMEC2YR.gt(0)
    d['fasting_weight_positive'] = np.isfinite(d.WTSAF2YR) & d.WTSAF2YR.gt(0)
    d['design_available'] = d[['SDMVSTRA', 'SDMVPSU']].notna().all(axis=1)
    d['cohort_a'] = (d.adult_eligible & d.no_prior_diabetes & ~d.known_pregnancy
                     & d.RIDSTATR.eq(2) & d.mec_weight_positive
                     & d.design_available & d.hba1c_available)
    d['cohort_b'] = d.cohort_a & d.fpg_available & d.fasting_weight_positive
    # Keep missing tests as NA rather than calling them negative.
    d['hba1c_positive'] = d.LBXGH.ge(6.5).astype('boolean').where(d.hba1c_available)
    d['fpg_positive'] = d.fpg_harmonized.ge(126).astype('boolean').where(d.fpg_available)
    d['phenotype'] = pd.Series('not_in_cohort_b', index=d.index)
    a = d.hba1c_positive.fillna(False)
    f = d.fpg_positive.fillna(False)
    for mask, label in [(~a & ~f, 'neither_positive'), (a & ~f, 'hba1c_only'),
                        (~a & f, 'fpg_only'), (a & f, 'both_positive')]:
        d.loc[d.cohort_b & mask, 'phenotype'] = label
    d['fpg_classification_changed'] = pd.Series(pd.NA, index=d.index, dtype='boolean')
    valid = d.fpg_available & d.LBXGLU.notna()
    d.loc[valid, 'fpg_classification_changed'] = (
        d.loc[valid, 'LBXGLU'].ge(126) != d.loc[valid, 'fpg_harmonized'].ge(126))
    d['data_version'] = VERSION
    return d


def build_flow(data):
    """Count exclusions in order so nobody is counted twice."""
    rows = []
    groups = list(data.groupby('cycle', sort=True))
    groups.append(('total', data))

    for cycle, d in groups:
        step1 = pd.Series(True, index=d.index)
        step2 = step1 & d.adult_eligible
        step3 = step2 & d.no_prior_diabetes
        step4 = step3 & ~d.known_pregnancy
        step5 = (
            step4 & d.RIDSTATR.eq(2)
            & d.mec_weight_positive & d.design_available
        )
        step6 = step5 & d.hba1c_available
        step7 = step6 & d.fpg_available
        step8 = step7 & d.fasting_weight_positive

        stages = [
            ('01_all_demo', 'Start with all DEMO records', step1),
            ('02_adults', 'Age missing or below 20', step2),
            ('03_no_prior_diabetes', 'DIQ010 is not explicitly No (2)', step3),
            ('04_pregnancy_rule', 'Known pregnancy (RIDEXPRG=1)', step4),
            ('05_mec_design',
             'Not MEC-examined, MEC weight not positive, or design fields missing',
             step5),
            ('06_cohort_a', 'HbA1c unavailable', step6),
            ('07_fpg_available', 'FPG unavailable', step7),
            ('08_cohort_b', 'No positive fasting weight', step8),
        ]
        previous_n = len(d)
        for stage, reason, keep in stages:
            retained_n = int(keep.sum())
            rows.append({
                'cycle': cycle,
                'stage': stage,
                'exclusion_reason': reason,
                'previous_n': previous_n,
                'excluded_n': previous_n - retained_n,
                'retained_n': retained_n,
            })
            previous_n = retained_n
    return pd.DataFrame(rows)


def build_master(root):
    """Read each cycle and join its files to DEMO."""
    manifest = pd.read_csv(root / 'download_manifest.csv')
    tables, inventory, zero_checks, join_checks, sources = [], [], [], [], []
    for cycle, suffix in CYCLES.items():
        modules = {}
        for module, fields in FIELDS.items():
            path = root / 'data' / 'raw' / cycle / f'{module}_{suffix}.xpt'
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            matching_file = manifest['file'].eq(path.stem)
            matching_cycle = manifest.cycle.eq(cycle)
            record = manifest.loc[matching_file & matching_cycle]
            if len(record) != 1:
                raise ValueError(f'Manifest entry not unique: {path}')
            if digest != record.iloc[0].sha256:
                raise ValueError(f'Checksum mismatch: {path}')
            raw = pd.read_sas(path, format='xport')
            if raw.SEQN.isna().any() or raw.SEQN.duplicated().any():
                raise ValueError(f'Invalid SEQN: {path}')
            selected = list(fields)
            if module == 'GLU' and suffix == 'H':
                selected += ['PHAFSTHR', 'PHAFSTMN']
            d = raw[selected].copy()
            # Fix the tiny values that represent zero in these XPT files.
            # This includes age 0, zero visits, and zero fasting minutes.
            for column in d.columns:
                if column in ['WTMEC2YR', 'WTSAF2YR']:
                    continue  # Preserve reader values until derive_variables().
                tiny = d[column].gt(0) & d[column].lt(XPT_ZERO_TOLERANCE)
                zero_checks.append(dict(cycle=cycle, module=module, variable=column,
                                        restored_zero_n=int(tiny.sum())))
                d.loc[tiny, column] = 0.0
            d['SEQN'] = d.SEQN.astype('int64')
            if module == 'GLU' and suffix == 'H':
                d = d.rename(columns={
                    'PHAFSTHR': 'PHAFSTHR_glu',
                    'PHAFSTMN': 'PHAFSTMN_glu',
                })
            if module != 'DEMO':
                d[f'{module.lower()}_record_present'] = True
            modules[module] = d
            inventory.append(dict(cycle=cycle, module=module, rows=len(raw),
                                  columns=len(raw.columns), sha256=digest, checksum_match=True))
            sources.append(dict(cycle=cycle, module=module, source_url=record.iloc[0].source_url,
                                codebook_url=record.iloc[0].codebook_url, sha256=digest))
        # DEMO lists everyone in the cycle. Left joins keep people who have
        # no lab or questionnaire record, so we can explain missing data later.
        merged = modules['DEMO']
        for module, d in modules.items():
            if module == 'DEMO':
                continue
            unmatched = int((~d.SEQN.isin(merged.SEQN)).sum())
            if unmatched > 0:
                raise ValueError(f'{cycle} {module}: participants outside DEMO')
            overlap = (set(merged) & set(d)) - {'SEQN'}
            if overlap:
                raise ValueError(f'Unresolved column collision: {overlap}')
            before = len(merged)
            merged = merged.merge(d, on='SEQN', how='left', validate='one_to_one', sort=False)
            if len(merged) != before:
                raise ValueError(f'Join changed DEMO count: {cycle} {module}')
            flag = f'{module.lower()}_record_present'
            merged[flag] = merged[flag].eq(True)
            join_checks.append(dict(cycle=cycle, module=module, before_n=before,
                                    after_n=len(merged), unmatched_source_n=unmatched,
                                    absent_from_module_n=int((~merged[flag]).sum())))
        merged.insert(0, 'cycle', cycle)
        tables.append(merged)
    master = derive_variables(pd.concat(tables, ignore_index=True))
    if len(master) != 29400:
        raise ValueError('Expected 29,400 DEMO participants')
    if master.duplicated(['cycle', 'SEQN']).any():
        raise ValueError('Duplicate master key')
    if (master.cohort_b & ~master.cohort_a).any():
        raise ValueError('B must be a subset of A')
    return master, inventory, zero_checks, join_checks, sources


def save_sample_flow(master, qc):
    """Save both a long table and a side-by-side cycle table."""
    flow = build_flow(master)
    flow.to_csv(qc / 'sample_flow_long.csv', index=False)
    measures = ['previous_n', 'excluded_n', 'retained_n']
    cycles = list(CYCLES) + ['total']
    wide = flow.pivot(
        index=['stage', 'exclusion_reason'], columns='cycle', values=measures
    )
    columns = pd.MultiIndex.from_product([measures, cycles])
    wide = wide.reindex(columns=columns)
    wide.columns = [f'{cycle}_{measure}' for measure, cycle in wide.columns]
    wide.to_csv(qc / 'sample_flow.csv')


def save_weight_checks(master, qc):
    """Check weight repairs and compare MEC weights with examination status."""
    rows = []
    for cycle, d in master.groupby('cycle', sort=True):
        # Equal totals alone could hide different people being misclassified.
        mismatch = d.WTMEC2YR.eq(0) != d.RIDSTATR.eq(1)
        if mismatch.any():
            raise ValueError(f'MEC weight / status mismatch: {cycle}')
        for weight in ['WTMEC2YR', 'WTSAF2YR']:
            read = d[f'{weight}_read']
            repaired = read.gt(0) & read.lt(WEIGHT_ZERO_TOLERANCE)
            positive = d[weight].gt(0)
            rows.append({
                'cycle': cycle,
                'weight': weight,
                'normalized_to_zero_n': int(repaired.sum()),
                'zero_n': int(d[weight].eq(0).sum()),
                'missing_n': int(d[weight].isna().sum()),
                'min_positive': float(d.loc[positive, weight].min()),
                'interview_only_n': int(d.RIDSTATR.eq(1).sum()),
                'mec_status_mismatch_n': int(mismatch.sum()),
            })
    pd.DataFrame(rows).to_csv(qc / 'weight_checks.csv', index=False)


def save_fpg_checks(master, qc):
    """Count glucose results without eligible weights and bridge changes."""
    rows = []
    for cycle, d in master.groupby('cycle', sort=True):
        glucose_without_weight = d.LBXGLU.notna() & d.WTSAF2YR.eq(0)
        changed = d.fpg_classification_changed.fillna(False)
        changed_in_b = changed & d.cohort_b
        rows.append({
            'cycle': cycle,
            'fpg_present_zero_fasting_weight_n': int(glucose_without_weight.sum()),
            'classification_changed_all_valid_n': int(changed.sum()),
            'classification_changed_cohort_b_n': int(changed_in_b.sum()),
            'cohort_a_n': int(d.cohort_a.sum()),
            'cohort_b_n': int(d.cohort_b.sum()),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(qc / 'fpg_checks.csv', index=False)

    changed = master.fpg_classification_changed.fillna(False)
    columns = ['cycle', 'SEQN', 'LBXGLU', 'fpg_harmonized', 'cohort_b']
    master.loc[changed, columns].to_csv(
        qc / 'fpg_classification_changes.csv', index=False
    )
    return summary


def save_fasting_checks(master, qc):
    """Compare the two fasting-time sources available in 2013-2014."""
    d = master.loc[master.cycle.eq('2013-2014')]
    rows = []
    any_difference = pd.Series(False, index=d.index)
    for field in ['PHAFSTHR', 'PHAFSTMN']:
        fastqx_value = d[field]
        glu_value = d[f'{field}_glu']
        both = fastqx_value.notna() & glu_value.notna()
        different = both & fastqx_value.ne(glu_value)
        one_missing = fastqx_value.isna() != glu_value.isna()
        agreement = fastqx_value.loc[both].eq(glu_value.loc[both]).mean()
        rows.append({
            'cycle': '2013-2014',
            'variable': field,
            'both_nonmissing_n': int(both.sum()),
            'different_value_n': int(different.sum()),
            'agreement_among_both_nonmissing': float(agreement),
            'one_side_missing_n': int(one_missing.sum()),
        })
        # A missing value in one source is different from two unequal answers.
        any_difference = any_difference | different | one_missing
    pd.DataFrame(rows).to_csv(qc / 'fasting_source_comparison.csv', index=False)
    columns = [
        'cycle', 'SEQN', 'PHAFSTHR', 'PHAFSTMN', 'PHAFSTHR_glu',
        'PHAFSTMN_glu', 'WTSAF2YR', 'cohort_b',
    ]
    d.loc[any_difference, columns].to_csv(
        qc / 'fasting_source_differences.csv', index=False
    )


def save_missingness_checks(master, qc):
    """Show missing values and the response categories behind sample exclusions."""
    missing = []
    categories = []
    category_columns = [
        'diabetes_history', 'pregnancy_status', 'insured',
        'usual_care_category', 'visit_frequency', 'phenotype',
    ]
    for cycle, d in master.groupby('cycle'):
        for column in master.columns:
            missing.append({
                'cycle': cycle,
                'variable': column,
                'missing_n': int(d[column].isna().sum()),
                'total_n': len(d),
            })
        for column in category_columns:
            counts = d[column].value_counts(dropna=False)
            for value, n in counts.items():
                categories.append({
                    'cycle': cycle, 'variable': column,
                    'category': value, 'n': int(n),
                })
    pd.DataFrame(missing).to_csv(qc / 'missingness.csv', index=False)
    pd.DataFrame(categories).to_csv(qc / 'category_counts.csv', index=False)


def save_metadata(root, output, master, master_path, sources):
    """Record which data, code, and dictionary produced this master table."""
    dictionary_path = root / 'data_cleaning' / 'variable_dictionary.csv'
    dictionary = pd.read_csv(dictionary_path)
    if not dictionary.variable.is_unique:
        raise ValueError('Dictionary variables must be unique')
    if set(dictionary.variable) != set(master.columns):
        raise ValueError('Dictionary/master columns differ')

    metadata = {
        'data_version': VERSION,
        'master_file': master_path.name,
        'rows': len(master),
        'columns': len(master.columns),
        'master_sha256': hashlib.sha256(master_path.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dictionary_sha256': hashlib.sha256(dictionary_path.read_bytes()).hexdigest(),
        'pandas_version': pd.__version__,
        'numpy_version': np.__version__,
        'weight_zero_tolerance': WEIGHT_ZERO_TOLERANCE,
        'xpt_zero_tolerance': XPT_ZERO_TOLERANCE,
        'sources': sources,
    }
    metadata_path = output / 'release_metadata.json'
    metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')


def write_outputs(root, output):
    master, inventory, zero_checks, join_checks, sources = build_master(root)
    output.mkdir(parents=True, exist_ok=True)
    qc = output / 'qc'
    qc.mkdir(exist_ok=True)

    master_path = output / f'master_v{VERSION}.csv'
    master.to_csv(master_path, index=False, float_format='%.17g', na_rep='NA')
    pd.DataFrame(inventory).to_csv(qc / 'file_inventory.csv', index=False)
    pd.DataFrame(zero_checks).to_csv(qc / 'zero_repairs.csv', index=False)
    pd.DataFrame(join_checks).to_csv(qc / 'join_checks.csv', index=False)

    save_sample_flow(master, qc)
    save_weight_checks(master, qc)
    fpg_summary = save_fpg_checks(master, qc)
    save_fasting_checks(master, qc)
    save_missingness_checks(master, qc)
    save_metadata(root, output, master, master_path, sources)

    print(f'Wrote {master_path}: {len(master):,} rows, {len(master.columns)} columns')
    print(fpg_summary.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output-dir', type=Path,
        default=ROOT / 'data' / 'processed' / f'v{VERSION}',
    )
    args = parser.parse_args()
    write_outputs(ROOT, args.output_dir)


if __name__ == '__main__':
    main()
