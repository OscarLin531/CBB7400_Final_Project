"""Check high-impact rules and the actual NHANES handoff data."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from clean_data import ROOT, build_flow, build_master, derive_variables, normalize_weight


def example_records():
    """Small examples with known screening and eligibility outcomes."""
    n = 9
    return pd.DataFrame({
        'cycle': ['2017-2018'] * n,
        'WTMEC2YR': [100.0] * n, 'WTSAF2YR': [200.0] * n,
        'DIQ010': [2.0] * n, 'diq_record_present': [True] * n,
        'RIDEXPRG': [np.nan] * n, 'HIQ011': [1.0] * n,
        'HUQ030': [3.0] * n, 'HUQ051': [0.0] * n,
        'LBXGLU': [100., 100., 126., 126., np.nan, 140., 140., 140., 140.],
        'LBXGH': [5., 6.5, 5., 6.5, 5., 5., 5., 5., 5.],
        'RIDAGEYR': [30.] * n, 'RIDSTATR': [2.] * n,
        'SDMVSTRA': [1.] * n, 'SDMVPSU': [1.] * n,
    })


class RuleTests(unittest.TestCase):
    def test_zero_weight_and_missing_are_distinct(self):
        values = pd.Series([5.397605346934028e-79, 0., 3000., np.nan])
        result = normalize_weight(values)
        self.assertEqual(result.iloc[:3].tolist(), [0., 0., 3000.])
        self.assertTrue(pd.isna(result.iloc[3]))
        with self.assertRaises(ValueError):
            normalize_weight(pd.Series([-1.]))

    def test_thresholds_missingness_and_sample_eligibility(self):
        d = example_records()
        d.loc[5, 'WTSAF2YR'] = 5.397605346934028e-79
        d.loc[6, 'DIQ010'] = 3
        d.loc[7, 'RIDEXPRG'] = 1
        d.loc[8, 'RIDEXPRG'] = 3
        out = derive_variables(d)
        self.assertEqual(out.phenotype.iloc[:4].tolist(),
                         ['neither_positive', 'hba1c_only', 'fpg_only', 'both_positive'])
        self.assertTrue(pd.isna(out.loc[4, 'fpg_positive']))
        self.assertTrue(out.loc[5, 'fpg_available'])
        self.assertFalse(out.loc[5, 'fasting_weight_positive'])
        self.assertFalse(out.loc[5, 'cohort_b'])
        self.assertFalse(out.loc[6, 'cohort_a'])
        self.assertFalse(out.loc[7, 'cohort_a'])
        self.assertTrue(out.loc[8, 'cohort_a'])
        self.assertTrue(out.loc[8, 'pregnancy_uncertain'])
        self.assertEqual(out.loc[0, 'usual_care'], 'yes')
        self.assertEqual(out.loc[0, 'visit_frequency'], 'none')

    def test_response_reasons_are_not_combined(self):
        d = example_records()
        d.loc[:4, 'DIQ010'] = [np.nan, np.nan, 7, 9, 3]
        d.loc[0, 'diq_record_present'] = False
        out = derive_variables(d)
        self.assertEqual(out.diabetes_history.iloc[:5].tolist(),
                         ['not_in_diq', 'item_missing', 'refused', 'dont_know', 'borderline'])
        self.assertFalse(out.cohort_a.iloc[:5].any())

    def test_bridge_only_applies_to_old_cycle(self):
        d = example_records()
        d.loc[0, ['cycle', 'LBXGLU']] = ['2013-2014', 124.]
        d.loc[1, 'LBXGLU'] = 124.
        out = derive_variables(d)
        self.assertAlmostEqual(out.loc[0, 'fpg_harmonized'], 126.3412)
        self.assertEqual(out.loc[1, 'fpg_harmonized'], 124.)
        self.assertTrue(out.loc[0, 'fpg_classification_changed'])
        self.assertFalse(out.loc[1, 'fpg_classification_changed'])
        self.assertEqual(out.loc[0, 'LBXGLU'], 124.)


class ActualDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.master, cls.inventory, _, cls.joins, _ = build_master(ROOT)

    def test_input_integrity_and_dictionary_coverage(self):
        self.assertEqual(len(self.inventory), 24)
        self.assertEqual(len(self.master), 29400)
        self.assertFalse(self.master.duplicated(['cycle', 'SEQN']).any())
        self.assertTrue(all(r['before_n'] == r['after_n'] and r['unmatched_source_n'] == 0
                            for r in self.joins))
        dictionary = pd.read_csv(ROOT / 'data_cleaning/variable_dictionary.csv')
        self.assertTrue(dictionary.variable.is_unique)
        self.assertEqual(set(dictionary.variable), set(self.master.columns))

    def test_weight_checks_and_fasting_exclusions(self):
        expected_mec_zero = [362, 427, 550]
        expected_infants = [405, 396, 357]
        expected_glucose_zero = [245, 229, 180]
        for i, (_, d) in enumerate(self.master.groupby('cycle', sort=True)):
            self.assertEqual(int(d.WTMEC2YR.eq(0).sum()), expected_mec_zero[i])
            self.assertEqual(int(d.RIDAGEYR.eq(0).sum()), expected_infants[i])
            self.assertTrue(d.WTMEC2YR.eq(0).eq(d.RIDSTATR.eq(1)).all())
            self.assertEqual(int((d.LBXGLU.notna() & d.WTSAF2YR.eq(0)).sum()), expected_glucose_zero[i])
            self.assertFalse((d.cohort_b & ~d.fasting_weight_positive).any())
            # HUQ051=0 is a valid category, not a refusal or missing result.
            self.assertTrue(d.loc[d.HUQ051.eq(0), 'visit_frequency'].eq('none').all())

    def test_flow_and_phenotypes_conserve_people(self):
        flow = build_flow(self.master)
        for _, rows in flow.groupby('cycle'):
            self.assertTrue((rows.previous_n - rows.excluded_n).eq(rows.retained_n).all())
            self.assertEqual(rows.previous_n.iloc[1:].tolist(), rows.retained_n.iloc[:-1].tolist())
        m = self.master
        self.assertEqual(int(m.cohort_a.sum()), 12724)
        self.assertEqual(int(m.cohort_b.sum()), 5658)
        self.assertEqual(int(m.phenotype.ne('not_in_cohort_b').sum()), int(m.cohort_b.sum()))
        self.assertFalse((m.cohort_b & ~m.cohort_a).any())
        self.assertEqual(int(flow.loc[(flow.cycle == 'total') & (flow.stage == '06_cohort_a'), 'retained_n'].iloc[0]), int(m.cohort_a.sum()))
        self.assertEqual(int(flow.loc[(flow.cycle == 'total') & (flow.stage == '08_cohort_b'), 'retained_n'].iloc[0]), int(m.cohort_b.sum()))

    def test_fasting_sources_agree_when_both_observed(self):
        d = self.master.loc[self.master.cycle == '2013-2014']
        for field in ['PHAFSTHR', 'PHAFSTMN']:
            both = d[field].notna() & d[field + '_glu'].notna()
            self.assertEqual(int(both.sum()), 3291)
            self.assertTrue(d.loc[both, field].eq(d.loc[both, field + '_glu']).all())


if __name__ == '__main__':
    unittest.main()
