import csv
import json
import itertools
import unittest
from pathlib import Path
import numpy as np
import predict

HERE = Path(__file__).resolve().parent


class PredictorTests(unittest.TestCase):
    def test_all_previously_saved_group_holdouts(self):
        names, x, y = predict.load_matrix(HERE / 'source_matrix.tsv')
        with (HERE / 'heldout_reference.tsv').open(newline='') as f:
            expected = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(expected), 36)
        for i, name in enumerate(names):
            actual = predict.predict(names, x, y, x[i], predict.group(name))
            for row in (r for r in expected if r['sleep_trait'] == name):
                result = actual[row['target']]
                for key in ('bounded_ridge', 'training_mean', 'knn3', 'lambda'):
                    self.assertAlmostEqual(result[key], float(row[key]), places=10, msg=(name, row['target'], key))

    def test_clinical_insomnia_common_marker_example(self):
        names, x, y = predict.load_matrix(HERE / 'source_matrix.tsv')
        features = json.loads((HERE / 'mvp_common_support_features.json').read_text())
        supplied = np.array([features[k] for k in predict.FEATURES])
        actual = predict.predict(names, x, y, supplied, 'insomnia')
        with (HERE / 'mvp_common_support_reference.tsv').open(newline='') as f:
            expected = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual({r['target'] for r in expected}, set(predict.TARGETS))
        for row in expected:
            result = actual[row['target']]
            self.assertAlmostEqual(result['bounded_ridge'], float(row['mvp_feature_bounded_prediction']), places=10)
            self.assertAlmostEqual(result['training_mean'], float(row['mean_prediction']), places=10)
            self.assertAlmostEqual(result['knn3'], float(row['knn3_prediction']), places=10)
            self.assertEqual(result['lambda'], float(row['lambda']))

    def test_input_box_matches_all_64_corners(self):
        names, x, y = predict.load_matrix(HERE / 'source_matrix.tsv')
        values = json.loads((HERE / 'mvp_common_support_features.json').read_text())
        errors = json.loads((HERE / 'mvp_common_support_ses.json').read_text())
        feature = np.array([values[k] for k in predict.FEATURES])
        ses = np.array([errors[k] for k in predict.FEATURES])
        result = predict.predict(names, x, y, feature, 'insomnia', ses)
        train = np.array([i for i, n in enumerate(names) if predict.group(n) != 'insomnia'])
        for j, target in enumerate(predict.TARGETS):
            lam = result[target]['lambda']
            corners = [float(predict.ridge(x[train], y[train, j], feature + np.array(signs) * ses, lam))
                       for signs in itertools.product((-1, 1), repeat=6)]
            box = result[target]['one_se_input_box']
            self.assertAlmostEqual(box['low'], min(corners), places=12)
            self.assertAlmostEqual(box['high'], max(corners), places=12)
            self.assertLessEqual(box['low'], result[target]['bounded_ridge'])
            self.assertGreaterEqual(box['high'], result[target]['bounded_ridge'])
        with self.assertRaises(ValueError):
            predict.input_box(x[train], y[train, 0], feature, np.ones(6), 1.)

    def test_related_groups_are_excluded_together(self):
        self.assertEqual({predict.group(x) for x in ('sleepdur', 'shortsleep', 'longsleep')}, {'dashti_2019_duration'})
        self.assertEqual({predict.group(x) for x in ('sleep_efficiency', 'accel_sleep_duration', 'sleep_timing')}, {'jones_accelerometer'})

    def test_feature_range_uses_only_remaining_training_rows(self):
        names = ['sleepdur', 'shortsleep', 'longsleep', 'other', 'third']
        x = np.array([[0.9]*6, [0.8]*6, [0.7]*6, [0.1]*6, [0.2]*6])
        result = predict.feature_range_diagnostic(names, x, np.array([0.5, 0.15, 0.1, 0.2, 0.05, 0.3]), 'dashti_2019_duration')
        self.assertEqual(result['training_rows'], 2)
        flags = [result['features'][key]['outside_training_range'] for key in predict.FEATURES]
        self.assertEqual(flags, [True, False, False, False, True, True])
        self.assertEqual(result['features']['adhd']['training_max'], 0.2)


if __name__ == '__main__':
    unittest.main()
