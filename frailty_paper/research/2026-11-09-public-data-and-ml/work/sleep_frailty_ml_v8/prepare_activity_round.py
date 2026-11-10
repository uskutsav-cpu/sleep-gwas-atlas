"""Development-only amendment before any fresh-study model performance."""
from v8_common import *


def main():
    assert_previous_intact()
    path = OUT/'measured_activity_round_plan.json'
    plan = {
        'amendment': 'After observing weak initial V8 study-to-study transfer, before any Matoba model performance.',
        'original_fits_and_data_preserved': True,
        'input_data': 'measured_training_arrays_with_activity_and_SE.npz',
        'input_data_sha256': sha(DATA/'measured_training_arrays_with_activity_and_SE.npz'),
        'recipes': ['activity', 'activity_precision', 'local_activity_precision'],
        'backbone': 'Real MPRA-LegNet WTC11 weights; tune original mapper and dense representation while retaining frozen convolutional features for the controlled initial comparison.',
        'activity_targets': 'Directly supervise both original measured allele activities, centered within original study/cell context using training rows only. Preserve log2 units. Also supervise their difference separately from original fitted magnitude.',
        'primary_objective': 'Unchanged raw absolute-magnitude MAE for all primary training examples, balanced by study/cell context. No significance label or test-outcome filtering.',
        'precision_use': 'When available, bounded reliability factors from measured SE affect auxiliary signed/activity-difference and local/cell-pair losses only. Main raw-magnitude loss stays unchanged. Unknown SE is explicitly unknown and receives default unmodified training weight, not a claim of certainty.',
        'local_features': 'Actual pretrained final convolutional activations around the known substituted position, averaged over the containing bin and its immediate valid neighbors; preserve both physical allele sequences and source orientation. Local branch starts with zero added input coefficients.',
        'additional_global_features': 'Predicted mean allele activity and global representation norm, to preserve activity intensity removed by per-row normalization.',
        'held_out_development_studies': ['none', 'Kosicki2025', 'IGVF240K', 'Salomon', 'Siraj'],
        'fixed_epochs': 8, 'seed': 42, 'learning_rate': .0003,
        'selection_rule': 'One common recipe and epoch maximizing macro raw-magnitude error skill averaged over same-cell-supported whole-study omissions; report withheld-locus, per-source, local and within-cell-contrast results separately. Do not optimize a different epoch for every held study and call it confirmatory validation.',
        'baseline_requirements': ['Training-cell median for magnitude', 'Zero difference and training-cell median difference for cell comparisons', 'Within-fixed-study-and-cell-contrast correlation to distinguish variant interactions from average cell shifts'],
        'fresh_study_policy': 'Do not read model performance on the 610 source-eligible Matoba effects until recipe, weights, calibration, AVI experiment and predictions are fixed.',
        'calibration_rows_not_used_for_development_selection': True,
        'AlphaGenome_training_targets': 'None. Only designated permissible composite AVI can be an additional input later, under a matched comparison.',
        'scope': 'A finite new development experiment, not a success claim or an alteration of the preserved V7 failure and AVI rejection.'
    }
    if path.exists():
        assert json.loads(path.read_text()) == plan
    else:
        write_json(path, plan)
    print('MEASURED ACTIVITY ROUND DECLARED BEFORE FRESH PERFORMANCE', flush=True)


if __name__ == '__main__':
    main()
