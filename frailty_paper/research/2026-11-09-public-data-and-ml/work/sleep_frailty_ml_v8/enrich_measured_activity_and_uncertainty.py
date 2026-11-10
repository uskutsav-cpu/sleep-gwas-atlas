"""Recover real allele activities and measured uncertainty, preserving old fits."""
import collections
from v8_common import *


def main():
    assert_previous_intact()
    with np.load(DATA/'measured_training_arrays.npz') as z:
        d = {k: z[k] for k in z.files}
    nc = len(d['cell_names']); lookup = np.full((int(d['bank_index'].max())+1, nc), -1, dtype='int32')
    lookup[d['bank_index'], d['cell']] = np.arange(len(d['effect']))
    cells = {s: i for i, s in enumerate(d['cell_names'])}
    d['measurement_standard_error'] = np.full(len(d['effect']), np.nan, dtype='float32')
    d['standard_error_kind'] = np.zeros(len(d['effect']), dtype='int8')
    # kind 1: empirical paired-mean SE across three published assay replicates.
    # kind 2: author's original GLM Skew_SE. These are explicitly distinct.
    before_activity = np.isfinite(d['activity_reference']).sum()
    oldmeta = {}
    for row in csv_rows(V6/'sequence_bank.csv.gz'):
        if row['sequence_source'] == 'siraj':
            oldmeta[int(row['bank_id'])] = int(row['original_row_index'])
    siraj_banks = {}; kosicki_banks = {}
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        bank = int(r['bank_id'])
        if r['context_feature_source'] == 'kosicki':
            kosicki_banks[int(r['context_feature_index'])] = bank
        elif r['context_feature_source'] == 'v6' and int(r['context_feature_index']) in oldmeta:
            siraj_banks[oldmeta[int(r['context_feature_index'])]] = bank
    needed = {}; names = set()
    reconstructed = V6/'outputs/kosicki_reconstructed_continuous_effects.csv'
    for r in csv_rows(reconstructed):
        bank = kosicki_banks.get(int(r['design_index']))
        if bank is None:
            continue
        i = lookup[bank, cells['induced_excitatory_neuron']]
        if i >= 0:
            needed[i] = (r['reference_name'], r['name'])
            names.update(needed[i])
    original = V6/'reserved_outcome_source/published_replicate_readouts.csv'
    prior = json.loads((V6/'outputs/kosicki_continuous_recovery_audit.json').read_text())
    assert sha(original) == prior['replicate_export_sha256']
    measurements = collections.defaultdict(dict)
    for r in csv_rows(original):
        if r['name'] in names:
            rep = int(r['rep']); assert rep not in measurements[r['name']]
            value = float(r['log2.ratio']); assert np.isfinite(value)
            measurements[r['name']][rep] = value
    replicate_log2 = []; replicate_ids = []
    for i, (ref, alt) in needed.items():
        assert set(measurements[ref]) == set(measurements[alt]) == {1, 2, 3}
        a = np.array([measurements[ref][k] for k in [1, 2, 3]])
        b = np.array([measurements[alt][k] for k in [1, 2, 3]])
        delta = b-a; mean = delta.mean()
        assert abs(mean-d['signed_effect'][i]) < 2e-6
        assert abs(abs(mean)-d['effect'][i]) < 2e-6
        d['activity_reference'][i] = a.mean(); d['activity_alternate'][i] = b.mean()
        d['measurement_standard_error'][i] = delta.std(ddof=1)/np.sqrt(3)
        d['standard_error_kind'][i] = 1
        if d['split'][i] != 2:
            replicate_log2.append(delta); replicate_ids.append(i)
    manifest = ROOT/'work/sleep_frailty_ml_v3/siraj_manifest_local.csv'
    mapping = {r['variant']: siraj_banks[i] for i, r in enumerate(csv_rows(manifest)) if i in siraj_banks}
    cr = json.loads((V6/'sources/siraj_additional_cell_receipt.json').read_text())
    br = json.loads((V6/'sources/siraj_pooled_K562_receipt.json').read_text())
    paths = [(BASE/'siraj2026_neural/OL42_SKNSH_emVAR_glm_paired_20221220.out.gz', 'neuronal_tumor')]
    paths += [(Path(r['file']), 'hepatocyte_tumor' if '_HepG2_' in r['file'] else 'colon_tumor') for r in cr['files']]
    paths += [(Path(br['file']), 'hematopoietic')]
    sources = []; duplicate = collections.Counter()
    for path, cell in paths:
        sources.append({'path': str(path), 'sha256': sha(path)}); seen = set()
        with gzip.open(path, 'rt') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if not row['comb'].endswith('_center_fwd_ref'):
                    continue
                bank = mapping.get('chr'+row['SNP'])
                if bank is None:
                    continue
                i = lookup[bank, cells[cell]]
                if i < 0:
                    continue
                if i in seen:
                    d['measurement_standard_error'][i] = np.nan; d['standard_error_kind'][i] = 0
                    duplicate[cell] += 1; continue
                seen.add(i)
                assert abs(abs(float(row['Log2Skew']))-d['effect'][i]) < 1e-5
                try:
                    se = float(row['Skew_SE'])
                except ValueError:
                    continue
                if np.isfinite(se) and se > 0:
                    d['measurement_standard_error'][i] = se; d['standard_error_kind'][i] = 2
    target = DATA/'measured_training_arrays_with_activity_and_SE.npz'
    np.savez_compressed(target, **d)
    repl = np.array(replicate_log2)
    np.savez_compressed(DATA/'development_Kosicki_actual_log2_replicate_differences.npz', observation_index=np.array(replicate_ids), log2_difference=repl)
    write_json(OUT/'measured_activity_uncertainty_enrichment_audit.json', {
        'original_training_arrays_sha256': sha(DATA/'measured_training_arrays.npz'),
        'new_arrays_sha256': sha(target), 'original_arrays_modified': False,
        'original_core_magnitude_targets_or_splits_changed': False,
        'before_paired_activity_rows': int(before_activity),
        'after_paired_activity_rows': int(np.isfinite(d['activity_reference']).sum()),
        'reconstructed_Kosicki_allele_activity_rows': len(needed),
        'signed_mean_difference_matches_original_target': True,
        'empirical_three_replicate_paired_mean_SE_rows': int(np.sum(d['standard_error_kind'] == 1)),
        'original_author_GLM_SE_rows': int(np.sum(d['standard_error_kind'] == 2)),
        'missing_SE_rows': int(np.sum(d['standard_error_kind'] == 0)),
        'SE_definitions': {'1': 'Sample SD of three paired ALT-minus-REF log2 RNA/DNA differences divided by sqrt(3). An empirical estimate assuming replicate independence, not a recovered source GLM SE or invented p value.', '2': 'Original author-provided GLM Skew_SE, exact mapped physical DNA and cell. No sign convention is inferred from the SE.'},
        'development_log2_replicate_rows_excluding_calibration': len(repl),
        'development_log2_replicate_Pearson_matrix': np.corrcoef(repl.T).tolist(),
        'source_replicate_export_sha256': sha(original), 'source_GLM_files': sources,
        'duplicate_source_rows_rejected': dict(duplicate),
        'fresh_test_outcomes_or_new_model_performance_used': False,
        'scope': 'Measured-data enrichment for the next development experiment. Repeatability is not a formal prediction ceiling. Original V8 fits keep their original inputs and hashes.'
    })
    assert_previous_intact()
    print('MEASURED ALLELE ACTIVITY AND UNCERTAINTY RECOVERED', len(needed), int(np.sum(d['standard_error_kind'] > 0)), flush=True)


if __name__ == '__main__':
    main()
