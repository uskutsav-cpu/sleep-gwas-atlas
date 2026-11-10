"""Pure DNA features for final-test and project inputs; no target values."""
from extract_variant_local_features import local_pool
from legnet_runtime import *


def main():
    assert_previous_intact()
    rows = []; bank = []; index = {}
    def add(r, role, identifier, cell, line):
        a, b = r['reference'], r['alternate']; pieces = r['variant_GRCh38'].split(':')
        why = ''
        if len(pieces) != 4 or len(pieces[2]) != 1 or len(pieces[3]) != 1:
            why = 'unsupported_non_SNV'
        elif len(a) != len(b) or len(a) < 150 or set(a+b)-set('ACGT'):
            why = 'invalid_or_insufficient_actual_DNA'
        changes = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        if not why and len(changes) != 1:
            why = 'requires_one_actual_substitution'
        bank_id = -1; changed = -1
        if not why:
            p = changes[0]; start = max(0, min(p-75, len(a)-150))
            a, b = a[start:start+150], b[start:start+150]; changed = p-start
            key = (a, b)
            if key not in index:
                index[key] = len(bank); bank.append(key)
            bank_id = index[key]
        rows.append({'id': identifier, 'role': role, 'variant_GRCh38': r['variant_GRCh38'], 'rsid': r.get('rsid', ''), 'cell_group': cell, 'cell_line': line, 'condition': 'baseline', 'bank_index': bank_id, 'unscorable_reason': why, 'changed_position150': changed, 'source_assay_activity_known': role == 'fresh_test', 'interpretation': 'Source active-element subset' if role == 'fresh_test' else 'Hypothetical genomic reporter context; actual construct activity not established'})
    for r in csv_rows(DATA/'reserved_matoba_physical_design.csv'):
        if r['source_baseline_allelic_coefficient_available'] == 'True':
            add(r, 'fresh_test', 'matoba__'+r['rsid'], 'neural_progenitor', 'HNP_primary')
    for r in csv_rows(V7/'all_project_variant_inputs.csv'):
        add(r, 'project', r['id'], r['cell_group'], r['cell_line'])
    assert sum(r['role'] == 'fresh_test' for r in rows) == 610
    assert all(not r['unscorable_reason'] for r in rows if r['role'] == 'fresh_test')
    write_csv(DATA/'blind_inference_row_metadata.csv', rows)
    model = base_model(); device = 'mps' if torch.backends.mps.is_available() else 'cpu'; model.to(device)
    features = np.empty((len(bank), 2, 2, 128), dtype='float32')
    for i in range(0, len(bank), 256):
        pairs = bank[i:i+256]; ref = [a for a, b in pairs]; alt = [b for a, b in pairs]
        positions = [next(j for j, (x, y) in enumerate(zip(a, b)) if x != y) for a, b in pairs]
        x = encode(ref+alt, device); pos = torch.tensor(positions+positions, device=device)
        with torch.no_grad():
            acts = model.main(model.stem(x)); global_values = acts.mean(-1).cpu().numpy(); local_values = local_pool(acts, pos).cpu().numpy()
        n = len(pairs)
        features[i:i+n, 0, 0] = global_values[:n]; features[i:i+n, 1, 0] = global_values[n:]
        features[i:i+n, 0, 1] = local_values[:n]; features[i:i+n, 1, 1] = local_values[n:]
    np.save(DATA/'blind_inference_pretrained_features.npy', features)
    assert np.isfinite(features).all()
    write_json(OUT/'blind_inference_feature_receipt.json', {
        'fresh_test_rows': 610, 'project_rows': sum(r['role'] == 'project' for r in rows),
        'project_unsupported_rows': sum(bool(r['unscorable_reason']) for r in rows if r['role'] == 'project'),
        'unique_physical_DNA_pairs': len(bank), 'feature_axes': ['physical pair', 'REF/ALT', 'global/local', '128 channels'],
        'input_row_metadata_sha256': sha(DATA/'blind_inference_row_metadata.csv'),
        'feature_sha256': sha(DATA/'blind_inference_pretrained_features.npy'),
        'source_design_sha256': sha(DATA/'reserved_matoba_physical_design.csv'),
        'source_project_inputs_sha256': sha(V7/'all_project_variant_inputs.csv'),
        'pretrained_checkpoint_sha256': sha(MODELS/'author_WTC11_original.safetensors'),
        'fresh_fitted_coefficients_pvalues_or_activity_values_read': False,
        'trained_research_model_test_predictions_generated': False,
        'scope': 'Independent frozen feature construction from actual DNA only, before final predictions. No outcome-based variant selection or claim of measured project activity.'
    })
    print('BLIND PURE-DNA INFERENCE FEATURES PREPARED', len(bank), flush=True)


if __name__ == '__main__':
    main()
