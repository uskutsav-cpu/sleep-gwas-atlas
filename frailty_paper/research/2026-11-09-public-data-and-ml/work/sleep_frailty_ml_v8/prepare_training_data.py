"""Measured targets, original DNA and strict new-test exclusions for V8."""
import collections
from v8_common import *
from legnet_runtime import crop_pair


def main():
    assert not (OUT/'full_model_search_started.json').exists()
    assert_previous_intact()
    audit = json.loads((OUT/'fresh_matoba_physical_design_audit.json').read_text())
    assert sha(DATA/'reserved_matoba_physical_design.csv') == audit['outputs_sha256']['reserved_matoba_physical_design.csv']
    bad = set(json.loads((OUT/'V7_development_groups_excluded_for_fresh_test.json').read_text()))
    with np.load(V7/'training_arrays.npz') as z:
        old = {k: z[k] for k in z.files if k not in ['local_pairs', 'cross_cell_pairs', 'position51']}
    old_bank_count = int(old['bank_index'].max())+1
    bad_codes = {i for i, name in enumerate(old['locus_names']) if name in bad}
    eligible = (old['unit'] == 0)&~np.isin(old['group'], sorted(bad_codes))
    ids = np.flatnonzero(eligible)
    keys = ['bank_index', 'effect', 'split', 'cell', 'study', 'experiment', 'group', 'context', 'line']
    d = {k: old[k][ids] for k in keys}
    for k in ['cell_names', 'cell_line_names', 'study_names', 'experiment_names', 'context_names', 'locus_names']:
        d[k] = old[k]
    d['original_V7_observation_index'] = ids
    n = len(ids); nc = len(d['cell_names'])
    codes = d['bank_index'].astype('int64')*nc+d['cell']
    assert len(np.unique(codes)) == n, 'Repeated bank/cell target requires explicit treatment'
    lookup = np.full((old_bank_count, nc), -1, dtype='int32')
    lookup[d['bank_index'], d['cell']] = np.arange(n)
    d['activity_reference'] = np.full(n, np.nan, dtype='float32')
    d['activity_alternate'] = np.full(n, np.nan, dtype='float32')
    d['signed_effect'] = np.full(n, np.nan, dtype='float32')
    ref = np.lib.format.open_memmap(DATA/'tokens_reference150.npy', mode='w+', dtype='uint8', shape=(old_bank_count, 150))
    alt = np.lib.format.open_memmap(DATA/'tokens_alternate150.npy', mode='w+', dtype='uint8', shape=(old_bank_count, 150))
    lut = np.full(256, 255, dtype='uint8'); lut[[ord(c) for c in 'AGCT']] = [0, 1, 2, 3]
    oldmeta = {}
    for row in csv_rows(V6/'sequence_bank.csv.gz'):
        if row['sequence_source'] in ['siraj', 'salomon']:
            oldmeta[int(row['bank_id'])] = (row['sequence_source'], int(row['original_row_index']))
    maps = collections.defaultdict(dict)
    kosicki_banks = {}
    for i, row in enumerate(csv_rows(V7/'sequence_bank.csv.gz')):
        assert int(row['bank_id']) == i
        a, b = crop_pair(row)
        ref[i] = lut[np.frombuffer(a.encode(), dtype='uint8')]
        alt[i] = lut[np.frombuffer(b.encode(), dtype='uint8')]
        if row['context_feature_source'] == 'v6' and int(row['context_feature_index']) in oldmeta:
            source, index = oldmeta[int(row['context_feature_index'])]
            maps[source][index] = i
        elif row['context_feature_source'] == 'kosicki':
            kosicki_banks[int(row['context_feature_index'])] = i
    ref.flush(); alt.flush()
    assert i+1 == old_bank_count and np.max(ref) < 4 and np.max(alt) < 4
    assert np.all(np.sum(ref != alt, axis=1) == 1)
    del oldmeta
    cell_codes = {s: i for i, s in enumerate(d['cell_names'])}
    sources = []
    def bind(path):
        sources.append({'path': str(path), 'sha256': sha(path)})
    # Direct measured ALT-minus-REF log2 ratios, already reconstructed and
    # physically verified independently in V6; do not infer a sign from abs(y).
    kp = V6/'outputs/kosicki_reconstructed_continuous_effects.csv'; bind(kp)
    for row in csv_rows(kp):
        bank = kosicki_banks.get(int(row['design_index']))
        if bank is None:
            continue
        j = lookup[bank, cell_codes['induced_excitatory_neuron']]
        if j >= 0:
            value = float(row['mean_log2_ratio_difference'])
            assert abs(abs(value)-d['effect'][j]) < 2e-6
            d['signed_effect'][j] = value
    # Per-allele ACTIVITY is a distinct auxiliary target; it is never silently
    # substituted for the author's fitted allelic contrast coefficient.
    manifest = {}
    for source in ['salomon', 'siraj']:
        p = ROOT/'work/sleep_frailty_ml_v3'/(source+'_manifest_local.csv'); bind(p)
        manifest[source] = {row['variant']: maps[source][i] for i, row in enumerate(csv_rows(p)) if i in maps[source]}
    sp = ROOT/'work/salomon_neuronal_variant_readout.tsv.gz'; bind(sp)
    with gzip.open(sp, 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            bank = manifest['salomon'].get(row['chr_pos_ref_alt'])
            if bank is None:
                continue
            j = lookup[bank, cell_codes['induced_excitatory_neuron']]
            if j < 0:
                continue
            a = float(row['NGN2_ref_region_logFC']); b = float(row['NGN2_alt_region_logFC'])
            if np.isfinite([a, b]).all():
                d['activity_reference'][j] = a; d['activity_alternate'][j] = b
    cr = json.loads((V6/'sources/siraj_additional_cell_receipt.json').read_text())
    br = json.loads((V6/'sources/siraj_pooled_K562_receipt.json').read_text())
    paths = [(BASE/'siraj2026_neural/OL42_SKNSH_emVAR_glm_paired_20221220.out.gz', 'neuronal_tumor')]
    paths += [(Path(r['file']), 'hepatocyte_tumor' if '_HepG2_' in r['file'] else 'colon_tumor') for r in cr['files']]
    paths += [(Path(br['file']), 'hematopoietic')]
    duplicate_auxiliary = collections.Counter()
    for path, cell in paths:
        bind(path); seen = set()
        with gzip.open(path, 'rt') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if not row['comb'].endswith('_center_fwd_ref'):
                    continue
                bank = manifest['siraj'].get('chr'+row['SNP'])
                if bank is None:
                    continue
                j = lookup[bank, cell_codes[cell]]
                if j < 0:
                    continue
                if j in seen:
                    d['activity_reference'][j] = np.nan; d['activity_alternate'][j] = np.nan
                    duplicate_auxiliary[cell] += 1; continue
                seen.add(j)
                a, b = float(row['A_log2FC']), float(row['B_log2FC'])
                if np.isfinite([a, b]).all():
                    d['activity_reference'][j] = a; d['activity_alternate'][j] = b
    assert np.isfinite(d['effect']).all() and np.all(d['effect'] >= 0)
    np.savez_compressed(DATA/'measured_training_arrays.npz', **d)
    write_json(OUT/'measured_training_data_audit.json', {
        'observations': n, 'primary_measured_log2_magnitude_only': True,
        'excluded_auxiliary_other_unit_observations': int(np.sum(old['unit'] != 0)),
        'excluded_whole_groups_for_fresh_test': len(bad),
        'excluded_primary_observations_for_fresh_test': int(np.sum((old['unit'] == 0)&~eligible)),
        'independently_signed_Kosicki_targets': int(np.isfinite(d['signed_effect']).sum()),
        'paired_allele_activity_auxiliary_targets': int(np.isfinite(d['activity_reference']).sum()),
        'duplicate_source_auxiliary_rows_rejected': dict(duplicate_auxiliary),
        'signed_conventions_for_other_sources_not_assumed': True,
        'magnitude_is_original_reported_contrast_not_reconstructed_from_activity': True,
        'crop_rule': '150 actual bases, centered on SNP where possible, clipped within original physical construct. Preserve source assay orientation. No new genomic/flanking bases or synthetic targets.',
        'pretrained_channel_order': 'A,G,C,T',
        'native_crop_is_shorter_than_230bp_pretraining_domain': True,
        'split': 'Reuse immutable V7 genomic-component train/validation/calibration roles after removing all new-test-overlapping whole components. Calibrating a model does not establish useful transfer.',
        'inventory': {str(s): {str(c): int(np.sum((d['study'] == i)&(d['cell'] == j))) for j,c in enumerate(d['cell_names']) if np.any((d['study'] == i)&(d['cell'] == j))} for i,s in enumerate(d['study_names'])},
        'source_files': sources, 'V7_targets_sha256': sha(V7/'training_arrays.npz'),
        'output_arrays_sha256': sha(DATA/'measured_training_arrays.npz'),
        'output_token_sha256': {p.name: sha(p) for p in [DATA/'tokens_reference150.npy', DATA/'tokens_alternate150.npy']},
        'fresh_Matoba_effect_performance_read': False,
    })
    assert_previous_intact()
    print('V8 MEASURED DEVELOPMENT DATA READY', n, flush=True)


if __name__ == '__main__':
    main()
