"""Use earlier public checks as development, never as fresh validation again."""
import bisect
from v9_common import *


def numeric(x):
    try:
        value = float(x)
        return value if np.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def main():
    assert_previous_sources()
    assert not (OUT/'model_search_started.json').exists()
    plan = {
        'why': 'All tested V8 neural recipes failed to beat the typical-effect baseline consistently across whole development studies. Relevant cells had insufficient study coverage despite already-observed public data being available.',
        'amendment_before_fresh_model_performance': True,
        'added_sources': ['McAfee2023 original 9,602 measured outcomes in primary neural progenitors and HEK293', 'Retallick-Townsley2026 all source-verified baseline biallelic finite effects, including previously excluded development-overlap regions after whole-component repartitioning'],
        'earlier_V7_validation_roles': 'Explicitly retired to revealed DEVELOPMENT in V9. Preserve their original frozen results; never call their V9 performance independent replication.',
        'new_test': 'The same 610 eligible Matoba baseline effects remain unused for model performance. Preserve the original physical-design eligibility and journal-source proof.',
        'calibration_scope': 'Preserve old calibration loci. If a new bridge joins them to old training/validation, exclude the entire connected component. Newly partitioned rows from earlier revealed studies can calibrate errors but are not a new independent validation.',
        'primary_targets': 'Original absolute fitted log2 effects, no significance filtering or imputed missing outcomes.',
        'cell_context': 'Retain actual cell-line identities, including HNP_primary and HEK293, separately from broader cell groups. Compare DNA models with both cell-group and line-specific median baselines.',
        'future_model_comparison': 'One finite line-aware paired-DNA architecture, frozen versus partially adapted pretrained representations; use measured activity/precision auxiliaries without feeding measured values as inference inputs. Select using development studies only, before a matched AVI comparison and fresh test.',
        'scope': 'Laboratory molecular prediction, not independent human genetic-score validation or a disease-risk model.'
    }
    plan_path = OUT/'revealed_source_amendment_plan.json'
    if plan_path.exists(): assert json.loads(plan_path.read_text()) == plan
    else: write_json(plan_path, plan)
    with np.load(V8/'measured_training_arrays_with_activity_and_SE.npz') as z: old = {k: z[k] for k in z.files}
    nold = len(old['effect']); b0 = int(old['bank_index'].max())+1
    names = {key: old[key].tolist() for key in ['cell_names', 'cell_line_names', 'study_names', 'experiment_names', 'context_names']}
    def code(key, value):
        if value not in names[key]: names[key].append(value)
        return names[key].index(value)
    newbanks = []; added = []; source_exclusions = collections.Counter()
    def bank(row, source):
        a, b = row['reference'], row['alternate']; pos = int(row['changed_position'])
        assert len(a) == len(b) and len(a) >= 150 and set(a+b) <= set('ACGT')
        assert [i for i,(x,y) in enumerate(zip(a,b)) if x != y] == [pos]
        i = b0+len(newbanks)
        newbanks.append({'bank_id': i, 'variant_GRCh38': row['variant_GRCh38'], 'chromosome': row['chromosome'], 'position_hg38': int(row['position_hg38']), 'reference': a, 'alternate': b, 'changed_position': pos, 'construct_strand': row['construct_strand'], 'source': source})
        return i
    def observation(bank_id, source, experiment, cell, line, value, a=None, b=None):
        value = float(value); assert np.isfinite(value) and value >= 0
        added.append({'bank_index': bank_id, 'effect': value, 'cell': code('cell_names', cell), 'line': code('cell_line_names', line), 'study': code('study_names', source), 'experiment': code('experiment_names', source+'|'+experiment+'|'+line), 'context': code('context_names', source+'|'+cell+'|absolute_log2_effect'), 'activity_reference': a if a is not None else np.nan, 'activity_alternate': b if b is not None else np.nan, 'signed_effect': np.nan, 'measurement_standard_error': np.nan, 'standard_error_kind': 0, 'original_V7_observation_index': -1})
    audit = json.loads((V7/'outputs/reserved_design_and_original_source_audit.json').read_text())
    assert sha(V7/'reserved_physical_design.csv') == audit['physical_design_sha256']
    assert sha(V7/'reserved_original_effects.csv.gz') == audit['original_effects_file_sha256']
    mc = list(csv_rows(V7/'reserved_original_effects.csv.gz')); used_design = {int(r['design_index']) for r in mc}; mapping = {}
    for r in csv_rows(V7/'reserved_physical_design.csv'):
        if int(r['design_index']) in used_design: mapping[int(r['design_index'])] = bank(r, 'McAfee2023')
    for r in mc:
        observation(mapping[int(r['design_index'])], 'McAfee2023', 'original_'+r['cell_line'], r['cell_group'], r['cell_line'], r['absolute_log2_effect'])
    assert len(mc) == 9602
    na = json.loads((V7/'outputs/neuroimmune_candidate_source_audit.json').read_text())
    assert sha(V7/'neuroimmune_verified_candidate_effects.csv') == na['verified_effect_file_sha256']
    raw = collections.defaultdict(list)
    for r in csv_rows(V7/'neuroimmune_source_tables/baseline.csv'): raw[r['ID']].append(r)
    baseline_count = 0
    for r in csv_rows(V7/'neuroimmune_physical_design_verified.csv'):
        rows = raw.get(r['source_ID'], []); byname = collections.defaultdict(list)
        for rr in rows: byname[rr['CRS']].append(rr)
        if len(byname) != 2 or any(k not in byname for k in [r['source_ref_construct'], r['source_alt_construct']]):
            source_exclusions['Retallick_missing_or_nonbiallelic_postQC_pair'] += 1; continue
        if any(len({rr.get('emVar.logFC', '') for rr in values}) > 1 for values in byname.values()):
            source_exclusions['Retallick_duplicate_effect_disagreement'] += 1; continue
        rr, aa = byname[r['source_ref_construct']][0], byname[r['source_alt_construct']][0]
        effect_a, effect_b = numeric(rr.get('emVar.logFC')), numeric(aa.get('emVar.logFC'))
        if effect_a is None or effect_b is None or abs(abs(effect_a)-abs(effect_b)) > 1e-8:
            source_exclusions['Retallick_missing_or_inconsistent_fitted_effect'] += 1; continue
        i = bank(r, 'Retallick2026')
        observation(i, 'Retallick2026', 'baseline_NGN2', 'induced_excitatory_neuron', 'Retallick_NGN2_two_donors', abs(effect_a), numeric(rr.get('quant.logFC')), numeric(aa.get('quant.logFC')))
        baseline_count += 1
    assert baseline_count == na['published_effect_measurement_inventory']['baseline']['verified_biallelic_measured_effects_before_overlap'] == 1016
    totalb = b0+len(newbanks)
    used = np.zeros(totalb, dtype=bool); used[old['bank_index']] = True; used[b0:] = True
    parent = np.arange(totalb, dtype='int32')
    def find(i):
        while parent[i] != i: parent[i] = parent[parent[i]]; i = int(parent[i])
        return i
    def union(i, j):
        a, b = find(i), find(j)
        if a != b: parent[max(a,b)] = min(a,b)
    oldgroup = np.full(b0, -1, dtype='int32'); oldrole = np.full(b0, -1, dtype='int8')
    oldgroup[old['bank_index']] = old['group']; oldrole[old['bank_index']] = old['split']
    first = {}
    for i in np.flatnonzero(used[:b0]):
        g = int(oldgroup[i])
        if g in first: union(i, first[g])
        else: first[g] = i
    chrom_names = {}; chromosome = np.zeros(totalb, dtype='int16'); position = np.zeros(totalb, dtype='int64')
    seen_core = {}; core_fresh = {}; fresh_points = collections.defaultdict(list)
    for r in csv_rows(V8/'reserved_matoba_physical_design.csv'):
        fresh_points[r['chromosome']].append(int(r['position_hg38']))
        for s in [r['reference'], r['alternate']]: core_fresh[s] = True; core_fresh[rc(s)] = True
    for c in fresh_points: fresh_points[c].sort()
    fresh_prefix = {s[:16] for s in core_fresh}; bad_nodes = set()
    def inspect_row(r):
        i = int(r['bank_id'])
        if not used[i]: return
        c = r['chromosome']; p = int(r['position_hg38']); chromosome[i] = chrom_names.setdefault(c, len(chrom_names)+1); position[i] = p
        ps = fresh_points[c]; j = bisect.bisect_left(ps, p)
        if any(0 <= k < len(ps) and abs(ps[k]-p) <= 2000 for k in [j-1,j]): bad_nodes.add(i)
        a, b = r['reference'], r['alternate']; at = int(r['changed_position']); start = max(0, min(at-75, len(a)-150))
        for s in [a,b]:
            core = s[start:start+150]; canonical = min(core, rc(core))
            if canonical in seen_core: union(i, seen_core[canonical])
            else: seen_core[canonical] = i
            for k in range(max(0, at-149), min(at, len(s)-150)+1):
                if s[k:k+16] in fresh_prefix and s[k:k+150] in core_fresh: bad_nodes.add(i)
    for r in csv_rows(V7/'sequence_bank.csv.gz'): inspect_row(r)
    for r in newbanks: inspect_row(r)
    ids = np.flatnonzero(used); order = ids[np.lexsort((position[ids], chromosome[ids]))]
    for i,j in zip(order[:-1], order[1:]):
        if chromosome[i] == chromosome[j] and position[j]-position[i] <= 2000: union(i,j)
    roots = {i: find(int(i)) for i in ids}; bad_roots = {roots[i] for i in bad_nodes}
    oldroles = collections.defaultdict(set)
    for i in np.flatnonzero(used[:b0]): oldroles[roots[i]].add(int(oldrole[i]))
    cal_conflicts = {r for r,roles in oldroles.items() if 2 in roles and len(roles) > 1}
    excluded_roots = bad_roots|cal_conflicts
    root_points = collections.defaultdict(list)
    for i in ids:
        if roots[i] not in excluded_roots: root_points[roots[i]].append((int(chromosome[i]), int(position[i])))
    root_roles = {}; hashes = {}
    for root, pts in root_points.items():
        digest = hashlib.sha256(('V9_revealed_sources_20261108|'+'|'.join(f'{c}:{p}' for c,p in sorted(set(pts)))).encode()).hexdigest()
        hashes[root] = digest[:24]
        roles = oldroles.get(root, set())
        if roles: root_roles[root] = max(roles)
        else:
            fold = int(digest[:8],16)%10; root_roles[root] = 2 if fold < 2 else 1 if fold < 4 else 0
    allkeys = ['bank_index','effect','cell','line','study','experiment','context','activity_reference','activity_alternate','signed_effect','measurement_standard_error','standard_error_kind','original_V7_observation_index']
    d = {k: np.r_[old[k], np.array([r[k] for r in added], dtype=old[k].dtype)] for k in allkeys}
    d['original_V8_observation_index'] = np.r_[np.arange(nold, dtype='int32'), np.full(len(added), -1, dtype='int32')]
    d['source_previously_external_in_V7'] = np.r_[np.zeros(nold, dtype='int8'), np.ones(len(added), dtype='int8')]
    keep = np.array([roots[int(i)] not in excluded_roots for i in d['bank_index']])
    for key in list(d): d[key] = d[key][keep]
    ordered_roots = sorted(root_points); gi = {r:i for i,r in enumerate(ordered_roots)}
    d['group'] = np.array([gi[roots[int(i)]] for i in d['bank_index']], dtype='int32')
    d['split'] = np.array([root_roles[roots[int(i)]] for i in d['bank_index']], dtype='int8')
    for key, values in names.items(): d[key] = np.array(values)
    d['locus_names'] = np.array([hashes[r] for r in ordered_roots])
    assert np.isfinite(d['effect']).all() and np.all(d['effect'] >= 0)
    for key, values in [('group',d['group']),('split',d['split'])]:
        mapping_check = {}
        for bank_id, value in zip(d['bank_index'], values):
            if int(bank_id) in mapping_check: assert mapping_check[int(bank_id)] == value
            else: mapping_check[int(bank_id)] = value
    np.savez_compressed(DATA/'training_arrays.npz', **d)
    with gzip.open(DATA/'extra_sequence_bank.csv.gz','wt') as f:
        writer = csv.DictWriter(f, fieldnames=list(newbanks[0])); writer.writeheader(); writer.writerows(newbanks)
    np.savez_compressed(DATA/'partition_bank_metadata.npz', bank_index=ids, chromosome=chromosome[ids], position=position[ids], component_root=np.array([roots[int(i)] for i in ids]), excluded=np.array([roots[int(i)] in excluded_roots for i in ids]), old_role=np.r_[oldrole, np.full(len(newbanks), -1, dtype='int8')][ids])
    inventory = {}
    for s in np.unique(d['study']):
        inventory[str(d['study_names'][s])] = {str(d['cell_names'][c]): {str(role):int(np.sum((d['study'] == s)&(d['cell'] == c)&(d['split'] == role))) for role in [0,1,2]} for c in np.unique(d['cell'][d['study'] == s])}
    report = {
        'observations': len(d['effect']), 'source_inventory': inventory,
        'original_V8_observations_before_repartitioning': nold,
        'added_public_measured_outcomes_before_genomic_exclusion': len(added),
        'McAfee_original_measured_outcomes': len(mc), 'Retallick_original_baseline_measured_outcomes': baseline_count,
        'new_physical_DNA_pairs': len(newbanks), 'old_feature_bank_offset': b0,
        'excluded_fresh_test_touching_components': len(bad_roots),
        'excluded_old_calibration_bridge_conflict_components': len(cal_conflicts),
        'excluded_observations': int((~keep).sum()),
        'previously_revealed_added_sources_are_development_not_new_validation': True,
        'original_V8_calibration_roles_preserved_or_whole_conflict_component_excluded': True,
        'newly_assigned_revealed_source_calibration_is_not_independent_confirmation': True,
        'source_significance_or_activity_hit_filter': False,
        'unverified_signed_conventions_or_missing_SE_imputed': False,
        'fresh_Matoba_design_sha256_unchanged': sha(V8/'reserved_matoba_physical_design.csv'),
        'fresh_Matoba_candidate_effects': 610, 'fresh_model_performance_calculated': False,
        'source_exclusions': dict(source_exclusions),
        'outputs': {p.name:sha(p) for p in [DATA/'training_arrays.npz',DATA/'extra_sequence_bank.csv.gz',DATA/'partition_bank_metadata.npz']},
        'source_hashes': {'V8_enriched_arrays': sha(V8/'measured_training_arrays_with_activity_and_SE.npz'), 'McAfee_original_effects': sha(V7/'reserved_original_effects.csv.gz'), 'Retallick_baseline_source_table': sha(V7/'neuroimmune_source_tables/baseline.csv'), 'Retallick_genome_verified_design': sha(V7/'neuroimmune_physical_design_verified.csv')},
        'scope': 'Augmented PUBLIC experimental development data with unchanged raw effect units, known cell lines and complete source-role disclosure. Not independent human genetic validation, new participants, causation or successful prediction.'
    }
    write_json(OUT/'revealed_source_augmentation_audit.json', report)
    write_json(V8/'outputs/external_evaluation_preservation_receipt.json', {'V8_neural_development_transfer_failed_to_beat_simple_baseline': True, 'V8_final_fits_preserved': True, 'fresh_test_performance_not_calculated': True, 'next_version': str(DATA), 'reason': 'Improve relevant cell/study coverage by legitimately using already revealed public validation experiments as development, before consuming the untouched Matoba test. No original result or model is overwritten.'})
    assert_previous_sources()
    print('REVEALED PUBLIC SOURCES ADDED WITH FRESH TEST PROTECTED',len(d['effect']), 'observations; excluded',int((~keep).sum()),flush=True)


if __name__ == '__main__':
    main()
