"""Identity/partition inventory for additional activity in training-only regions.

This is a new candidate scope, not a change to the original 3,345-row package.
All source members participate before measurement availability. This stage
rechecks spatial roles; complete new shared-core closure/audits remain required.
No row is declared training-eligible here.
"""
import bisect
import collections
import csv
import datetime
import gzip
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'work/sleep_frailty_ml_resumable'))
from readonly_native_inputs import sha
sys.path.append(str(ROOT / 'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver

SOURCE = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
DEST = SOURCE / 'training_region_extension'
CHROMS = {'chr' + str(i): i for i in range(1, 23)} | {'chrX': 23, 'chrY': 24}
FACTOR = 10**9


def rows(path):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        yield from csv.DictReader(handle)


def key(chrom, position):
    return CHROMS[chrom] * FACTOR + int(position)


def reduce_coordinates(keys, roles, studies):
    order = np.argsort(keys, kind='stable')
    keys, roles, studies = keys[order], roles[order], studies[order]
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
    return keys[starts], np.bitwise_or.reduceat(roles, starts), np.bitwise_or.reduceat(studies, starts)


def main():
    receipt_path = DEST / 'spatial_training_region_inventory_receipt.json'
    assert not receipt_path.exists()
    DEST.mkdir(exist_ok=True)
    spatial = json.loads((SOURCE / 'additional_activity_spatial_and_reference_receipt.json').read_text())
    closed = json.loads((SOURCE / 'additional_activity_closed_core_receipt.json').read_text())
    package = json.loads((SOURCE / 'audited_additional_neuron_activity_package_receipt.json').read_text())
    components = SOURCE / 'all_neuronal_source_identity_components.csv.gz'
    candidates = SOURCE / 'additional_activity_spatial_and_reference_screen.csv.gz'
    assert sha(components) == closed['all_source_components_sha256']
    assert sha(candidates) == spatial['screen_table_sha256']
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    native_plan = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/finite_stable_native_plan.json').read_text())
    assert sha(data) == native_plan['data_sha256']
    with np.load(data, allow_pickle=False) as archive:
        d = {name: archive[name] for name in ['bank_index', 'split', 'study', 'group', 'study_names']}
    n = int(d['bank_index'].max()) + 1
    masks = np.zeros(n, dtype='uint8')
    np.bitwise_or.at(masks, d['bank_index'], (1 << d['split']).astype('uint8'))
    used = np.flatnonzero(masks)
    assert np.isin(masks[used], [1, 2, 4]).all()
    group_studies = np.zeros(int(d['group'].max()) + 1, dtype='uint16')
    np.bitwise_or.at(group_studies, d['group'], (1 << d['study']).astype('uint16'))
    bank_studies = np.zeros(n, dtype='uint16')
    np.bitwise_or.at(bank_studies, d['bank_index'], group_studies[d['group']])
    coordinates = np.zeros(n, dtype='int64')
    started = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/development/frozen_reader_none_42_started.json').read_text())
    current_sources = started['physical_source_receipt']['sources']
    for source in current_sources:
        path = Path(source['path']); assert sha(path) == source['sha256']
        for row in rows(path):
            bank = int(row['bank_id'])
            if bank < n and masks[bank]:
                assert coordinates[bank] == 0
                coordinates[bank] = key(row['chromosome'], int(row['position_hg38']) - 1)
    assert np.all(coordinates[used] > 0)
    current, current_roles, current_studies = reduce_coordinates(coordinates[used], masks[used], bank_studies[used])
    assert np.isin(current_roles, [1, 2, 4]).all()
    del coordinates, masks, bank_studies, group_studies
    # Keep all former independent test designs protected, even if they happen
    # to coincide with an otherwise recognized current training coordinate.
    forced = {'matoba_full_verified_physical_design.csv',
              'Tardaguila_full_reference_verified_design.csv', 'full_reference_and_shared_core_roles.csv'}
    key_parts, role_parts, study_parts, history_profile = [], [], [], []
    for entry in spatial['historical_sources']:
        path = Path(entry['path']); assert sha(path) == entry['sha256']
        values = np.fromiter((key(row['chromosome'], int(row['position_hg38']) - 1)
                              for row in rows(path) if row['chromosome'] in CHROMS), dtype='int64')
        ix = np.searchsorted(current, values)
        known = ix < len(current)
        known[known] &= current[ix[known]] == values[known]
        role = np.full(len(values), 8, dtype='uint8')  # Unassigned historical regions stay protected.
        studies = np.zeros(len(values), dtype='uint16')
        role[known], studies[known] = current_roles[ix[known]], current_studies[ix[known]]
        if path.name in forced:
            role |= 8
        key_parts.append(values); role_parts.append(role); study_parts.append(studies)
        history_profile.append({**entry, 'canonical_positions': len(values),
                                'forced_previous_test_design_exclusion': path.name in forced})
    history, history_roles, history_studies = reduce_coordinates(
        np.concatenate(key_parts), np.concatenate(role_parts), np.concatenate(study_parts))
    del key_parts, role_parts, study_parts, current, current_roles, current_studies
    reserve = ROOT / 'outputs/sleep_frailty_ml_channels/source_candidates/GSE87711_design.txt.gz'
    assert sha(reserve) == spatial['reserved_original_design_sha256']
    chain = ROOT / 'work/hg19ToHg38.over.chain.gz'
    assert sha(chain) == spatial['chain_sha256'][chain.name]
    forward = LiftOver(str(chain))
    with gzip.open(reserve, 'rt') as handle:
        variants = {identifier.split('_')[0] for identifier, sequence in csv.reader(handle, delimiter='\t')}
    reserved = set(); unmapped = 0
    for variant in variants:
        chrom, position, ref, alt = variant.split(':')
        mapped = forward.convert_coordinate('chr' + chrom.removeprefix('chr'), int(position) - 1)
        if not mapped:
            unmapped += 1
        for c, p, strand, score in mapped:
            if c in CHROMS:
                reserved.add(key(c, int(p)))
    reserved = np.array(sorted(reserved), dtype='int64')
    size = spatial['source_rows_before_measurement_filter']
    component_roles = np.zeros(size, dtype='uint8')
    component_studies = np.zeros(size, dtype='uint16')
    component_reserved = np.zeros(size, dtype=bool)
    selected_fields = ['name', 'source_component', 'passes_spatial_and_reference_screen',
                       'screen_exclusions', 'full_actual_reference_match']
    screened = [{name: row[name] for name in selected_fields} for row in rows(candidates)]
    needed = {row['name'] for row in screened}
    roots, source_count = {}, 0
    for row in rows(components):
        source_count += 1
        root = int(row['component_after_shared_145bp_closure'])
        if row['name'] in needed:
            roots[row['name']] = root
        if not row['mapped_GRCh38_interval']:
            continue
        chrom, begin, end = row['mapped_GRCh38_interval'].split(':')
        a, b = int(begin), int(end)
        lo, hi = key(chrom, max(0, a - 2001)), key(chrom, min(FACTOR - 1, b + 2001))
        left, right = np.searchsorted(history, [lo, hi], side='left')
        right = np.searchsorted(history, hi, side='right')
        if right > left:
            component_roles[root] |= np.bitwise_or.reduce(history_roles[left:right])
            component_studies[root] |= np.bitwise_or.reduce(history_studies[left:right])
        lo, hi = key(chrom, max(0, a - 2146)), key(chrom, min(FACTOR - 1, b + 2146))
        left = np.searchsorted(reserved, lo, side='left')
        if left < len(reserved) and reserved[left] <= hi:
            component_reserved[root] = True
    assert source_count == size and set(roots) == needed
    records, counts, fold_counts = [], collections.Counter(), collections.Counter()
    for row in screened:
        root = roots[row['name']]
        role, studies = int(component_roles[root]), int(component_studies[root])
        reasons = set(row['screen_exclusions'].split('|')) - {''}
        fixed_reasons = reasons - {'whole_source_component_near_any_historical_region',
                                   'whole_source_component_near_entire_reserved_study'}
        if 'whole_source_component_near_any_historical_region' in reasons:
            assert role != 0
        if 'whole_source_component_near_entire_reserved_study' in reasons:
            assert component_reserved[root]
        passes = not fixed_reasons and role in [0, 1] and not component_reserved[root]
        kind = ('previous_novel_region_scope' if role == 0 else 'additional_training_only_region_scope') if passes else 'excluded'
        counts[kind] += 1
        excluded_sources = [name for i, name in enumerate(d['study_names']) if (studies & (1 << i)) or name == 'Salomon']
        if passes:
            for name in d['study_names']:
                if name != 'Kircher2019' and name not in excluded_sources:
                    fold_counts[str(name)] += 1
        records.append({'source_construct_name': row['name'], 'existing_source_component': root,
                        'historical_partition_mask': role, 'touches_reserved_study_envelope': bool(component_reserved[root]),
                        'passes_only_this_new_spatial_role_stage': passes, 'candidate_scope': kind,
                        'excluded_whole_source_roles_at_this_stage': '|'.join(map(str, excluded_sources)),
                        'training_eligible': False,
                        'remaining': 'Renew full exact shared-DNA closure and all historical/reserved DNA checks for the expanded scope, duplicates, independent replay, then fix a separate activity experiment.'})
    packaged = SOURCE / 'audited_additional_neuron_activity_examples.csv.gz'
    assert sha(packaged) == package['artifact_sha256']
    old_names = {row['source_construct_name'] for row in rows(packaged)}
    new_names = {row['source_construct_name'] for row in records if row['passes_only_this_new_spatial_role_stage']}
    assert old_names <= new_names
    target = DEST / 'spatial_training_region_candidates.csv.gz'
    with gzip.open(target, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    result = {'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'spatial_and_partition_inventory_complete': True, 'counts': dict(counts),
              'all_source_rows_participated_before_activity_availability': source_count,
              'source_components_are_cached_prior_joins_not_full_new_scope_core_closure': True,
              'historical_sources': history_profile, 'historical_unique_positions': len(history),
              'entire_reserved_source_variants': len(variants), 'reserved_variants_without_position_mapping': unmapped,
              'all_previous_3345_packaged_names_retained_at_this_spatial_stage': True,
              'upper_bound_rows_available_when_other_sources_omitted_before_new_core_audit': dict(fold_counts),
              'entire_new_activity_source_always_omitted_for_Salomon_test': True,
              'policy': 'Relax only the previous blanket training-region exclusion: allow recognized current training-only coordinates. Validation, calibration, unassigned historical regions, previous test-design sources and the entire reserved-study envelope remain protected. Whole original groups carry all source-membership bits into future source-omission exclusions.',
              'original_data_targets_partitions_models_or_3345_package_changed': False,
              'new_training_dataset_ready': False, 'new_models_trained': False, 'reserved_outcomes_read': False,
              'genome_reference_criteria_or_target_magnitude_filters_changed': False,
              'remaining': 'The new spatial candidates are NOT training-eligible. Rebuild complete shared-145bp components for the expanded candidate scope and propagate partition/source protections; scan all original and reserved DNA, handle duplicates and independently replay. Freeze an activity-data experiment only after those checks.',
              'inputs_sha256': {str(path): sha(path) for path in [data, components, candidates, reserve, chain, packaged]},
              'candidate_table_sha256': sha(target), 'code_sha256': sha(Path(__file__))}
    assert sha(data) == native_plan['data_sha256'] and sha(packaged) == package['artifact_sha256']
    receipt_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['counts', 'all_source_rows_participated_before_activity_availability',
                    'upper_bound_rows_available_when_other_sources_omitted_before_new_core_audit', 'new_training_dataset_ready', 'new_models_trained']}), flush=True)


if __name__ == '__main__':
    main()
