#!/usr/bin/env python3
"""Assemble all 25 Brain6 candidates and 20 regions without changing source analyses."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BOUNDED = ROOT / 'brain6/results/brain6_bounded_manuscript_v1'
ALT = ROOT / 'brain6/results/brain6_alternative_local_validation_v1'
FM = ROOT / 'brain6/results/brain6_exploratory_finemap_coloc_v2'
FUNCTION = ROOT / 'brain6/results/brain6_exploratory_functional_v1'
ENRICH = ROOT / 'brain6/results/brain6_exploratory_enrichment_v3'
CANON = ROOT / 'work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv'
PROJECTION = ROOT / 'brain6/results/confirmatory_source_rescue_20260927/phase5_family_projection_v2/summary.json'
OUT = ROOT / 'brain6/paper/final_package_v1'
EVIDENCE_CLASSES = ('PRIMARY_LAVA_CONFIRMED', 'ALTERNATIVE_LOCAL_VALIDATION_SUPPORTED',
                    'PLACO_PLUS_REPLICATION_SUPPORTED', 'EXPLORATORY_MULTIOMIC_SUPPORTED',
                    'UNSUPPORTED')
FROZEN_CANON_SHA256 = 'ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352'
FROZEN_GLOBAL_SHA256 = '0aa63e5fe446c6828b473069a2c75317131f76cb8f0f1fb19d83cd330b68e28d'
FROZEN_CANDIDATE_SHA256 = 'fffc0e23ab8c111385d78d00168d8dbdbf20b935e7a8fa4b7b6a77e754f3180c'
FROZEN_CANON_CONFIG_SHA256 = 'f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417'
FROZEN_SHARED_RULE_SHA256 = '66a9c923ca23c2349e1a690466d3897312b42891b9c2947715ef6bf52610b1d1'
FROZEN_PLACO_LOCK_SHA256 = 'e11a96bbb77ac6ad3e42426d3abdd927a7cf28c6d1150875d646d20a0d02c8b7'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def source_label(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f'Empty final table: {path.name}')
    fields = fields or list(rows[0])
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=fields, lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def unique(rows: list[dict[str, str]], key: str, expected: int) -> dict[str, dict[str, str]]:
    values = {r[key]: r for r in rows}
    if len(rows) != expected or len(values) != expected:
        raise ValueError(f'{key} must have {expected} unique rows, found {len(rows)}/{len(values)}')
    return values


def numeric(value: str | None) -> float | None:
    if value in (None, '', 'NA', 'N/A'):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def main() -> None:
    if OUT.exists():
        raise FileExistsError('Final package already exists; refuse overwrite')
    source_paths = {
        'bounded_candidates': BOUNDED / 'candidate_evidence_25.tsv',
        'protected_candidates': ROOT / 'brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv',
        'canonical_lava_config': ROOT / 'brain6/config/lava_family_canonical_v3.json',
        'shared_locus_rule': ROOT / 'brain6/config/shared_locus_rule_v1.json',
        'placo_five_track_lock': ROOT / 'brain6/config/placo_family_v3/family_lock.json',
        'bounded_regions': BOUNDED / 'region_evidence_20.tsv',
        'bounded_replication': BOUNDED / 'replication_25.tsv',
        'bounded_regulatory': BOUNDED / 'regulatory_gene_25.tsv',
        'bounded_tissue': BOUNDED / 'tissue_context_25.tsv',
        'bounded_pathways': BOUNDED / 'pathway_20.tsv',
        'bounded_literature': BOUNDED / 'prior_literature_records.tsv',
        'recent_literature_context': ROOT / 'brain6/results/novelty/literature_update_20260926.tsv',
        'gwas_source_manifest': ROOT / 'brain6/manifests/gwas_external_availability.tsv',
        'signed_directions': ROOT / 'brain6/results/brain6_candidate_signed_direction_v1/candidate_directions_25.tsv',
        'signed_direction_receipt': ROOT / 'brain6/results/brain6_candidate_signed_direction_v1/receipt.json',
        'canonical_lava': CANON,
        'locked_global_map': ROOT / 'brain6/results/global/brain6_72_locked.tsv',
        'family_projection': PROJECTION,
        'insomnia_v4_audit': ROOT / 'brain6/results/confirmatory_source_rescue_20260927/other_traits/insomnia_native_n_pilot_v4_audit/summary.json',
        'recovery_map_receipt': ROOT / 'brain6/results/confirmatory_source_rescue_20260927/phase2_recovery_map_v2/receipt.json',
        'alternative_candidates': ALT / 'candidate_crosswalk.tsv',
        'alternative_regions': ALT / 'region_crosswalk.tsv',
        'alternative_blocks': ALT / 'genome_wide_blocks.tsv',
        'alternative_provenance': ALT / 'run_provenance.json',
        'alternative_source_eligibility': ALT / 'source_eligibility_lock_v1.json',
        'fine_mapping': FM / 'fine_mapping_50.tsv',
        'fine_mapping_credible_sets': FM / 'fine_mapping_credible_sets.tsv',
        'fine_mapping_integrity': FM / 'integrity_provenance.json',
        'trait_coloc': FM / 'trait_coloc_25.tsv',
        'trait_coloc_prior_sensitivity': FM / 'trait_coloc_prior_sensitivity.tsv',
        'functional_candidates': FUNCTION / 'candidate_functional_status_25.tsv',
        'functional_regions': FUNCTION / 'geographic_region_status_20.tsv',
        'functional_coloc': FUNCTION / 'coloc_results.tsv',
        'functional_gate': FUNCTION / 'coloc_gate.tsv',
        'functional_provenance': FUNCTION / 'provenance.json',
        'functional_qtl_input_receipts': FUNCTION / 'coloc_input_receipts.json',
        'enrichment_matching': ENRICH / 'matching_status_20.tsv',
        'enrichment_tissues': ENRICH / 'tissue_bulk_gtex_v8_enrichment_18.tsv',
        'enrichment_pathways': ENRICH / 'reactome_positional_enrichment_18.tsv',
        'enrichment_provenance': ENRICH / 'run_summary.json',
        'enrichment_file_provenance': ENRICH / 'provenance.json',
        'enrichment_frozen_binding': ENRICH / 'frozen_binding.json',
    }
    for name, path in source_paths.items():
        if not path.is_file():
            raise FileNotFoundError(f'Required completed source {name}: {path}')
    for name, expected in (('canonical_lava', FROZEN_CANON_SHA256),
                           ('locked_global_map', FROZEN_GLOBAL_SHA256),
                           ('protected_candidates', FROZEN_CANDIDATE_SHA256),
                           ('canonical_lava_config', FROZEN_CANON_CONFIG_SHA256),
                           ('shared_locus_rule', FROZEN_SHARED_RULE_SHA256),
                           ('placo_five_track_lock', FROZEN_PLACO_LOCK_SHA256)):
        if sha(source_paths[name]) != expected:
            raise ValueError(f'Frozen anchor hash changed: {name}')
    base = unique(table(source_paths['bounded_candidates']), 'candidate_locus_id', 25)
    protected = unique(table(source_paths['protected_candidates']), 'locus_id', 25)
    geography = unique(table(source_paths['bounded_regions']), 'geographic_region_grch37', 20)
    alt = unique(table(source_paths['alternative_candidates']), 'candidate_locus_id', 25)
    alt_region = unique(table(source_paths['alternative_regions']), 'region_group', 20)
    alt_blocks = table(source_paths['alternative_blocks'])
    if len(alt_blocks) != 8465 or Counter(r['pair_id'] for r in alt_blocks) != {pair: 1693 for pair in
        ('insomnia__adhd', 'insomnia__mdd', 'longsleep__scz', 'longsleep__bipolar', 'longsleep__parkinson')}:
        raise ValueError('Secondary five-pair by 1,693-block family incomplete')
    eligibility_lock = json.loads(source_paths['alternative_source_eligibility'].read_text())
    eligibility = eligibility_lock['pair_status']
    if eligibility_lock['family_denominator_pairs'] != 5 or eligibility_lock['blocks_per_pair'] != 1693 or abs(
        eligibility_lock['alpha_block'] - 0.05 / 8465) > 1e-18:
        raise ValueError('Frozen secondary multiple-testing family changed')
    for block in alt_blocks:
        if block['pair_id'].startswith('longsleep__') and (eligibility[block['pair_id']] != 'METHOD_INAPPLICABLE'
               or block['analysis_status'] != 'METHOD_INAPPLICABLE'):
            raise ValueError('Long-sleep secondary method inapplicability changed')
    fm = table(source_paths['fine_mapping'])
    if len(fm) != 50 or len({(r['candidate_locus_id'], r['trait_id']) for r in fm}) != 50:
        raise ValueError('Expected 50 unique candidate-trait fine-map rows')
    fm_by = {(r['candidate_locus_id'], r['trait_id']): r for r in fm}
    coloc = unique(table(source_paths['trait_coloc']), 'candidate_locus_id', 25)
    prior = table(source_paths['trait_coloc_prior_sensitivity'])
    low_prior = unique([r for r in prior if r['p12'] == '1e-06'], 'candidate_locus_id', 6)
    functional = unique(table(source_paths['functional_candidates']), 'candidate_locus_id', 25)
    signed = unique(table(source_paths['signed_directions']), 'candidate_locus_id', 25)
    functional_region = unique(table(source_paths['functional_regions']), 'region_grch37', 20)
    replication = unique(table(source_paths['bounded_replication']), 'candidate_locus_id', 25)
    regulatory = unique(table(source_paths['bounded_regulatory']), 'candidate_locus_id', 25)
    tissue = unique(table(source_paths['bounded_tissue']), 'candidate_locus_id', 25)
    pathways = unique(table(source_paths['bounded_pathways']), 'geographic_region_grch37', 20)
    qtl = table(source_paths['functional_coloc'])
    qtl_by = defaultdict(list)
    for result in qtl:
        if result['candidate_locus_id'] not in base or result['qtl_context'].lower().split('_')[-1] not in {'eqtl', 'sqtl'}:
            raise ValueError('QTL result has unknown candidate or context')
        qtl_by[(result['candidate_locus_id'], result['qtl_context'].lower().split('_')[-1])].append(result)
    matching = unique(table(source_paths['enrichment_matching']), 'region_group', 20)
    enrichment_tissues = table(source_paths['enrichment_tissues'])
    enrichment_pathways = table(source_paths['enrichment_pathways'])
    if len(enrichment_tissues) != 54 or len(enrichment_pathways) != 1680:
        raise ValueError('18-region sensitivity family size changed')
    if Counter(r['v3_subset'] for r in matching.values()) != {'INCLUDED_18': 18, 'EXCLUDED_MATCHING_GATE': 2}:
        raise ValueError('18-region subset gate changed')
    lava = table(CANON)
    global_map = table(source_paths['locked_global_map'])
    if len(global_map) != 72 or sum(r['significance_under_original_396_family'] == 'True' for r in global_map) != 35:
        raise ValueError('Locked 12-by-6 global map changed')
    if len(lava) != 17465 or Counter(r['status'] for r in lava) != {'TESTED': 13745, 'NOT_RUN': 3720}:
        raise ValueError('Frozen canonical LAVA changed')
    projection = json.loads(PROJECTION.read_text())
    if projection['validated_projected_not_run'] != 3720 or projection['frozen_maximum_not_run'] != 873:
        raise ValueError('Frozen LAVA family gate changed')
    if set(base) != set(protected) or set(base) != set(alt) or set(base) != set(coloc) or set(base) != set(functional) or set(base) != set(signed):
        raise ValueError('One of the 25-candidate analyses changed membership')
    if set(geography) != set(functional_region) or set(geography) != set(alt_region):
        raise ValueError('One of the 20-region analyses changed geography')
    if set(geography) != set(matching):
        raise ValueError('Enrichment geography changed')
    region_members = [key for row in geography.values() for key in row['candidate_locus_ids'].split(';')]
    if Counter(region_members) != Counter({key: 1 for key in base}):
        raise ValueError('25 candidates are not partitioned once across 20 regions')
    for key, row in base.items():
        pair, region, leads = row['pair_id'], row['geographic_region_grch37'], row['lead_variants']
        if protected[key]['pair_id'] != pair or protected[key]['lead_variants'] != leads:
            raise ValueError(f'Protected identity mismatch: {key}')
        if pair.startswith('longsleep__') and alt[key]['status'] != 'METHOD_INAPPLICABLE':
            raise ValueError(f'Long-sleep candidate unexpectedly estimated: {key}')
        if alt[key]['source_eligibility'] != eligibility[pair]:
            raise ValueError(f'Secondary source-eligibility mismatch: {key}')
        if alt[key]['status'] == 'SUPPORTED_FWER' and (eligibility[pair] != 'ELIGIBLE_SECONDARY'
                or numeric(alt[key]['best_p']) is None
                or numeric(alt[key]['best_p']) > eligibility_lock['alpha_block']):
            raise ValueError(f'Unsupported secondary FWER promotion: {key}')
        for other, region_key in ((functional[key], 'region_grch37'), (signed[key], None),
                                  (alt[key], None), (coloc[key], None)):
            if other['pair_id'] != pair or (region_key and other[region_key] != region):
                raise ValueError(f'Analysis identity mismatch: {key}')
        if signed[key]['lead_variants'] != leads or functional[key]['lead_variants'] != leads:
            raise ValueError(f'Lead identity mismatch: {key}')
        if key not in geography[region]['candidate_locus_ids'].split(';'):
            raise ValueError(f'Geographic identity mismatch: {key}')
        if replication[key]['independent_two_trait_locus_replication'] != functional[key]['independent_two_trait_locus_replication']:
            raise ValueError(f'Replication source conflict: {key}')
    if any({(r['candidate_locus_id'], trait) for trait in r['pair_id'].split('__')} - set(fm_by) for r in base.values()):
        raise ValueError('Fine-map rows missing a candidate component trait')
    # The primary family failed and cannot be replaced by a secondary method.
    final_lava = [dict(r, family_decision='FAILED_QC_NOT_PROMOTED', rescue_status='NOT_ADMITTED') for r in lava]
    candidates = []
    for key, row in base.items():
        ar, cr, fr = alt[key], coloc[key], functional[key]
        left, right = row['pair_id'].split('__')
        left_fm, right_fm = fm_by[(key, left)], fm_by[(key, right)]
        eqtl_tests = qtl_by[(key, 'eqtl')]
        sqtl_tests = qtl_by[(key, 'sqtl')]
        eqtl_status = 'EXPLORATORY_COMPONENT_ABF_ONLY' if eqtl_tests else 'NOT_ESTIMATED_NO_ELIGIBLE_EQTL_CONTEXT'
        sqtl_status = 'EXPLORATORY_COMPONENT_ABF_ONLY' if sqtl_tests else 'NOT_ESTIMATED_NO_ELIGIBLE_SQTL_CONTEXT'
        low_pp = low_prior.get(key, {}).get('PP.H4', 'NA')
        prior_robust = ('ROBUST_TO_LOW_P12' if cr.get('descriptive_support') == 'ABF_MODEL_SUPPORT'
                        and numeric(low_pp) is not None and numeric(low_pp) >= 0.8
                        else 'PRIOR_SENSITIVE' if cr.get('descriptive_support') == 'ABF_MODEL_SUPPORT'
                        else 'NO_DEFAULT_ABF_MODEL_SUPPORT')
        if ar['status'] == 'SUPPORTED_FWER':
            evidence_class = 'ALTERNATIVE_LOCAL_VALIDATION_SUPPORTED'
        elif fr['independent_two_trait_locus_replication'] in {'INDEPENDENT_TWO_TRAIT_LOCUS_REPLICATION_SUPPORTED'}:
            evidence_class = 'PLACO_PLUS_REPLICATION_SUPPORTED'
        elif fr['gwas_qtl_coloc_status'] == 'SUPPORTED_SHARED_VARIANT_QTL':
            evidence_class = 'EXPLORATORY_MULTIOMIC_SUPPORTED'
        else:
            evidence_class = 'UNSUPPORTED'
        candidates.append({
            'candidate_locus_id': key, 'geographic_region_grch37': row['geographic_region_grch37'],
            'pair_id': row['pair_id'], 'lead_variants': row['lead_variants'],
            'placo_lead_p': row['placo_lead_p'],
            'n_candidate_variants': protected[key]['n_candidate_variants'],
            'ld_reference': protected[key]['reference'],
            'ld_r2_threshold': protected[key]['r2_threshold'],
            'clump_window_kb': protected[key]['window_kb'],
            'placo_lead_signed_relations': signed[key]['direction_relations'],
            'placo_lead_signed_z_products': signed[key]['signed_z_products'],
            'canonical_lava_family': 'FAILED_QC_NOT_PROMOTED',
            'canonical_trait1_status': row['canonical_trait1_status'],
            'canonical_trait2_status': row['canonical_trait2_status'],
            'rescued_lava': 'NOT_ADMITTED',
            'alternative_local_status': ar['status'], 'alternative_local_min_p': ar['best_p'],
            'alternative_local_direction': ar['covariance_sign'],
            'left_fine_mapping_status': left_fm['status'], 'right_fine_mapping_status': right_fm['status'],
            'trait_coloc_status': cr['status'], 'trait_coloc_pp_h4': cr.get('PP.H4', 'NA'),
            'trait_coloc_descriptive_support': cr.get('descriptive_support', 'NA'),
            'trait_coloc_low_prior_pp_h4': low_pp,
            'trait_coloc_prior_sensitivity': prior_robust,
            'eqtl_coloc_status': eqtl_status,
            'eqtl_n_component_tests': len(eqtl_tests),
            'max_component_qtl_pp_h4': fr['max_exploratory_pp_h4'],
            'sqtl_coloc_status': sqtl_status,
            'sqtl_n_component_tests': len(sqtl_tests),
            'pair_level_gwas_qtl_coloc_status': 'NOT_ESTIMATED',
            'published_eqtl_cs_memberships': fr['published_eqtl_cs_memberships'],
            'published_sqtl_cs_memberships': fr['published_sqtl_cs_memberships'],
            'exact_lead_regulatory_overlaps': fr['exact_lead_regulatory_overlaps'],
            'independent_two_trait_locus_replication': fr['independent_two_trait_locus_replication'],
            'global_pair_replication_class': fr['global_pair_replication_class'],
            'tissue_cell_type_enrichment': fr['tissue_cell_type_enrichment'],
            'pathway_enrichment': fr['pathway_enrichment'],
            'subset_positional_enrichment_status': matching[row['geographic_region_grch37']]['v3_subset'],
            'nearest_positional_gene_symbols': regulatory[key]['nearest_positional_gene_symbols'],
            'prior_trait_relevant_catalog_rows': regulatory[key]['prior_trait_relevant_catalog_rows'],
            'final_evidence_class': evidence_class,
            'class_interpretation': 'HIGHEST_ELIGIBLE_PAIR_SPECIFIC_TIER;UNSUPPORTED_MEANS_NO_A_TO_D_TIER_MET',
            'limitation_flags': 'PRIMARY_LAVA_FAILED_QC;NO_FORMAL_PROMOTION;NO_CAUSAL_VARIANT_OR_GENE_CLAIM',
        })
    out_by = {r['candidate_locus_id']: r for r in candidates}
    regions = []
    precedence = EVIDENCE_CLASSES
    for region, row in geography.items():
        ids = row['candidate_locus_ids'].split(';')
        if len(ids) != int(row['pair_specific_candidates']) or any(key not in out_by for key in ids):
            raise ValueError(f'Geographic region member mismatch: {region}')
        members = [out_by[key] for key in ids]
        klass = next(k for k in precedence if any(m['final_evidence_class'] == k for m in members))
        pp_trait = [numeric(m['trait_coloc_pp_h4']) for m in members]
        pp_qtl = [numeric(m['max_component_qtl_pp_h4']) for m in members]
        regions.append({'geographic_region_grch37': region, 'pair_specific_candidates': len(ids),
                        'pairs': row['pairs'], 'candidate_locus_ids': row['candidate_locus_ids'],
                        'lead_variants': row['lead_variants'], 'best_placo_p': row['best_placo_p'],
                        'candidate_variant_rows_pair_counted': row['candidate_variant_rows'],
                        'placo_lead_signed_relations': ';'.join(m['placo_lead_signed_relations'] for m in members),
                        'canonical_lava_family': 'FAILED_QC_NOT_PROMOTED',
                        'rescued_lava': 'NOT_ADMITTED',
                        'alternative_local_statuses': ';'.join(sorted({m['alternative_local_status'] for m in members})),
                        'fine_mapping_statuses': ';'.join(sorted({m['left_fine_mapping_status'] for m in members} | {m['right_fine_mapping_status'] for m in members})),
                        'trait_coloc_statuses': ';'.join(sorted({m['trait_coloc_status'] for m in members})),
                        'n_abf_model_supported_trait_coloc_candidates': sum(m['trait_coloc_descriptive_support'] == 'ABF_MODEL_SUPPORT' for m in members),
                        'n_low_prior_robust_abf_candidates': sum(m['trait_coloc_prior_sensitivity'] == 'ROBUST_TO_LOW_P12' for m in members),
                        'max_exploratory_trait_pp_h4': max((x for x in pp_trait if x is not None), default='NA'),
                        'max_exploratory_component_qtl_pp_h4': max((x for x in pp_qtl if x is not None), default='NA'),
                        'eqtl_coloc_statuses': ';'.join(sorted({m['eqtl_coloc_status'] for m in members})),
                        'sqtl_coloc_statuses': ';'.join(sorted({m['sqtl_coloc_status'] for m in members})),
                        'pair_level_gwas_qtl_coloc_status': 'NOT_ESTIMATED',
                        'independent_two_trait_locus_replication': functional_region[region]['independent_two_trait_locus_replication'],
                        'published_eqtl_cs_membership_rows_pair_counted': functional_region[region]['published_eqtl_cs_membership_rows_pair_counted'],
                        'published_sqtl_cs_membership_rows_pair_counted': functional_region[region]['published_sqtl_cs_membership_rows_pair_counted'],
                        'exact_lead_regulatory_feature_rows_pair_counted': functional_region[region]['exact_lead_regulatory_feature_rows_pair_counted'],
                        'nearest_positional_gene_symbols': ';'.join(sorted({regulatory[key]['nearest_positional_gene_symbols'] for key in ids if regulatory[key]['nearest_positional_gene_symbols'] not in ('', 'NONE', 'NA')})) or 'NA',
                        'prior_trait_relevant_catalog_rows_pair_counted': row['prior_trait_relevant_catalog_rows'],
                        'tissue_cell_type_pathway_inference': functional_region[region]['tissue_cell_type_pathway_inference'],
                        'subset_positional_enrichment_status': matching[region]['v3_subset'],
                        'final_evidence_class': klass,
                        'highest_pair_specific_evidence_class': klass,
                        'class_unit': 'HIGHEST_MEMBER_PAIR_SPECIFIC_CLASS_NOT_ALL_MEMBERS',
                        'limitation_flags': 'PRIMARY_LAVA_FAILED_QC;GEOGRAPHIC_MERGE_DESCRIPTIVE;NO_CAUSAL_CLAIM'})
    if len(regions) != 20 or len(candidates) != 25:
        raise ValueError('Final geographic/candidate count mismatch')
    if any(numeric(r['q_bh_54']) is None for r in enrichment_tissues) or any(numeric(r['q_bh_1680']) is None for r in enrichment_pathways):
        raise ValueError('Enrichment BH q values missing')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.final_package_v1_staging_', dir=OUT.parent))
    write(stage / 'BRAIN6_FINAL_REGIONS.tsv', regions)
    write(stage / 'BRAIN6_FINAL_CANDIDATES.tsv', candidates)
    write(stage / 'BRAIN6_FINAL_LAVA.tsv', final_lava)
    write(stage / 'BRAIN6_FINAL_GLOBAL.tsv', global_map)
    write(stage / 'BRAIN6_FINAL_ALT_LOCAL_VALIDATION.tsv', list(alt.values()))
    write(stage / 'BRAIN6_FINAL_ALT_LOCAL_BLOCKS.tsv', alt_blocks)
    write(stage / 'BRAIN6_FINAL_ALT_LOCAL_REGIONS.tsv', list(alt_region.values()))
    write(stage / 'BRAIN6_FINAL_FINE_MAPPING.tsv', fm)
    write(stage / 'BRAIN6_FINAL_TRAIT_COLOC.tsv', list(coloc.values()))
    write(stage / 'BRAIN6_FINAL_TRAIT_COLOC_PRIOR_SENSITIVITY.tsv', prior)
    # A status row for every candidate is retained alongside any numerical
    # component-GWAS molecular test. They are different record types.
    for qtype, label, path in (('eqtl', 'EQTL', 'BRAIN6_FINAL_EQTL_COLOC.tsv'),
                               ('sqtl', 'SQTL', 'BRAIN6_FINAL_SQTL_COLOC.tsv')):
        rows = []
        for key, fr in functional.items():
            rows.append({'record_type': 'CANDIDATE_STATUS', 'candidate_locus_id': key,
                         'pair_id': fr['pair_id'], 'gwas_trait': 'NA', 'qtl_context': qtype,
                         'molecular_trait_id': 'NA', 'gene_id': 'NA', 'n_shared': 'NA',
                         'pp_h0': 'NA', 'pp_h1': 'NA', 'pp_h2': 'NA', 'pp_h3': 'NA', 'pp_h4': 'NA',
                         'status': out_by[key][f'{qtype}_coloc_status'],
                         'limitation': 'COMPONENT_GWAS_ONLY_WHERE_ESTIMATED;NO_PAIR_LEVEL_QTL_COLOC'})
        for r in qtl:
            if qtype in r['qtl_context'].lower():
                rows.append({'record_type': 'NUMERICAL_COMPONENT_TEST',
                             'candidate_locus_id': r['candidate_locus_id'],
                             'pair_id': base[r['candidate_locus_id']]['pair_id'],
                             'gwas_trait': r['gwas_trait'], 'qtl_context': r['qtl_context'],
                             'molecular_trait_id': r['molecular_trait_id'], 'gene_id': r['gene_id'],
                             'n_shared': r['n_shared'], 'pp_h0': r['pp_h0'], 'pp_h1': r['pp_h1'],
                             'pp_h2': r['pp_h2'], 'pp_h3': r['pp_h3'], 'pp_h4': r['pp_h4'],
                             'status': r['status'], 'limitation': 'SINGLE_SIGNAL_COMPONENT_GWAS_EXPLORATORY'})
        write(stage / path, rows)
    final_replication = [dict(replication[k], independent_two_trait_locus_replication=functional[k]['independent_two_trait_locus_replication']) for k in base]
    final_regulatory = [dict(regulatory[k], published_eqtl_cs_memberships=functional[k]['published_eqtl_cs_memberships'], published_sqtl_cs_memberships=functional[k]['published_sqtl_cs_memberships'], exact_lead_regulatory_overlaps=functional[k]['exact_lead_regulatory_overlaps']) for k in base]
    final_tissue = ([dict(record_type='CANDIDATE_CONTEXT_STATUS', **tissue[k], current_enrichment_status=functional[k]['tissue_cell_type_enrichment'],
                          subset_positional_enrichment_status=matching[base[k]['geographic_region_grch37']]['v3_subset']) for k in base]
                    + [dict(record_type='GTEX_BULK_TISSUE_18_REGION_TEST', **r) for r in enrichment_tissues])
    final_pathways = ([dict(record_type='REGION_CONTEXT_STATUS', **pathways[k], current_enrichment_status=functional_region[k]['tissue_cell_type_pathway_inference'],
                            subset_positional_enrichment_status=matching[k]['v3_subset']) for k in geography]
                      + [dict(record_type='REACTOME_POSITIONAL_18_REGION_TEST', **r) for r in enrichment_pathways])
    write(stage / 'BRAIN6_FINAL_REPLICATION.tsv', final_replication)
    write(stage / 'BRAIN6_FINAL_REGULATORY.tsv', final_regulatory)
    write(stage / 'BRAIN6_FINAL_TISSUE_CELLTYPE.tsv', final_tissue,
          ['record_type'] + sorted({f for r in final_tissue for f in r if f != 'record_type'}))
    write(stage / 'BRAIN6_FINAL_PATHWAYS.tsv', final_pathways,
          ['record_type'] + sorted({f for r in final_pathways for f in r if f != 'record_type'}))
    literature = ([dict(record_type='REGION_CATALOG_CONTEXT', **r) for r in table(source_paths['bounded_literature'])]
                  + [dict(record_type='RECENT_SOURCE_CONTEXT', **r) for r in table(source_paths['recent_literature_context'])])
    literature_fields = ['record_type'] + sorted({field for r in literature for field in r if field != 'record_type'})
    write(stage / 'BRAIN6_FINAL_LITERATURE.tsv', literature, literature_fields)
    outputs = sorted(stage.glob('*.tsv'))
    receipt = {'schema_version': 1, 'package_status': 'BOUNDED_SECONDARY_VALIDATION_PAPER',
               'canonical_family_decision': 'FAILED_QC_NOT_PROMOTED',
               'candidate_count': 25, 'geographic_region_count': 20,
               'candidate_class_counts': {klass: sum(r['final_evidence_class'] == klass for r in candidates)
                                          for klass in EVIDENCE_CLASSES},
               'source_sha256': {source_label(p): sha(p) for p in source_paths.values()},
               'output_sha256': {p.name: sha(p) for p in outputs},
               'builder_sha256': sha(Path(__file__))}
    with (stage / 'BRAIN6_FINAL_PROVENANCE.json').open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    os.replace(stage, OUT)
    print(json.dumps({'candidate_count': 25, 'region_count': 20,
                      'class_counts': receipt['candidate_class_counts'], 'tables': len(outputs)}, sort_keys=True))


if __name__ == '__main__':
    main()
