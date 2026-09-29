#!/usr/bin/env python3
"""Independently check the Brain6 paper package against its frozen sources."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'brain6/paper/final_package_v1'
DOCUMENTS = ('BRAIN6_FINAL_RESULTS.md', 'BRAIN6_FINAL_METHODS.md',
             'BRAIN6_FINAL_DISCUSSION.md', 'BRAIN6_FINAL_LIMITATIONS.md',
             'BRAIN6_FINAL_CLAIMS.md', 'BRAIN6_FINAL_READINESS.md')


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def audit() -> dict[str, object]:
    provenance = json.loads((PACKAGE / 'BRAIN6_FINAL_PROVENANCE.json').read_text())
    assert provenance['candidate_count'] == 25
    assert provenance['geographic_region_count'] == 20
    assert provenance['canonical_family_decision'] == 'FAILED_QC_NOT_PROMOTED'
    for name, expected in provenance['source_sha256'].items():
        path = Path(name) if name.startswith('/') else ROOT / name
        assert path.is_file() and sha(path) == expected, f'Source hash mismatch: {name}'
    for name, expected in provenance['output_sha256'].items():
        path = PACKAGE / name
        assert path.is_file() and sha(path) == expected, f'Output hash mismatch: {name}'
    assert sha(ROOT / 'brain6/scripts/build_brain6_final_package_v1.py') == provenance['builder_sha256']
    upstream = (
        (ROOT / 'brain6/results/brain6_alternative_local_validation_v1',
         'run_provenance.json', 'output_sha256'),
        (ROOT / 'brain6/results/brain6_exploratory_finemap_coloc_v2',
         'integrity_provenance.json', 'result_sha256'),
        (ROOT / 'brain6/results/brain6_exploratory_functional_v1',
         'provenance.json', 'output_sha256'),
        (ROOT / 'brain6/results/brain6_exploratory_enrichment_v3',
         'provenance.json', 'outputs_sha256'),
    )
    for directory, receipt_name, key in upstream:
        receipt = json.loads((directory / receipt_name).read_text())
        for filename, expected in receipt[key].items():
            assert sha(directory / filename) == expected, f'Upstream receipt mismatch: {directory.name}/{filename}'
    gwas_sources = rows(ROOT / 'brain6/manifests/gwas_external_availability.tsv')
    primary_traits = {'insomnia', 'longsleep', 'adhd', 'mdd', 'scz', 'bipolar', 'parkinson'}
    assert primary_traits <= {r['trait'] for r in gwas_sources}
    assert all(r['source_URL'].startswith('https://') and len(r['raw_file_SHA256']) == 64
               for r in gwas_sources if r['trait'] in primary_traits)
    enrichment_sources = rows(ROOT / 'brain6/results/brain6_exploratory_enrichment_v2/source_receipts.tsv')
    assert len(enrichment_sources) >= 3
    assert sum(r['source_url'].startswith('https://') for r in enrichment_sources) >= 3
    assert all((r['source_url'].startswith('https://') or r['source_url'].startswith('repository:'))
               and len(r['sha256']) == 64 for r in enrichment_sources)

    candidates = rows(PACKAGE / 'BRAIN6_FINAL_CANDIDATES.tsv')
    regions = rows(PACKAGE / 'BRAIN6_FINAL_REGIONS.tsv')
    lava = rows(PACKAGE / 'BRAIN6_FINAL_LAVA.tsv')
    global_map = rows(PACKAGE / 'BRAIN6_FINAL_GLOBAL.tsv')
    alternative = rows(PACKAGE / 'BRAIN6_FINAL_ALT_LOCAL_VALIDATION.tsv')
    alternative_blocks = rows(PACKAGE / 'BRAIN6_FINAL_ALT_LOCAL_BLOCKS.tsv')
    alternative_regions = rows(PACKAGE / 'BRAIN6_FINAL_ALT_LOCAL_REGIONS.tsv')
    fine = rows(PACKAGE / 'BRAIN6_FINAL_FINE_MAPPING.tsv')
    trait = rows(PACKAGE / 'BRAIN6_FINAL_TRAIT_COLOC.tsv')
    prior = rows(PACKAGE / 'BRAIN6_FINAL_TRAIT_COLOC_PRIOR_SENSITIVITY.tsv')
    eqtl = rows(PACKAGE / 'BRAIN6_FINAL_EQTL_COLOC.tsv')
    sqtl = rows(PACKAGE / 'BRAIN6_FINAL_SQTL_COLOC.tsv')
    tissue = rows(PACKAGE / 'BRAIN6_FINAL_TISSUE_CELLTYPE.tsv')
    pathways = rows(PACKAGE / 'BRAIN6_FINAL_PATHWAYS.tsv')
    assert (len(candidates), len(regions), len(lava), len(global_map), len(fine), len(trait)) == (25, 20, 17465, 72, 50, 25)
    assert len({r['candidate_locus_id'] for r in candidates}) == 25
    assert len({r['geographic_region_grch37'] for r in regions}) == 20
    assert len(alternative) == 25
    assert len(alternative_blocks) == 8465 and len(alternative_regions) == 20
    assert Counter(r['pair_id'] for r in alternative_blocks) == {
        pair: 1693 for pair in ('insomnia__adhd', 'insomnia__mdd', 'longsleep__scz',
                               'longsleep__bipolar', 'longsleep__parkinson')}
    assert sum(r['analysis_status'] == 'METHOD_INAPPLICABLE' for r in alternative_blocks) == 5079
    assert Counter(r['status'] for r in lava) == {'TESTED': 13745, 'NOT_RUN': 3720}
    assert all(r['family_decision'] == 'FAILED_QC_NOT_PROMOTED' and r['rescue_status'] == 'NOT_ADMITTED' for r in lava)
    assert sum(r['significance_under_original_396_family'] == 'True' for r in global_map) == 35
    assert len(prior) == 12
    assert Counter(r['record_type'] for r in eqtl) == {'CANDIDATE_STATUS': 25, 'NUMERICAL_COMPONENT_TEST': 3}
    assert Counter(r['record_type'] for r in sqtl) == {'CANDIDATE_STATUS': 25, 'NUMERICAL_COMPONENT_TEST': 1}
    assert Counter(r['record_type'] for r in tissue) == {'CANDIDATE_CONTEXT_STATUS': 25, 'GTEX_BULK_TISSUE_18_REGION_TEST': 54}
    assert Counter(r['record_type'] for r in pathways) == {'REGION_CONTEXT_STATUS': 20, 'REACTOME_POSITIONAL_18_REGION_TEST': 1680}
    assert all(float(r['q_bh_54']) >= 0.05 for r in tissue if r['record_type'] == 'GTEX_BULK_TISSUE_18_REGION_TEST')
    assert all(float(r['q_bh_1680']) >= 0.05 for r in pathways if r['record_type'] == 'REACTOME_POSITIONAL_18_REGION_TEST')

    alternative_by = {r['candidate_locus_id']: r for r in alternative}
    assert set(alternative_by) == {r['candidate_locus_id'] for r in candidates}
    allowed = {'PRIMARY_LAVA_CONFIRMED', 'ALTERNATIVE_LOCAL_VALIDATION_SUPPORTED',
               'PLACO_PLUS_REPLICATION_SUPPORTED', 'EXPLORATORY_MULTIOMIC_SUPPORTED', 'UNSUPPORTED'}
    assert all(r['final_evidence_class'] in allowed for r in candidates)
    assert set(provenance['candidate_class_counts']) == allowed
    assert sum(provenance['candidate_class_counts'].values()) == 25
    assert all(provenance['candidate_class_counts'][klass] == sum(r['final_evidence_class'] == klass for r in candidates)
               for klass in allowed)
    assert all(r['final_evidence_class'] != 'PRIMARY_LAVA_CONFIRMED' for r in candidates)
    assert all(r['canonical_lava_family'] == 'FAILED_QC_NOT_PROMOTED' and r['rescued_lava'] == 'NOT_ADMITTED' for r in candidates)
    assert all(r['pair_level_gwas_qtl_coloc_status'] == 'NOT_ESTIMATED' for r in candidates)
    assert all(r['final_evidence_class'] != 'ALTERNATIVE_LOCAL_VALIDATION_SUPPORTED'
               or alternative_by[r['candidate_locus_id']]['status'] == 'SUPPORTED_FWER' for r in candidates)
    assert all(r['final_evidence_class'] != 'PLACO_PLUS_REPLICATION_SUPPORTED'
               or r['independent_two_trait_locus_replication'] == 'INDEPENDENT_TWO_TRAIT_LOCUS_REPLICATION_SUPPORTED'
               for r in candidates)
    assert all(r['final_evidence_class'] != 'EXPLORATORY_MULTIOMIC_SUPPORTED'
               for r in candidates), 'Current source has no supported pair-level molecular shared-variant result'
    assert Counter(r['trait_coloc_descriptive_support'] for r in candidates)['ABF_MODEL_SUPPORT'] == 3
    assert Counter(r['trait_coloc_prior_sensitivity'] for r in candidates)['ROBUST_TO_LOW_P12'] == 1
    candidate_ids = Counter(key for r in regions for key in r['candidate_locus_ids'].split(';'))
    assert candidate_ids == Counter({r['candidate_locus_id']: 1 for r in candidates})
    assert Counter(r['subset_positional_enrichment_status'] for r in regions) == {'INCLUDED_18': 18, 'EXCLUDED_MATCHING_GATE': 2}

    document_hashes = {}
    for name in DOCUMENTS:
        path = PACKAGE / name
        assert path.is_file() and path.stat().st_size > 0, f'Missing document {name}'
        document_hashes[name] = sha(path)
    results_text = (PACKAGE / 'BRAIN6_FINAL_RESULTS.md').read_text()
    required_results = ('17,465', '13,745', '3,720', '873', '25 pair-specific', '20 overlapping',
                        '47 `TESTED`', '41 `NOT_RUN`', '3,304', '5,079', '24 model-conditional',
                        'maximum PP.H4 was **0.0824**', '0/25', '54 GTEx', '1,680 Reactome',
                        'A `PRIMARY_LAVA_CONFIRMED`: 0', 'E `UNSUPPORTED`: 25')
    assert all(value in results_text for value in required_results), 'Results text missing a required audited count'
    claims_text = (PACKAGE / 'BRAIN6_FINAL_CLAIMS.md').read_text()
    assert all(status in claims_text for status in ('PRIMARY_CONFIRMATORY_SUPPORTED',
        'SECONDARY_VALIDATION_SUPPORTED', 'EXPLORATORY_SUPPORTED', 'DESCRIPTIVE_ONLY',
        'UNSUPPORTED_DO_NOT_CLAIM'))
    return {'status': 'PASS', 'canonical_family_decision': 'FAILED_QC_NOT_PROMOTED',
            'candidate_count': 25, 'region_count': 20,
            'candidate_class_counts': provenance['candidate_class_counts'],
            'alternative_status_counts': dict(Counter(r['status'] for r in alternative)),
            'output_table_count': len(provenance['output_sha256']),
            'document_sha256': document_hashes,
            'audit_script_sha256': sha(Path(__file__)),
            'provenance_sha256': sha(PACKAGE / 'BRAIN6_FINAL_PROVENANCE.json')}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true', help='write immutable audit receipt after passing checks')
    args = parser.parse_args()
    result = audit()
    if args.write:
        path = PACKAGE / 'BRAIN6_FINAL_AUDIT.json'
        with path.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
