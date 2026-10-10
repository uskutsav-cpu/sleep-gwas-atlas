#!/usr/bin/env python3
"""Evidence-only validation. No genetics, source-body downloads, or old-file writes.

Run from any directory with Python 3. Outputs stay in this review's v4 directories.
New public-paper receipt hashes certify the read bodies, not absent cached bodies.
"""
from pathlib import Path
import csv
import hashlib
import json
import urllib.parse

BASE = Path(__file__).resolve().parents[2]
V4 = BASE / 'sleep_unified_research_v4'
REVIEW = V4 / 'reviews'
PROVENANCE = V4 / 'source_provenance'
OLD_CACHE = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/goa/work/sleep-gwas-atlas/discovery_extension/sleep_submission_evidence_v1/sources')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

expected = {
    BASE / 'discovery_extension/provenance/prior_screens/morrison_2024_supplementary_tables_1_15.xlsx': '649860bba073fe5c449fa670976759162dce13493f1a26cd915dc546cc6fb6ec',
    BASE / 'discovery_extension/provenance/prior_screens/goodman_2025_supplementary_data_1_27.xlsx': 'e204c09fc6b5448df1795fa3bd62fd3b6b786c1539f7345f4fc66c8901cab8e0',
    OLD_CACHE / 'literature_sleepchart_supplement.xlsx': '3488344573bf588d0ad748c83b5a86f5fbd487d711390b6a601661de56587dd7',
    OLD_CACHE / 'literature_multiorgan_MOESM12.xlsx': '84b832a3fba7235a712d528aa87ddd089fb68831fcf709245c0c405c2e654711',
    OLD_CACHE / 'literature_sleepchart_page.html': 'dbc134c3e94e67de95a030c1b59476e7d2c48fdc52cd3a2fffd674ce4114b97c',
    OLD_CACHE / 'literature_multiorgan_page.html': 'dd7ae2f1e7a0b648b4085b4940c829965e0f57928d997f7ff002f053ee3d99d0',
    OLD_CACHE / 'literature_multiorgan_MOESM3.xlsx': '134bb87af34768b3245055ab4e31d514ecea200f24216192f2b9ed5648f846e0',
}
inputs = []
for path, wanted in expected.items():
    actual = sha(path)
    assert actual == wanted, (str(path), actual, wanted)
    inputs.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': actual, 'expected_match': True})

audit_path = BASE / 'sleep_unified_research_v1/reviews/novelty_v1_supplement_inspection.tsv'
audit = read_tsv(audit_path)
total = sum(int(row['rows_inspected']) for row in audit)
assert total == 14568
contract_path = BASE / 'sleep_unified_research_v2/source_definition_review/provenance_v2_definition_contract.json'
contract = json.loads(contract_path.read_text())
assert contract['gwas_accession'] == 'GCST90475826'
assert contract['analysis_name'] == 'Phe_327_4.EUR.GIA'
assert contract['n_cases'] == 78566 and contract['n_controls'] == 329572
assert contract['case_rule']['mapped_icd_instances_minimum'] == 2
assert contract['control_rule']['related_phecode_exclusions_applied'] is False
counts_path = BASE / 'sleep_unified_research_v3/tables/MVP_PREPROCESSING_COUNTS.tsv'
counts = read_tsv(counts_path)
assert sum(int(row['count']) for row in counts) == 19703815
assert next(int(row['count']) for row in counts if row['metric'] == 'retained_HM3') == 826027
gates_path = BASE / 'sleep_unified_research_v3/tables/MVP_DOWNSTREAM_GATES_V3.tsv'
gates = read_tsv(gates_path)
assert next(row['status'] for row in gates if row['gate'] == 'Pair admission') == 'ZERO'

sd2_path = PROVENANCE / 'sleep_geneticist_mounier_SD2_phenotype_coverage_v4.json'
sd2 = json.loads(sd2_path.read_text())
labels = sd2['sleep_related_pair_labels']
assert sd2['label_count'] == 71
assert len(labels) == 139
insomnia_pairs = [row for row in labels if 'Insomnia' in (row['condition1'], row['condition2'])]
apnea_pairs = [row for row in labels if 'Sleep Apnoea' in (row['condition1'], row['condition2'])]
assert len(insomnia_pairs) == len(apnea_pairs) == 70
assert next(row['excel_row'] for row in insomnia_pairs if row['condition1'] == 'Chronic Obstructive Pulmonary Disease') == 729
official_path = PROVENANCE / 'sleep_geneticist_mounier_official_asset_metadata_v4.json'
official = json.loads(official_path.read_text())['metadata']
fetch_path = PROVENANCE / 'sleep_geneticist_mounier_fulltext_fetch_v4.json'
fetch = json.loads(fetch_path.read_text())
url_md5 = urllib.parse.parse_qs(urllib.parse.urlparse(official['xml_url']).query)['md5'][0]
xml_record = next(row for row in fetch['records'] if row['url'].endswith('.xml'))
assert xml_record['md5'] == url_md5
for name in ['sleep_geneticist_mounier_SD1_inspection_v4.json', 'sleep_geneticist_mounier_SD2_inspection_v4.json']:
    receipt = json.loads((PROVENANCE / name).read_text())
    matching = next(item for item in official['media_urls'] if urllib.parse.urlparse(item).path == urllib.parse.urlparse(receipt['url']).path)
    assert receipt['md5'] == urllib.parse.parse_qs(urllib.parse.urlparse(matching).query)['md5'][0]
    assert receipt['official_asset_md5_match'] is True

matrix_path = REVIEW / 'sleep_geneticist_primary_admission_matrix_v4.tsv'
matrix = read_tsv(matrix_path)
assert len(matrix) == 5
assert {row['candidate_id'] for row in matrix[:3]} == {'Q1', 'Q2', 'Q3'}
assert all(row['primary_admission'].startswith('NO_GO_CURRENT') for row in matrix)
for path in [audit_path, contract_path, counts_path, gates_path, sd2_path, official_path, fetch_path,
             BASE / 'sleep_unified_research_v1/reviews/novelty_v1_report.md',
             BASE / 'sleep_unified_research_v1/reviews/novelty_v1_sources.tsv',
             BASE / 'sleep_unified_research_v1/reviews/novelty_v1_prior_numeric_examples.tsv',
             BASE / 'sleep_unified_research_v1/reviews/novelty_v1_feasibility_matrix.tsv',
             BASE / 'sleep_unified_research_v1/reviews/sleep_phenotyping_v1.md',
             BASE / 'sleep_unified_research_v1/reviews/cohort_independence_v1.md',
             BASE / 'sleep_unified_research_v2/source_definition_review/provenance_v2_source_definition_addendum.md',
             BASE / 'sleep_unified_research_v3/reviews/statistical_mvp_preprocessing_postexecution_v3.md',
             BASE / 'sleep_unified_research_v1/tables/phenotype_and_source_metadata.tsv',
             REVIEW / 'genomicsem_specialist_independent_review_v4.md']:
    inputs.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)})

outputs = sorted(list(REVIEW.glob('sleep_geneticist*')) + list(PROVENANCE.glob('sleep_geneticist*')))
outputs = [path for path in outputs if path.name not in {'sleep_geneticist_evidence_validation_v4.json', 'sleep_geneticist_v4_SHA256SUMS'}]
result = {
    'review_date': '2026-10-09', 'status': 'PASS_EVIDENCE_VALIDATION',
    'new_primary_admitted': False, 'candidate_count': 3, 'retained_alternative_count': 2,
    'completed_prior_displayed_rows_reused': total, 'prior_rows_reextracted_this_review': 0,
    'mounier_displayed_pair_rows': 2485, 'mounier_insomnia_pair_labels': len(insomnia_pairs),
    'mounier_sleep_apnoea_pair_labels': len(apnea_pairs), 'mounier_unique_sleep_related_pairs': len(labels),
    'v2_definition_objection': 'DOCUMENTED_GIA_RULE_RESOLVED_QUALIFIED',
    'v3_hm3_retained': 826027, 'v3_gates_snapshot_only': True,
    'estimators_run': False, 'new_gwas_body_downloads': 0,
    'verification_limits': ['New public article and supplements read in memory; hashes and official MD5 references recorded, original bodies not redistributed.', 'Published partialLDSC block alignment is not certified by this reviewer.', 'No exhaustive negative-priority certificate for uninspected supplements or later unindexed work.'],
    'input_hashes': inputs,
    'output_hashes': [{'path': str(path.relative_to(BASE)), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in outputs],
}
out_path = REVIEW / 'sleep_geneticist_evidence_validation_v4.json'
out_path.write_text(json.dumps(result, indent=2) + '\n')
outputs.append(out_path)
(REVIEW / 'sleep_geneticist_v4_SHA256SUMS').write_text(''.join(sha(path) + '  ' + str(path.relative_to(BASE)) + '\n' for path in sorted(outputs)))
print(json.dumps({key: result[key] for key in ['status','new_primary_admitted','completed_prior_displayed_rows_reused','mounier_insomnia_pair_labels','mounier_sleep_apnoea_pair_labels','v3_hm3_retained']}, indent=2))
