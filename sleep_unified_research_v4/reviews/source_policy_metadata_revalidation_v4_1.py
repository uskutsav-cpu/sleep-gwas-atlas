#!/usr/bin/env python3
"""Revalidate only policy metadata/cache and exact membership; no network/data reads."""
import collections
import csv
import datetime
import hashlib
import json
import re
import stat
from pathlib import Path

P = Path(__file__).resolve().parents[1]
ROOT = P.parent
S = P / 'source_provenance'
OUT = P / 'reviews/source_policy_metadata_revalidation_receipt_v4_1.json'


def regular(p):
    p = Path(p)
    assert not p.is_symlink() and all(not q.is_symlink() for q in p.parents), str(p)
    assert stat.S_ISREG(p.stat().st_mode), str(p)
    return p


def sha(p):
    h = hashlib.sha256()
    with regular(p).open('rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def load(p):
    return json.loads(regular(p).read_text())


def table(p):
    with regular(p).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


evidence_path = S / 'source_policy_review_evidence_receipt_v4_1.json'
evidence = load(evidence_path)
original = evidence['input_metadata_sha256_before']
assert original == evidence['input_metadata_sha256_after']
observed = {path: sha(path) for path in original}
assert observed == original
assert sha(S / 'source_rights_evidence_gaps_v4.tsv') == 'e218575e8ffc66e26b9a16ef9d43e7b7a830c66a01839c2491386e2a68b2b31e'
assert sha(P / 'reviews/validation13_acquisition_independent_review_seal_v3.json') == 'd14acb535b93f9d590caa47ed7ee78d1cbd9ed914b8d33741ebd1dce03a30ea4'

receipt_specs = [
    ('source_policy_microfetch_receipt_v4_2.json', 'source_policy_microfetch_plan_v4_2.json', 'source_policy_microfetch_v4_2.py'),
    ('source_policy_followup_microfetch_receipt_v4_1.json', 'source_policy_followup_microfetch_plan_v4_1.json', 'source_policy_followup_microfetch_v4_1.py'),
    ('source_policy_followup_microfetch_receipt_v4_2.json', 'source_policy_followup_microfetch_plan_v4_2.json', 'source_policy_followup_microfetch_v4_2.py'),
]
responses = {}
diagnostic_bytes = 0
executed_plans = {}
for receipt_name, plan_name, code_name in receipt_specs:
    receipt = load(S / receipt_name)
    plan = load(S / plan_name)
    assert sha(S / plan_name) == receipt['plan_sha256']
    assert sha(S / code_name) == receipt['script_sha256']
    assert len(receipt['responses']) == plan['request_count'] == len(plan['requests'])
    assert plan['workers'] == 1 and plan['automatic_retry'] is False
    assert plan['TLS_verification'] is True and plan['HTTPS_only'] is True
    assert plan['per_response_body_cap_bytes'] <= 1048576 and plan['total_retained_cap_bytes'] <= 20971520
    assert receipt['GWAS_or_reference_body_bytes'] == 0
    assert receipt['form_submission_or_investigator_contact'] is False
    assert receipt['copyrighted_fulltext_committed'] is False
    for request, r in zip(plan['requests'], receipt['responses']):
        assert request['id'] == r['id'] and request['url'] == r['url']
        assert r['id'] not in responses
        cache = Path(plan['cache_namespace'])
        path = regular(r['path'])
        assert path.parent == cache and path.name == r['id'] + '.response.txt'
        assert path.stat().st_size == r['bytes'] <= plan['per_response_body_cap_bytes']
        assert sha(path) == r['sha256'] == evidence['response_identity_sha256'][str(path)]
        assert r['metadata_only'] is True and r['GWAS_or_reference_body'] is False
        assert not re.search(r'\.(gz|bgz|zip|vcf|bed|bim|fam)(\?|$)', r['url'])
        diagnostic_bytes += len(r['curl_writeout'].encode()) + len(r['stderr'].encode())
        assert len(r['curl_writeout'].encode()) <= plan['per_response_header_cap_bytes']
        assert len(r['stderr'].encode()) <= plan['per_response_header_cap_bytes']
        responses[r['id']] = r
    assert sum(r['bytes'] for r in receipt['responses']) == receipt['actual_response_body_bytes_retained']
    executed_plans[str(S / plan_name)] = receipt['plan_sha256']

assert len(responses) == 42 and len(evidence['response_identity_sha256']) == 42
response_bytes = sum(r['bytes'] for r in responses.values())
assert response_bytes == 301489
assert response_bytes + diagnostic_bytes < 20971520
failed = [r['id'] for r in responses.values() if r['returncode'] != 0 or r['http_code'] != '200']
assert sorted(failed) == ['bcac_breast', 'edinburgh_mdd', 'ldsc_pinned_license']

facts = table(S / 'source_policy_public_facts_v4_1.tsv')
assert len(facts) == 23 and len({f['policy_id'] for f in facts}) == 23
for fact in facts:
    expected = fact['compact_fact_sha256']
    payload = {k: v for k, v in fact.items() if k != 'compact_fact_sha256'}
    assert hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest() == expected
fact_ids = {f['policy_id'] for f in facts}
rows = table(S / 'source_policy_exact_source_ledger_v4_1.tsv')
assert collections.Counter(r['collection'] for r in rows) == {'core45': 45, 'extension100': 100, 'validation13': 13, 'new_Finn_R13_insomnia': 1}
assert len({(r['collection'], r['index']) for r in rows}) == 159
for r in rows:
    assert sha(r['membership_manifest']) == r['membership_manifest_sha256']
    assert set(r['policy_ids'].split(';')) <= fact_ids
    assert r['project_raw_body_publication_status'] == 'PROHIBITED_UNTIL_SOURCE_RIGHTS_AND_HUMAN_RELEASE_AUTHORITY_CLEARED'
    assert r['body_hash_not_a_rights_grant'] == 'True'

core = load(P / 'reviews/independent_core_pipeline_replay_manifest_v4.json')['rows']
registry = {r['source_id']: r for r in table(ROOT / 'config/public_gwas_sources.tsv')}
core_rows = [r for r in rows if r['collection'] == 'core45']
assert len({r['source_id'] for r in core_rows}) == 43
for expected, row in zip(core, core_rows):
    q = registry[expected['source_id']]
    assert row['source_id'] == expected['source_id'] and row['trait_id'] == expected['trait_id']
    assert row['source_page_url'] == q['source_page_url'] and row['frozen_source_url'] == q['download_url']
    assert int(row['frozen_size_bytes']) == int(q['archive_bytes'])
    assert row['historical_or_registry_body_sha256'] == q['archive_sha256']
    if row['accession']:
        rr = responses[row['accession']]
        assert rr['returncode'] == 0 and rr['http_code'] == '200'
        meta = load(rr['path'])
        assert meta['accession_id'] == row['accession']
        assert row['observed_license_mark'] == meta['terms_of_license']
catalog = {r['accession']: r for r in core_rows if r['accession']}
license_counts = collections.Counter(r['observed_license_mark'] for r in catalog.values())
assert license_counts == {'https://creativecommons.org/publicdomain/zero/1.0/': 6, 'https://www.ebi.ac.uk/about/terms-of-use/': 17, 'Please Refer to ReadMe File': 1}

ext = load(P / 'manifests/extension_raw_acquisition_plan_v4_8.json')['members']
for expected, row in zip(ext, [r for r in rows if r['collection'] == 'extension100']):
    assert row['source_id'] == expected['extension_trait_id'] and row['index'] == str(expected['index'])
    assert row['frozen_source_url'] == expected['url'] and row['historical_or_registry_body_sha256'] == expected['expected_sha256']
    assert row['historical_or_metadata_body_md5'] == expected['expected_md5']
    assert int(row['frozen_size_bytes']) == expected['expected_bytes'] and row['policy_ids'] == 'PAN_UKB'
val = load(P / 'statistical_validation/validation_raw_replay_design_v1.json')['members']
for expected, row in zip(val, [r for r in rows if r['collection'] == 'validation13']):
    h = expected['historical_source_body']
    assert row['source_id'] == expected['source_id'] and row['index'] == str(expected['index'])
    assert row['frozen_source_url'] == h['source'] and row['historical_or_registry_body_sha256'] == h['observed_sha256']
    assert row['historical_or_metadata_body_md5'] == h['observed_md5'] and int(row['frozen_size_bytes']) == h['observed_size_bytes']
    if row['accession']:
        meta = load(responses[row['accession']]['path'])
        assert meta['terms_of_license'] == row['observed_license_mark'] == 'https://www.ebi.ac.uk/about/terms-of-use/'
assert sum(not r['accession'] for r in rows if r['collection'] == 'validation13') == 11
fg = load(P / 'manifests/finngen_preprocessing_plan_v4_3_2.json')
new_fg = next(r for r in rows if r['collection'] == 'new_Finn_R13_insomnia')
assert fg['source_generation'] == '1777989563097164'
assert new_fg['historical_or_registry_body_sha256'] == fg['source_identity']['sha256']
assert new_fg['historical_or_metadata_body_md5'] == 'f16c21acbf8b9ebfeb0f26a19f0ecfe3'
assert new_fg['frozen_source_url'].endswith('?generation=1777989563097164')
assert int(new_fg['frozen_size_bytes']) == 809346932

gaps = table(S / 'source_policy_human_gaps_v4_1.tsv')
crosswalk = table(S / 'source_policy_original8group_crosswalk_v4_1.tsv')
assert [g['gap_id'] for g in gaps] == ['H%02d' % i for i in range(1, 13)]
assert len(crosswalk) == 8 and all(r['original_ledger_preserved'] == 'True' for r in crosswalk)
assert not (ROOT / 'LICENSE').exists()
compiled_code_sha256 = {}
for name in ['source_policy_ledger_builder_v4_1.py', *[spec[2] for spec in receipt_specs]]:
    path = S / name
    compile(regular(path).read_bytes(), str(path), 'exec')
    compiled_code_sha256[str(path)] = sha(path)
after = {path: sha(path) for path in original}
assert after == observed
receipt = {
    'schema': 'public_source_policy_metadata_revalidation_receipt_v4_1',
    'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'verdict': 'PASS_METADATA_IDENTITY_AND_MEMBERSHIP_ONLY',
    'script_sha256': sha(__file__),
    'input_metadata_sha256_before': observed,
    'input_metadata_sha256_after': after,
    'compiled_without_executing_code_sha256': compiled_code_sha256,
    'executed_metadata_plans_sha256': executed_plans,
    'source_membership_counts': dict(collections.Counter(r['collection'] for r in rows)),
    'core_unique_sources': 43,
    'policy_facts': len(facts), 'human_gaps': len(gaps), 'original_group_crosswalk_rows': len(crosswalk),
    'exact_catalog_license_counts': dict(license_counts),
    'metadata_requests_revalidated': len(responses), 'response_bytes_revalidated': response_bytes,
    'retained_curl_status_and_diagnostic_string_bytes': diagnostic_bytes,
    'full_HTTP_headers_not_retained_or_reconstructed': True,
    'three_failed_HTTP_attempts_preserved': failed,
    'maximum_policy_response_bytes': max(r['bytes'] for r in responses.values()),
    'policy_response_sha256': {r['path']: r['sha256'] for r in responses.values()},
    'new_Finn_insomnia_body_identity_from_existing_metadata_only': True,
    'GWAS_or_reference_body_reads_or_downloads': 0,
    'network_requests_or_contacts_or_forms': 0,
    'native_acquisition_audits_repeated': False,
    'legal_permission_or_raw_body_publication_granted': False,
    'output_namespace_exclusive': True,
}
with OUT.open('x') as f:
    json.dump(receipt, f, indent=2)
    f.write('\n')
print(json.dumps({'receipt': str(OUT), 'sha256': sha(OUT), 'verdict': receipt['verdict'], 'membership_rows': len(rows), 'metadata_requests': len(responses)}))
