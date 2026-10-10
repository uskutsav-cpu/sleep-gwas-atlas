#!/usr/bin/env python3
"""Read-only locked100 identity and static guard audit; no download/body read."""
import ast
import csv
import datetime
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse,parse_qs

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
SEEN = {}


def sha(q):
    h = hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
    return h.hexdigest()


def read(q):
    SEEN[str(q)] = {'sha256':sha(q),'bytes':q.stat().st_size}
    return q.read_text()


def rows(q):
    return list(csv.DictReader(read(q).splitlines(),delimiter='\t'))


def main():
    read(Path(__file__))
    source = P/'scripts/45_acquire_extension_raw_sources.py'
    body = read(source)
    snapshot = P/'reviews/independent_45_original_source_snapshot_v4.txt'
    assert not snapshot.exists();snapshot.write_text(body)
    helper_path = P/'scripts/30_prepare_and_run_ssd_native_campaign.py'
    helper = read(helper_path)
    plan_path = P/'manifests/extension_raw_acquisition_plan_v4.json'
    plan = json.loads(read(plan_path))
    panel = rows(ROOT/'discovery_extension/config/candidate_traits.tsv')
    remote_rows = rows(ROOT/'discovery_extension/provenance/panukbb/remote_object_snapshot.tsv')
    raw = [r for r in remote_rows if r['object_role'] == 'phenotype_sumstats']
    remote = {r['extension_trait_id']:r for r in raw}
    assert len(panel) == len(raw) == len(remote) == len(plan['members']) == 100
    assert len({r['extension_trait_id'] for r in panel}) == len({r['extension_trait_id'] for r in plan['members']}) == 100
    assert len({r['filename'] for r in plan['members']}) == 100
    results = []
    for i,(p,m) in enumerate(zip(panel,plan['members']),1):
        q = remote[p['extension_trait_id']]
        u = urlparse(m['url'])
        checks = {'panel_order_and_index':m['index'] == i and m['extension_trait_id'] == p['extension_trait_id'],
                  'exact_filename':m['filename'] == p['source_filename'] == q['filename'],
                  'original_exact_byte_count':m['expected_bytes'] == int(p['source_file_size_bytes']) == int(q['expected_size_bytes']) == int(q['observed_content_length']),
                  'original_md5':p['checksum'] == q['expected_checksum'] == 'md5:'+m['expected_md5'],
                  'original_S3_version':m['expected_s3_version_id'] == q['s3_version_id'] == parse_qs(u.query)['versionId'][0],
                  'original_versioned_URL':m['url'] == q['versioned_url'],
                  'original_ETag':m['expected_etag'] == q['observed_etag'],
                  'historical_HEAD_status':q['http_status'] == '200' and q['head_verification_status'].startswith('HEAD_SIZE_VERSION'),
                  'build_ancestry_phenotype_release':(m['build'],m['ancestry'],m['phenotype'],m['source_release']) == (p['build'],p['ancestry'],p['phenotype_name'],p['study_accession']),
                  'https_original_bucket':u.scheme == 'https' and u.netloc == 'pan-ukb-us-east-1.s3.amazonaws.com',
                  'body_path_new_namespace':Path(m['body_path']) == Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/extension_raw_replay/raw')/m['filename'],
                  'filename_basename':Path(m['filename']).name == m['filename']}
        results.append({'trait_id':m['extension_trait_id'],'checks':checks,'all_pass':all(checks.values())})
    assert sum(m['expected_bytes'] for m in plan['members']) == plan['compressed_network_bytes'] == 227610388647
    assert max(m['expected_bytes'] for m in plan['members']) == plan['largest_body_bytes'] == 2765036925
    for s,h in plan['bound_sources'].items():
        read(Path(s));assert SEEN[s]['sha256'] == h
    # Evaluate only literal arithmetic assignments from the monitor; no import.
    constants = {}
    for node in ast.parse(helper).body:
        if isinstance(node,ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ['INTERNAL_FLOOR','SSD_FLOOR','RSS_LIMIT','STAGE_SECONDS','OUTPUT_LIMIT']:
            constants[node.targets[0].id] = eval(compile(ast.Expression(node.value),'<literal monitor constant>','eval'),{'__builtins__':{}},{})
    guard_match = {'internal_floor':constants['INTERNAL_FLOOR'] == plan['internal_floor_bytes'],
                   'ssd_floor':constants['SSD_FLOOR'] == plan['ssd_floor_bytes'],
                   'rss':constants['RSS_LIMIT'] == plan['maximum_owned_transfer_rss_bytes'],
                   'largest_body_below_helper_output_limit':plan['largest_body_bytes'] < constants['OUTPUT_LIMIT']}
    unchanged = all(sha(Path(s)) == m['sha256'] for s,m in SEEN.items());assert unchanged
    receipt = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Read-only current original100 plan identity and static operational review; no download,HEAD network call or source-body read.',
               'locked100_identity_all_pass':all(r['all_pass'] for r in results),'members':results,'compressed_network_bytes':plan['compressed_network_bytes'],
               'compressed_network_GiB':plan['compressed_network_bytes']/1024**3,'SSD_budget_GiB':plan['SSD_reservation_bytes']/1024**3,
               'current_helper_constants':constants,'current_helper_guards_match_plan':guard_match,'source_snapshot_path':str(snapshot),'source_snapshot_sha256':sha(snapshot),
               'body_admission_gates_present':['HTTP200','exact Content-Length','S3 versionId','snapshot ETag','original MD5','actual byte count','computed new SHA256'],
               'review_verdict':'GUARD_AND_IMMUTABILITY_AMENDMENT_REQUIRED','findings':[
                   {'id':'PLAN_AND_EXECUTOR_BINDING','severity':'P1','finding':'45 loads PLAN once but hashes the live file again for receipts without checking equality to a frozen initial plan hash. Plan bound_sources omits45 and dynamically imported30. Current metadata matches the historical panel, but mid-run changes are not explicitly rejected or bound as execution dependencies.','remedy':'Version/freeze plan plus executor/helper hashes; compare at launch,before and after every body and completion.'},
                   {'id':'CLEANUP_RECEIPT_FAILURE','severity':'P1','finding':'Per-source finally calls terminate_owned without a cleanup catch. If group teardown raises, control skips per-source receipt writing, unlike corrected native monitor30. A second cleanup attempt can also replace the first proof.','remedy':'Always retain first teardown proof and cleanup error,write per-source FAILED receipt even when teardown fails,then stop family.'},
                   {'id':'PINNED_SIZE_HARD_LIMIT','severity':'P2','finding':'curl lacks exact --max-filesize and version validation; only a two-second sampled size comparison guards the partial. The imported native output cap is4GiB rather than the pinned per-file size.','remedy':'Pin compatible curl and add per-member --max-filesize expected_bytes; retain the sampled/post-exit gates.'},
                   {'id':'BUDGET_CONFIGURATION_BINDING','severity':'P2','finding':'Resource checks import monitor constants rather than explicitly consuming or asserting plan fields. Current3GiB internal,5GiB SSD,2GiB transfer RSS match, but later plan/helper edits are not fail-closed.','remedy':'Hash-bind helper and assert constants against frozen plan or use the frozen plan fields directly.'}],
               'original_source_sha256_available':False,'new_actual_sha_is_not_archival_SHA_identity':True,'gzip_full_stream_verified':False,
               'historical_exact_admission_not_yet_independently_body_verified':True,'scientific_source_quality_or_reproduction_claimed':False,
               'network_or_native_executor_called':False,'GWAS_bodies_opened':False,'all_consumed_hashes_unchanged':unchanged,'inputs':SEEN}
    out = P/'reviews/independent_extension_acquisition_identity_receipt_v4.json';assert not out.exists()
    out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['members','inputs','findings']},indent=2))
    assert receipt['locked100_identity_all_pass'] and all(guard_match.values())


if __name__ == '__main__':
    main()
