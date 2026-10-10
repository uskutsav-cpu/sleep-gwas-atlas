#!/usr/bin/env python3
"""Metadata-only original13 design: never open a GWAS/reference body or launch a worker."""
import csv
import hashlib
import json
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas')
P = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
ARCHIVE = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
RECOVERED = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension')
HERE = P / 'statistical_validation'
QUEUE = ROOT / 'discovery_extension/results/replication/replication_source_queue.tsv'
META_HASHES = {}


def meta(path):
    """Only explicitly named small JSON/TSV/Python/text metadata, never compressed bodies."""
    path = Path(path)
    if path.suffix not in {'.py', '.json', '.tsv', '.txt'} or path.stat().st_size > 32 << 20:
        raise ValueError('Not a bounded metadata file: ' + str(path))
    body = path.read_bytes()
    META_HASHES[str(path)] = hashlib.sha256(body).hexdigest()
    return body.decode('utf-8')


def tsv(path):
    return list(csv.DictReader(meta(path).splitlines(), delimiter='\t'))


def js(path):
    return json.loads(meta(path))


def stat_only(path):
    path = Path(path)
    try:
        st = path.stat()
        return {'path': str(path), 'exists': True, 'bytes': st.st_size,
                'is_file': path.is_file(), 'sha256_freshly_verified': False}
    except FileNotFoundError:
        return {'path': str(path), 'exists': False, 'sha256_freshly_verified': False}


def main():
    meta(Path(__file__))
    rows = tsv(QUEUE)
    curated = [r for r in rows if r['source_curation_status'] == 'COMPLETE_BEFORE_RESULTS']
    unique = {}
    for r in curated:
        unique.setdefault(r['replication_source_id'], r)
    assert len(rows) == 217 and len(curated) == 58 and len(unique) == 13
    ledger = tsv(ROOT / 'sleep_unified_research_v1/tables/archived_extension_recovery.tsv')
    inventory = tsv(ROOT / 'sleep_unified_research_v1/SSD_INPUT_MANIFEST.tsv')
    input_checks = tsv(ROOT / 'sleep_unified_research_v1/tables/native_input_hash_checks.tsv')
    jobs = tsv(ROOT / 'discovery_extension/results/replication/replication_rg_jobs.tsv')
    result_meta = tsv(ROOT / 'discovery_extension/results/replication/replication_results.tsv')
    h2_meta = tsv(ROOT / 'discovery_extension/results/replication/replication_source_h2.tsv')
    source_code = ROOT / 'discovery_extension/scripts/47_stream_replication_sources.py'
    stream_code = ROOT / 'discovery_extension/scripts/streaming_io.py'
    meta(source_code); meta(stream_code); meta(ROOT / 'environment/tool_versions.tsv')
    head_path = HERE / 'validation_source_head_probes_v1.json'
    heads = {r['source_id']: r for r in js(head_path)['receipts']}
    official_receipt = P / 'source_provenance/sleep_clinician_additional_metadata_receipt_v4.json'
    official = js(official_receipt)
    official_text = {}
    for d in official:
        f = official_receipt.parent / d['cache']
        official_text[d['source_id']] = meta(f)
        assert META_HASHES[str(f)] == d['sha256']
    ref = RECOVERED / 'data/reference/panukbb_hm3_variant_reference.tsv.gz'
    alleles = ARCHIVE / 'work/track_b_completion/local_dependency_copies/ref/w_hm3.snplist'
    fixed_fields = ['replication_study_accession', 'replication_source_url', 'replication_checksum',
                    'replication_source_generation', 'replication_etag', 'replication_content_length_bytes',
                    'replication_storage_mode', 'replication_receipt_path', 'replication_munged_path',
                    'replication_phenotype_definition', 'ancestry', 'build', 'sample_size', 'cases', 'controls',
                    'source_identity_status', 'schema_status', 'effect_allele_status']
    raw_roots = [ROOT/'raw', ROOT/'data/raw', ROOT/'discovery_extension/data/replication/raw',
                 ARCHIVE/'raw', ARCHIVE/'data/raw', ARCHIVE/'discovery_extension/data/replication/raw',
                 Path('/Volumes/Extreme SSD/Codex Archive/2026-08-26-sleep-gwas-atlas/raw'),
                 RECOVERED/'data/replication/raw', SSD/'validation_raw_replay_v1/raw']
    member_rows = []
    members = []
    for index, (sid, r) in enumerate(unique.items(), 1):
        associated = [(i+2, x) for i, x in enumerate(rows)
                      if x['source_curation_status'] == 'COMPLETE_BEFORE_RESULTS' and x['replication_source_id'] == sid]
        for field in fixed_fields:
            assert len({x[field] for _, x in associated}) == 1, (sid, field)
        receipt_path = ROOT / r['replication_receipt_path']
        receipt = js(receipt_path)
        ver = receipt['source_verification']
        assert receipt['pipeline_status'] == 'STREAM_HARMONIZE_DIRECT_MUNGE_PASS'
        assert ver['verification_status'] == 'PASS' and receipt['source_url'] == r['replication_source_url']
        assert ver['expected_md5'] == ver['observed_md5'] == r['replication_checksum'].split(':')[1]
        assert ver['expected_size_bytes'] == ver['observed_size_bytes'] == int(r['replication_content_length_bytes'])
        assert receipt['queue_sha256'] == META_HASHES[str(QUEUE)]
        assert receipt['pipeline_code_sha256'][str(source_code.relative_to(ROOT))] == META_HASHES[str(source_code)]
        assert receipt['pipeline_code_sha256'][str(stream_code.relative_to(ROOT))] == META_HASHES[str(stream_code)]
        qc_path = ROOT / receipt['harmonization_qc']
        qc = {x['metric']: x['value'] for x in tsv(qc_path)}
        assert receipt['harmonization_qc_sha256'] == META_HASHES[str(qc_path)]
        archived = [x for x in ledger if x['archive_member'].endswith('/'+receipt['munged_output'])]
        assert len(archived) == 1
        a = archived[0]
        assert a['status'] == 'EXACT_RECEIPT_HASH_RECOVERED'
        assert a['actual_sha256'] == a['expected_sha256'] == receipt['munged_output_sha256'] == qc['munged_output_sha256']
        archived_stat = stat_only(a['path'])
        assert archived_stat['exists'] and archived_stat['bytes'] == int(a['bytes'])
        cases, controls = int(r['cases']), int(r['controls'])
        n = 4.0 / (1.0/cases + 1.0/controls)
        assert cases + controls == int(r['sample_size'])
        assert format(n, '.12g') == qc['effective_sample_size']
        terminal_metrics = ['not_nonambiguous_hm3', 'duplicate_hm3_rsid', 'extended_MHC',
                            'invalid_beta_or_se', 'invalid_or_p_ci_or_r2_le_0_9',
                            'maf_at_or_below_0_01_or_invalid', 'allele_mismatch_or_ambiguous', 'invalid_z', 'output_rows']
        terminal_sum = sum(int(qc.get(k, 0)) for k in terminal_metrics)
        assert terminal_sum == int(qc['input_rows'])
        filename = Path(urlsplit(r['replication_source_url']).path).name
        candidates = [stat_only(base / name) for base in raw_roots
                      for name in dict.fromkeys([filename, sid+'.gz', sid+'.tsv.gz'])]
        local_meta_matches = [x for x in inventory if x.get('actual_sha256') == ver['observed_sha256']]
        check_matches = [x for x in input_checks if ver['observed_sha256'] in x.values()]
        admitted_local = [x for x in local_meta_matches if x.get('verification_status') == 'MATCH_EXPECTED_SHA256']
        assert not admitted_local and not check_matches
        if sid.startswith('gwas_catalog_'):
            accession = r['replication_study_accession']
            yaml_text = official_text[accession+'_yaml']
            yaml_md5 = next(line.split(':', 1)[1].strip() for line in yaml_text.splitlines() if line.startswith('data_file_md5sum:'))
            assert yaml_md5 == ver['observed_md5']
        else:
            yaml_md5 = None
        m = {'index': index, 'source_id': sid, 'first_original_queue_line': associated[0][0],
             'frozen_source_metadata': {k:r[k] for k in fixed_fields},
             'original_receipt_path': str(receipt_path), 'original_receipt_sha256': META_HASHES[str(receipt_path)],
             'historical_source_body': ver, 'historical_output_sha256': receipt['munged_output_sha256'],
             'archived_processed_evidence': a, 'archived_processed_current_stat_only': archived_stat,
             'original_qc_path': str(qc_path), 'original_qc_sha256': META_HASHES[str(qc_path)],
             'original_qc': qc, 'terminal_qc_partition_sum': terminal_sum,
             'effective_N_independent': n, 'effective_N_serialized_12g': format(n, '.12g'),
             'effective_N_over_total_N': n / (cases+controls),
             'source_reference_sha256': receipt['reference_sha256'],
             'source_allele_reference_sha256': receipt['hm3_alleles_sha256'],
             'candidate_pairs': [{k:x[k] for k in ['pair_id','sleep_trait','extension_trait_id','phenotype_match_status',
                                                 'participant_overlap_status','participant_overlap_evidence','discovery_cohort_relation']}
                                 | {'original_queue_line':i} for i,x in associated],
             'HEAD_probe': heads[sid], 'current_official_yaml_MD5_if_MVP': yaml_md5,
             'local_body_evidence': {'ledger_matches':local_meta_matches,'native_hash_check_matches':check_matches,
                                     'bounded_direct_candidate_path_stats':candidates,
                                     'verdict':'NO_EXACT_HASH_ADMITTED_LOCAL_RAW_BODY_IN_BOUNDED_EVIDENCE'},
             'proposed_body_path':str(SSD/'validation_raw_replay_v1/raw'/filename),
             'proposed_private_workspace':str(SSD/'validation_pipeline_replay_v1/workspace'),
             'proposed_output_path':str(SSD/'validation_pipeline_replay_v1/workspace'/receipt['munged_output'])}
        members.append(m)
        member_rows.append({'index':index,'source_id':sid,'queue_line':associated[0][0],
                            'curated_candidate_pairs':len(associated),'URL':r['replication_source_url'],
                            'generation':r['replication_source_generation'],'historical_ETag':r['replication_etag'],
                            'bytes':ver['observed_size_bytes'],'MD5':ver['observed_md5'],'SHA256':ver['observed_sha256'],
                            'cases':cases,'controls':controls,'N_total':cases+controls,'N_effective_12g':format(n,'.12g'),
                            'raw_input_rows':qc['input_rows'],'processed_rows':qc['output_rows'],
                            'historical_processed_SHA256':a['actual_sha256'],'archived_processed_path':a['path'],
                            'HEAD_current_ETag':heads[sid].get('headers',{}).get('etag',''),
                            'HEAD_verdict':heads[sid]['verdict'],'current_MVP_yaml_MD5_matches':yaml_md5==ver['observed_md5'] if yaml_md5 else 'NA'})
    assert len({m['source_reference_sha256'] for m in members}) == 1
    assert len({m['source_allele_reference_sha256'] for m in members}) == 1
    raw_plan = js(P/'manifests/extension_raw_acquisition_plan_v4_4.json')
    ext_plan_path = SSD/'extension_pipeline_replay_v2/extension_pipeline_replay_plan_v2.json'
    ext_plan = js(ext_plan_path)
    relocation = js(P/'logs/authorized_archive_relocation_receipt_v4.json')
    finn_plan = js(P/'manifests/finngen_preprocessing_plan_v4_2.json')
    meta(P/'scripts/58_run_finngen_feasibility_stage_v2.py')
    meta(Path(finn_plan['worker_code']))
    component = dict(ext_plan['reservation_arithmetic'])
    base = component.pop('total_reserved_component_bytes'); ceiling = component.pop('frozen_ceiling_bytes')
    assert base == sum(component.values()) and ceiling == raw_plan['SSD_reservation_bytes'] == 300 << 30
    component['authorized_relocated_body_bytes'] = relocation['logical_bytes_relocated']
    component['new_Finn_namespace_including_raw_cap_bytes'] = finn_plan['guard']['new_output_limit_bytes']
    component['proposed_core_derivative_cap_bytes'] = 16 << 30
    component['proposed_validation13_raw_bytes'] = sum(m['historical_source_body']['observed_size_bytes'] for m in members)
    component['proposed_validation_pipeline_including_failure_cap_bytes'] = 2 << 30
    total = sum(component.values())
    assert total < ceiling
    job_pairs = [p for j in jobs for p in j['pair_ids'].split(',')]
    assert len(job_pairs) == len(set(job_pairs)) == 41 and len(h2_meta) == 13 and len(result_meta) == 217
    design = {'schema':'original13_metadata_only_proposed_replay_design_v1',
              'status':'PROPOSED_METHOD_RUNTIME_TRANSPORT_QUALIFICATION_REQUIRED_NOT_OPERATIONALLY_FROZEN',
              'recorded_utc':datetime.now(timezone.utc).isoformat(),
              'collector_runtime':{'executable':sys.executable,'python_version':platform.python_version()},
              'scope':{'GWAS_body_bytes_read_or_downloaded':0,'reference_body_bytes_read':0,'decompressions':0,
                       'scientific_workers':0,'estimators':0,'new_pair_outcomes':0},
              'cardinalities':{'candidate_family':217,'curated_candidate_pairs':58,'unique_original_sources':13,
                               'historical_rg_estimates':41,'historical_h2_sources':13,
                               'historical_result_classes':dict(Counter(r['replication_class'] for r in result_meta)),
                               'raw_compressed_bytes':component['proposed_validation13_raw_bytes'],
                               'archived_processed_compressed_bytes':sum(int(m['archived_processed_evidence']['bytes']) for m in members),
                               'processed_rows_expected':sum(int(m['original_qc']['output_rows']) for m in members)},
              'members':members,
              'reference_proposal':{'reference_path':str(ref),'reference_stat_only':stat_only(ref),
                                    'reference_sha256_from_all13_original_receipts':members[0]['source_reference_sha256'],
                                    'allele_list_path':str(alleles),'allele_list_stat_only':stat_only(alleles),
                                    'allele_sha256_from_all13_original_receipts':members[0]['source_allele_reference_sha256'],
                                    'new_coordinate_or_INFO_filters':False,'raw_build_GRCh38_projection':'ORIGINAL_RSID_AND_NONAMBIGUOUS_ALLELE_ONLY'},
              'budget_proposal':{'ceiling_bytes':ceiling,'component_bytes':component,'reserved_total_bytes':total,
                                 'remaining_margin_bytes':ceiling-total,'remaining_margin_GiB':(ceiling-total)/(1<<30),
                                 'Finn2GiB_inclusive_evidence':'58_run_finngen_feasibility_stage_v2.py:109 sets supervisor.SSD to namespace; sensitivity_executor_v4_3_1.py:191 totals every file there',
                                 'proposal_requires_consistent_root_ledger_and_global_actual_namespace_guard':True},
              'metadata_dependencies_sha256':META_HASHES,
              'operational_execution_or_acquisition_admitted':False}
    out = HERE/'validation_raw_replay_design_v1.json'
    with out.open('x') as h:json.dump(design,h,indent=2,sort_keys=True);h.write('\n')
    table = HERE/'validation_raw_replay_design_v1.tsv'
    with table.open('x',newline='') as h:
        writer=csv.DictWriter(h,delimiter='\t',fieldnames=list(member_rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(member_rows)
    print(json.dumps({'JSON_path':str(out),'JSON_SHA256':hashlib.sha256(out.read_bytes()).hexdigest(),
                      'TSV_SHA256':hashlib.sha256(table.read_bytes()).hexdigest(),'cardinalities':design['cardinalities'],
                      'budget':design['budget_proposal'],'metadata_dependencies':len(META_HASHES)},indent=2))


if __name__ == '__main__':
    main()
