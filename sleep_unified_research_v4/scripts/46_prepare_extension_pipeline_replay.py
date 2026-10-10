#!/usr/bin/env python3
"""Freeze the unchanged 100-source preprocessing replay; never launch workers."""
import csv
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common import sha,write_new,physical_mount

ROOT=Path(__file__).resolve().parents[2];P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v1'
ACQ=P/'manifests/extension_raw_acquisition_plan_v4_4.json'
ACQ_SHA='c318a18bd4912dd15ba3ffcdbd8fe9bf5c97aa535591834a2c065238bd23bc9c'
SENS=SSD/'sensitivities/sensitivity_operational_plan_v4_3_1.json'
SENS_SHA='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'


def main():
    if sha(ACQ)!=ACQ_SHA or sha(SENS)!=SENS_SHA:raise RuntimeError('FROZEN_ACQUISITION_OR_BASELINE_CONTROL_PLAN_CHANGED')
    a=json.loads(ACQ.read_text());b=json.loads(SENS.read_text())
    panel=ROOT/'discovery_extension/config/candidate_traits.tsv'
    with panel.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    traits={r['extension_trait_id']:r for r in rows}
    ledger=ROOT/'sleep_unified_research_v1/tables/archived_extension_recovery.tsv'
    with ledger.open() as f:recovered=[r for r in csv.DictReader(f,delimiter='\t') if '/data/munged/' in r['archive_member'] and r['archive_member'].endswith('.sumstats.gz')]
    archived={Path(r['path']).name.removesuffix('.sumstats.gz'):r for r in recovered}
    if len(rows)!=100 or len(traits)!=100 or len(a['members'])!=100 or set(traits)!=set(archived) or set(traits)!={m['extension_trait_id'] for m in a['members']}:raise RuntimeError('LOCKED_100_MEMBERSHIP_DIFFERS')
    ref=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz')
    w=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/w_hm3.snplist')
    pinned={ROOT/'discovery_extension/scripts/10_harmonize_panukbb.py':'873cbd98f8f41c2bea4d743618ba0363f8c43707b4c1eb6a6280c954ba3f63d6',
            ROOT/'discovery_extension/scripts/streaming_io.py':'28af3379171194ace67c2c09ffe8f3a46e83a54da1875be72a69535fbb566ac8',
            ROOT/'discovery_extension/scripts/11_munge_extension.sh':'598ce18c07cb7d49d64c9064884e1e9e4a0d78961503a704cc7e1e81c376ac97',
            panel:'0a04a3074d798800f47671de39012cbfaf142b9173466e9a309400b7d5553329',
            ROOT/'discovery_extension/config/extension_harmonization_policy.json':'8f3e28277161d3ebba9ee487c3be10d49117264e6b95070386e9db8e7b7f0f7b',
            ref:'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8',w:'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed',
            Path('/usr/sbin/diskutil'):sha('/usr/sbin/diskutil'),Path(b['ldsc_dir'])/'munge_sumstats.py':'65783b746c36ea13a4c1bb64fe210acec9c9cf910a5706081aee3a76ad102be6'}
    deps=dict(b['dependencies_sha256']);deps.update(a['bound_sources'])
    for path,digest in pinned.items():
        if sha(path)!=digest:raise RuntimeError('ORIGINAL_PIPELINE_IDENTITY_DIFFERS: '+str(path))
        deps[str(path)]=digest
    for path in [ACQ,SENS,ledger,ROOT/'sleep_unified_research_v1/tables/extension_processed_source_counts_v1.tsv',ROOT/'sleep_unified_research_v1/logs/extension_processed_source_count_receipt_v1.json',OUT/'METHODS_REVIEW_v1.md',SSD/'sensitivities/proofs/sensitivity_correction_controls_v4_3_1.json',ROOT/'discovery_extension/scripts/08_validate_source_contracts.py',ROOT/'discovery_extension/scripts/05_validate_extension_panel.py',ROOT/'discovery_extension/scripts/09_build_panukbb_hm3_reference.py',ROOT/'discovery_extension/provenance/panukbb/hm3_variant_reference_build.json',ROOT/'scripts/_common.sh',ROOT/'discovery_extension/config/extension_gwas_sources.tsv',ROOT/'discovery_extension/config/extension_gwas_schemas.tsv',ROOT/'discovery_extension/config/extension_variant_mapping_plan.tsv']:
        deps[str(path)]=sha(path)
    ref_build=json.loads((ROOT/'discovery_extension/provenance/panukbb/hm3_variant_reference_build.json').read_text())
    if ref_build['output_sha256']!=pinned[ref] or ref_build['hm3_map_sha256']!='6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4':raise RuntimeError('ORIGINAL_REFERENCE_BUILD_RECEIPT_BINDING_DIFFERS')
    names=['46_prepare_extension_pipeline_replay.py','47_run_extension_pipeline_replay.py','extension_replay_common.py','extension_replay_validate.py','extension_replay_fault_controls.py']
    for name in names:deps[str(P/'scripts'/name)]=sha(P/'scripts'/name)
    for folder in ['harmonized','munged','qc','receipts','logs_v4','receipts_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    members=[];archive_hashes={};python=b['python']
    for m in a['members']:
        t=m['extension_trait_id'];r=json.loads(Path(m['historical_streaming_receipt']).read_text());arc=archived[t]
        if r['pipeline_status']!='STREAM_HARMONIZE_MUNGE_PASS' or r['source_verification']['observed_sha256']!=m['expected_sha256'] or r['reference_sha256']!=pinned[ref]:raise RuntimeError('HISTORICAL_SOURCE_RECEIPT_DIFFERS')
        for rel,h in r['pipeline_code_sha256'].items():
            if str(ROOT/rel) not in deps or deps[str(ROOT/rel)]!=h:raise RuntimeError('HISTORICAL_PIPELINE_CODE_DIFFERS')
        if arc['status']!='EXACT_RECEIPT_HASH_RECOVERED' or arc['expected_sha256']!=arc['actual_sha256'] or arc['expected_sha256']!=r['munged_output_sha256'] or sha(arc['path'])!=arc['expected_sha256']:raise RuntimeError('ARCHIVED_MUNGED_IDENTITY_DIFFERS')
        deps[arc['path']]=arc['expected_sha256'];archive_hashes[arc['path']]=arc['expected_sha256']
        harmonized=str(OUT/'harmonized'/(t+'.txt.gz'));prefix=str(OUT/'munged'/t)
        member=dict(extension_trait_id=t,index=m['index'],acquisition_member=m,acquisition_receipt=str(Path(m['body_path']).parent.parent/'receipts'/(t+'.json')),
                    raw=m['body_path'],raw_sha256=m['expected_sha256'],raw_md5=m['expected_md5'],raw_bytes=m['expected_bytes'],
                    archived_munged=arc['path'],archived_munged_sha256=arc['expected_sha256'],archived_munged_bytes=int(arc['bytes']),
                    harmonized=harmonized,harmonization_qc=str(OUT/'qc'/(t+'.tsv')),harmonization_receipt=str(OUT/'receipts'/(t+'.harmonization.json')),
                    munged=prefix+'.sumstats.gz',munged_prefix=prefix,source_gate_receipt=str(OUT/'receipts'/(t+'.source_gate.json')),comparison_receipt=str(OUT/'receipts'/(t+'.comparison.json')))
        member['harmonize_command']=[python,'-u','-B',str(ROOT/'discovery_extension/scripts/10_harmonize_panukbb.py'),'--trait-id',t,'--source',m['body_path'],'--panel',str(panel),'--reference',str(ref),'--out',harmonized,'--qc-out',member['harmonization_qc'],'--receipt-out',member['harmonization_receipt']]
        member['munge_command']=[python,'-u','-B',str(Path(b['ldsc_dir'])/'munge_sumstats.py'),'--sumstats',harmonized,'--merge-alleles',str(w),'--snp','SNP','--a1','A1','--a2','A2','--frq','FRQ','--p','P','--N-col','N','--signed-sumstats','BETA,0','--chunksize','500000','--out',prefix]
        members.append(member)
    for path,h in deps.items():
        if sha(path)!=h:raise RuntimeError('BOUND_DEPENDENCY_CHANGED: '+path)
    prefix=a['additional_preserved_prefix_bytes'];reserved=a['compressed_network_bytes']+32*(1<<30)+4*(1<<30)+4*(1<<30)+8*(1<<30)+prefix
    if reserved>=a['SSD_reservation_bytes']:raise RuntimeError('PIPELINE_EXCEEDS_FROZEN_RAW_PLUS_DERIVATIVE_RESERVATION')
    qc_roots=[ROOT/'discovery_extension/results/qc/harmonization',ROOT/'sleep_unified_research_v1/sources/recovered/discovery_extension/results/qc/harmonization',Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/discovery_extension/results/qc/harmonization'),Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/results/qc/harmonization')]
    qc_search=[dict(path=str(root),exists=root.exists(),known_panel_QC_files=[str(root/(t+'.tsv')) for t in traits if (root/(t+'.tsv')).is_file()]) for root in qc_roots]
    baseline={k:b[k] for k in ['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs','core_precision_adjudication','python','ldsc_dir','reference_prefix','environment']}
    plan=dict(baseline,schema='historical_extension_raw_pipeline_replay_plan_v1',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='PREPARED_ONLY_NEVER_LAUNCHED',
              members=members,member_count=100,heavy_preprocessing_commands=200,owned_stdlib_validation_commands=200,total_owned_commands=400,estimator_calls=0,
              panel=str(panel),reference=str(ref),reference_sha256=pinned[ref],w_hm3=str(w),template_rows=1217311,
              dependencies_sha256=deps,archived_input_sha256=archive_hashes,acquisition_plan=str(ACQ),acquisition_plan_sha256=ACQ_SHA,
              acquisition_operational_identity=dict(source_folder=str(Path(a['members'][0]['body_path']).parent.parent),curl=a['curl_path'],per_body_seconds_limit=a['per_body_seconds_limit'],resume_prefix_bytes=a['explicit_resume']['prefix_bytes']),
              guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3*(1<<30),SSD_floor_bytes=5*(1<<30),observed_aggregate_worker_RSS_limit_bytes=2*(1<<30),new_output_limit_bytes=32*(1<<30),deadline_seconds=96*3600,per_worker_deadline_seconds=7200,poll_seconds=2,global_reservation_bytes=300*(1<<30),shared_heavy_worker_lock=str(SSD/'native_heavy_worker.lock')),
              reservation_arithmetic=dict(raw_bytes=a['compressed_network_bytes'],pipeline_bytes=32*(1<<30),historical_native_bytes=4*(1<<30),sensitivity_bytes=4*(1<<30),canonical_bytes=8*(1<<30),preserved_prefix_bytes=prefix,total_reserved_component_bytes=reserved,frozen_ceiling_bytes=a['SSD_reservation_bytes']),
              original_harmonization_QC_bounded_search=qc_search,original_filter_count_concordance='UNESTABLISHED_WHEN_QC_TSV_UNAVAILABLE',
              execution_preconditions='Root hash-bound operational admission, exact SHA-bound100 successful acquisition receipts, completed reviewed190 campaign, shared exclusive heavy-worker mutex and unchanged resource guards.',
              compressed_harmonized_or_munged_byte_identity_required=False,full_decompressed_munged_byte_and_field_identity_required=True,automatic_retry=False,scientific_membership_or_threshold_changes=False,raw_GWAS_bodies_read_in_preparation=False)
    mount=physical_mount()
    plan['physical_SSD_preflight']=mount
    state=dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,SSD_free_bytes=shutil.disk_usage(SSD).free)
    if state['internal_free_bytes']<3*(1<<30) or state['SSD_free_bytes']<5*(1<<30):raise RuntimeError('ORIGINAL_RESOURCE_GUARD_FAILS')
    plan['resource_preflight']=state
    out=OUT/'extension_pipeline_replay_plan_v1.json';write_new(out,plan)
    print(json.dumps(dict(plan=str(out),sha256=sha(out),members=100,prepared_owned_commands=400,workers_launched=0,reservation_bytes=reserved),indent=2))


if __name__=='__main__':main()
