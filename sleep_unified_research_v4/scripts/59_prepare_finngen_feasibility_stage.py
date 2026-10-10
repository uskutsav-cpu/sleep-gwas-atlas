#!/usr/bin/env python3
"""Freeze one source-diagnostic stage; no decompression, worker, fit or network."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys

from canonical_calibration_common_v4_3 import check_hashes, sha, utc, write_new

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
NAMESPACE = SSD/'new_source_feasibility/finngen_R13_F5_INSOMNIA'
BASE = SSD/'sensitivities/sensitivity_operational_plan_v4_3_1.json'
BASE_SHA = '05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'


def prepare(stage):
    if sha(BASE) != BASE_SHA:
        raise RuntimeError('REVIEWED_BASELINE_GATE_PLAN_CHANGED')
    baseline = json.loads(BASE.read_text())
    acquisition = P/'logs/finngen_R13_insomnia_acquisition_verified_v4_2.json'
    source = NAMESPACE/'finngen_R13_F5_INSOMNIA.gz'
    acquired = json.loads(acquisition.read_text())
    identity = dict(bytes=809346932, md5='f16c21acbf8b9ebfeb0f26a19f0ecfe3',
                    sha256='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049')
    if acquired['source_path'] != str(source) or acquired['actual_sha256'] != identity['sha256'] or acquired['actual_md5'] != identity['md5'] or acquired['actual_size'] != identity['bytes']:
        raise RuntimeError('ACQUISITION_IDENTITY_DIFFERS')
    if source.is_symlink() or not source.is_file() or source.stat().st_size != identity['bytes'] or sha(source) != identity['sha256']:
        raise RuntimeError('FULL_SOURCE_SHA_OR_SIZE_CHANGED')
    chain = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/hg38ToHg19.over.chain.gz')
    reference = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz')
    alleles = chain.parent/'w_hm3.snplist'
    expected = {str(chain):'14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1',
                str(reference):'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8',
                str(alleles):'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed'}
    check_hashes(expected)
    dependencies = dict(baseline['dependencies_sha256'])
    dependencies.update(expected)
    for path in [BASE, acquisition, ROOT/'scripts/liftover_chain.py',
                 P/'FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md',
                 P/'FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md',
                 P/'scripts/56_preprocess_finngen_insomnia_feasibility.py',
                 P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py',
                 P/'scripts/57_verify_finngen_row_controls.py',
                 P/'scripts/57_verify_finngen_row_controls_v2.py',
                 P/'logs/finngen_row_controls_v4.json',P/'logs/finngen_row_controls_v4_2.json',
                 P/'scripts/58_run_finngen_feasibility_stage.py', Path(__file__),
                 P/'scripts/canonical_calibration_common_v4_3.py']:
        dependencies[str(path)] = sha(path)
    worker_code = P/'scripts/sensitivity_executor_v4_3_1.py'
    if str(worker_code) not in dependencies:
        raise RuntimeError('REVIEWED_WORKER_CODE_NOT_BOUND')
    for directory in ['tmp','cache','logs_v4','receipts_v4','derived','fits']:
        (NAMESPACE/directory).mkdir(parents=True, exist_ok=True)
    derivative = NAMESPACE/'derived/insomnia.sumstats.gz'
    preprocessing_receipt = NAMESPACE/'derived/insomnia.preprocessing.json'
    plan_path = P/'manifests'/('finngen_'+stage+'_plan_v4.json')
    prefix = NAMESPACE/('derived/insomnia' if stage=='preprocessing' else 'fits/insomnia_observed_h2')
    plan = dict(schema='frozen_qualified_finngen_feasibility_stage_v1',prepared_utc=utc(),stage=stage,
                scope='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY',new_rg_commands=0,namespace=str(NAMESPACE),
                source=str(source),source_identity=identity,derivative=str(derivative),
                preprocessing_receipt=str(preprocessing_receipt),source_generation='1777989563097164',
                assumed_effective_N=4.0/(1.0/51643+1.0/446273),
                reference=str(reference),allele_list=str(alleles),chain=str(chain),chain_bytes=1246411,
                chain_sha256=expected[str(chain)],liftover_code=str(ROOT/'scripts/liftover_chain.py'),
                dependencies_sha256=dependencies,input_sha256={str(source):identity['sha256']},
                worker_code=str(worker_code),monitor_code=str(P/'scripts/30_prepare_and_run_ssd_native_campaign.py'),
                baseline_gate_plan=str(BASE),baseline_gate_plan_sha256=BASE_SHA,
                shared_heavy_worker_lock=str(SSD/'native_heavy_worker.lock'),
                admission=str(P/'manifests'/('finngen_'+stage+'_admission_v4.json')),
                stage_receipt=str(P/'logs'/('finngen_'+stage+'_stage_receipt_v4.json')),
                guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,
                           observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=2<<30,
                           deadline_seconds=7200,poll_seconds=2),
                scientific_source_admitted=False, independent_replication_established=False,
                all_190_original_commands_and_numerical_reviews_required=True,
                jobs=[],status='PREPARED_ONLY_NO_WORKER_LAUNCHED')
    job = dict(job_id='finngen_insomnia_'+stage,output_prefix=str(prefix),
               worker_receipt=str(NAMESPACE/'receipts_v4'/(stage+'.execution_receipt.json')))
    if stage == 'preprocessing':
        job['result_receipt'] = str(preprocessing_receipt)
        job['command_template'] = [sys.executable,'-B',str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py'),
                                   '--plan',str(plan_path),'--plan-sha256','{PLAN_SHA256}',
                                   '--inherited-heavy-lock-fd','{HEAVY_LOCK_FD}']
    else:
        prior_plan = P/'manifests/finngen_preprocessing_plan_v4.json'
        prior_stage = P/'logs/finngen_preprocessing_stage_receipt_v4.json'
        prior = json.loads(prior_stage.read_text());result=json.loads(preprocessing_receipt.read_text())
        if prior['status'] != 'QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED' or prior['plan_sha256'] != sha(prior_plan) or result['plan_sha256'] != sha(prior_plan):
            raise RuntimeError('PRIOR_PREPROCESSING_STAGE_NOT_EXACT_SUCCESS')
        check_hashes(prior['result_sha256'])
        if result['status'] != 'QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS' or sha(derivative) != result['derivative_sha256']:
            raise RuntimeError('PREPROCESSING_DERIVATIVE_NOT_VERIFIED')
        for path in [prior_plan,prior_stage,preprocessing_receipt]:plan['dependencies_sha256'][str(path)] = sha(path)
        plan['input_sha256'] = {str(derivative):result['derivative_sha256']}
        plan['python'],plan['ldsc_dir'],plan['environment'] = baseline['python'],baseline['ldsc_dir'],baseline['environment']
        job.update(kind='h2',inputs=[str(derivative)],estimates=1,out_prefix=str(prefix),
                   ldsc_args=['--h2',str(derivative),'--ref-ld-chr',baseline['reference_prefix'],
                              '--w-ld-chr',baseline['reference_prefix'],'--n-blocks','200','--print-delete-vals','--out',str(prefix)])
        job['result_receipt'] = str(prefix)+'.full_precision.json'
        job['command_template'] = [baseline['python'],'-u',str(P/'scripts/38_sensitivity_ldsc_capture.py'),
                                   '--plan',str(plan_path),'--plan-sha','{PLAN_SHA256}','--job-id',job['job_id']]+job['ldsc_args']
    plan['jobs'] = [job]
    plan['resource_preflight'] = dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,ssd_free_bytes=shutil.disk_usage(NAMESPACE).free)
    if plan['resource_preflight']['internal_free_bytes'] < 3<<30 or plan['resource_preflight']['ssd_free_bytes'] < 5<<30:
        raise RuntimeError('FROZEN_RESOURCE_FLOORS_FAILED')
    check_hashes(plan['dependencies_sha256'])
    write_new(plan_path,plan)
    print(json.dumps(dict(stage=stage,plan=str(plan_path),sha256=sha(plan_path),workers_launched=0)))


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['preprocessing','observed_h2'],required=True)
    prepare(p.parse_args().stage)


if __name__ == '__main__':
    main()
