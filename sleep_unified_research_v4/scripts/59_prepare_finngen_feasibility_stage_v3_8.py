#!/usr/bin/env python3
"""Freeze one source-diagnostic stage; no decompression, worker, fit or network."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys

from canonical_calibration_common_v4_5 import check_hashes, sha, utc, write_new

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
SOURCE_NAMESPACE = SSD/'new_source_feasibility/finngen_R13_F5_INSOMNIA'
NAMESPACE = SOURCE_NAMESPACE/'pipeline_replay_v3_8'
BASE = SSD/'sensitivities/sensitivity_operational_plan_v4_3_1.json'
BASE_SHA = '05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'


def prepare(stage):
    if sha(BASE) != BASE_SHA:
        raise RuntimeError('REVIEWED_BASELINE_GATE_PLAN_CHANGED')
    if stage!='preprocessing':raise RuntimeError('THIS_VERSION_ADMITS_PREPROCESSING_ONLY')
    baseline = json.loads(BASE.read_text())
    acquisition = P/'logs/finngen_R13_insomnia_acquisition_verified_v4_2.json'
    source = SOURCE_NAMESPACE/'finngen_R13_F5_INSOMNIA.gz'
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
                 P/'scripts/58_run_finngen_feasibility_stage.py', P/'scripts/59_prepare_finngen_feasibility_stage.py',
                 P/'scripts/58_run_finngen_feasibility_stage_v2.py', Path(__file__),
                 P/'manifests/finngen_preprocessing_plan_v4.json',
                 P/'scripts/canonical_calibration_common_v4_5.py']:
        dependencies[str(path)] = sha(path)
    worker_code = P/'scripts/sensitivity_executor_v4_4.py'
    dependencies[str(worker_code)] = sha(worker_code)
    for path in [P/'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py',P/'scripts/57_verify_finngen_row_controls_v3.py',P/'logs/finngen_row_controls_v4_3_7.json',P/'scripts/58_run_finngen_feasibility_stage_v3_2.py',P/'scripts/59_prepare_finngen_feasibility_stage_v3_2.py',P/'manifests/finngen_preprocessing_plan_v4_3_2.json',P/'manifests/finngen_preprocessing_admission_v4_3_2.json',P/'logs/finngen_preprocessing_stage_receipt_v4_3_2.json',P/'logs/finngen_preprocessing_performance_stop_v4_3_2.json',P/'logs/finngen_retry_preparation_harness_correction_v4_3_7.json']:
        dependencies[str(path)]=sha(path)
    stopped=json.loads((P/'logs/finngen_preprocessing_stage_receipt_v4_3_2.json').read_text())
    if stopped['status']!='FAILED_PRESERVED' or stopped['owned_cleanup_verified'] is not True:raise RuntimeError('PRIOR_FAILED_ATTEMPT_MUST_BE_PRESERVED_AND_REAPED')
    stop=json.loads((P/'logs/finngen_preprocessing_performance_stop_v4_3_2.json').read_text())
    for path,digest in stop['preserved_artifact_sha256'].items():
        if sha(path)!=digest:raise RuntimeError('PRIOR_FAILED_ATTEMPT_EVIDENCE_CHANGED')
        dependencies[path]=digest
    for path in [Path(__file__),P/'scripts/58_run_finngen_feasibility_stage_v3_8.py',P/'scripts/terminal_commit_common_v2.py',P/'scripts/native_stage_completion_v4_3.py',P/'scripts/canonical_calibration_common_v4_5.py',P/'scripts/extension_replay_common_v4.py',P/'scripts/canonical_calibration_common_v4_3.py',P/'scripts/58_run_finngen_feasibility_stage_v3_1.py',P/'scripts/59_prepare_finngen_feasibility_stage_v3_1.py',P/'manifests/finngen_preprocessing_plan_v4_3_1.json',P/'scripts/58_run_finngen_feasibility_stage_v3.py',P/'scripts/59_prepare_finngen_feasibility_stage_v3.py',P/'manifests/finngen_preprocessing_plan_v4_3.json',P/'manifests/global_SSD_resource_reservation_v4_4.json',P/'logs/native190_root_independent_completion_addendum_v4.json']:
        dependencies[str(path)]=sha(path)
    for name in ['independent_whole_extension_v4.sha256','independent_whole_validation_v4.sha256','independent_core_numerical_adjudication_v4.sha256','terminal_commit_common_review_seal_v2.json']:
        path=P/'reviews'/name;dependencies[str(path)]=sha(path)
    preserved={'scripts/58_run_finngen_feasibility_stage_v3_7.py':'94183216b0db8e8be4c1d89b6b882b0b7453fdf97e3677fc7130c15e5760ef9c',
        'scripts/59_prepare_finngen_feasibility_stage_v3_7.py':'c8f1e6136257f56c00339e500fc0feeba9e90d323be8167855909e34a7f1b11a',
        'manifests/finngen_preprocessing_plan_v4_3_7.json':'d75fb3f96795f38080081be89da10be04c82d74aa67d11821a7eac0d6e7a8046'}
    for relative,digest in preserved.items():
        path=P/relative
        if sha(path)!=digest:raise RuntimeError('PRESERVED_REJECTED3_7_METADATA_CHANGED')
        dependencies[str(path)]=digest
    rejected=P/'reviews/independent_finngen_operational_prelaunch_seal_v4_3_7.json'
    rejected_sha='8e02e674783178c4765148cac16c00b114554b18c781af51276955c384e15af2'
    if sha(rejected)!=rejected_sha:raise RuntimeError('EXACT_REJECTED3_7_REVIEW_SEAL_REQUIRED')
    record=json.loads(rejected.read_text())
    artifacts=record['artifact_sha256']
    if len(artifacts)!=476:raise RuntimeError('ALL476_REJECTED3_7_ARTIFACTS_REQUIRED')
    check_hashes(artifacts);dependencies.update(artifacts)
    if sha(rejected)!=rejected_sha:raise RuntimeError('REJECTED3_7_SEAL_CHANGED_DURING_CONSUMPTION')
    dependencies[str(rejected)]=rejected_sha
    for directory in ['tmp','cache','logs_v4','receipts_v4','derived','fits']:
        (NAMESPACE/directory).mkdir(parents=True, exist_ok=True)
    derivative = NAMESPACE/'derived/insomnia.sumstats.gz'
    preprocessing_receipt = NAMESPACE/'derived/insomnia.preprocessing.json'
    plan_path = P/'manifests'/('finngen_'+stage+'_plan_v4_3_8.json')
    prefix = NAMESPACE/('derived/insomnia' if stage=='preprocessing' else 'fits/insomnia_observed_h2')
    plan = dict(schema='frozen_qualified_finngen_feasibility_stage_v1',prepared_utc=utc(),stage=stage,
                terminal_pending=str(NAMESPACE/(stage+'_pending_v3_8.json')),terminal_seal=str(NAMESPACE/(stage+'_terminal_seal_v3_8.json')),source_namespace=str(SOURCE_NAMESPACE),global_reservation_bytes=300*(1<<30),
                executor_path=str(P/'scripts/58_run_finngen_feasibility_stage_v3_8.py'),
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
                admission=str(P/'manifests'/('finngen_'+stage+'_admission_v4_3_8.json')),
                stage_receipt=str(P/'logs'/('finngen_'+stage+'_stage_receipt_v4_3_8.json')),
                guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,
                           observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=2<<30,
                           deadline_seconds=7200,poll_seconds=2),
                required_independent_review_paths=[str(P/'reviews'/name) for name in [
                    'independent_finngen_preprocessing_science_review_v4.md',
                    'independent_finngen_row_science_receipt_v4.json',
                    'independent_finngen_operational_preflight_v4_2.md',
                    'independent_finngen_operational_binding_receipt_v4_2.json',
                    'independent_finngen_operational_prelaunch_v4_3_8.md',
                    'independent_finngen_operational_prelaunch_v4_3_8.json',
                    'independent_finngen_operational_prelaunch_seal_v4_3_8.json']],
                scientific_source_admitted=False, independent_replication_established=False,
                all_190_original_commands_and_numerical_reviews_required=True,
                jobs=[],status='PREPARED_ONLY_NO_WORKER_LAUNCHED')
    job = dict(job_id='finngen_insomnia_'+stage,output_prefix=str(prefix),
               worker_receipt=str(NAMESPACE/'receipts_v4'/(stage+'.execution_receipt_v3_8.json')))
    if stage == 'preprocessing':
        job['result_receipt'] = str(preprocessing_receipt)
        job['command_template'] = [baseline['python'],'-B',str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py'),
                                   '--plan',str(plan_path),'--plan-sha256','{PLAN_SHA256}',
                                   '--inherited-heavy-lock-fd','{HEAVY_LOCK_FD}']
    plan['jobs'] = [job]
    plan['resource_preflight'] = dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,ssd_free_bytes=shutil.disk_usage(NAMESPACE).free)
    if plan['resource_preflight']['internal_free_bytes'] < 3<<30 or plan['resource_preflight']['ssd_free_bytes'] < 5<<30:
        raise RuntimeError('FROZEN_RESOURCE_FLOORS_FAILED')
    check_hashes(plan['dependencies_sha256'])
    write_new(plan_path,plan)
    print(json.dumps(dict(stage=stage,plan=str(plan_path),sha256=sha(plan_path),workers_launched=0)))


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['preprocessing'],required=True)
    prepare(p.parse_args().stage)


if __name__ == '__main__':
    main()
