#!/usr/bin/env python3
"""Freeze unchanged sensitivity commands with terminal output identity gates."""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v3 import sha, write_new, physical_mount

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'sensitivities'


def main():
    prior=OUT/'sensitivity_operational_plan_v4_3_1.json'
    prior_sha='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
    if sha(prior)!=prior_sha:raise RuntimeError('PRESERVED_SENSITIVITY_PLAN_CHANGED')
    plan=json.loads(prior.read_text())
    for group in ['jobs','audit_jobs']:
        for job in plan[group]:
            for key in ['out_prefix','ldsc_args']:
                old=job[key]
                def route(s):return s.replace('/fits_v4_3/','/fits_v4_5/').replace('/intersections_v4_3/','/intersections_v4_5/')
                job[key]=route(old) if isinstance(old,str) else [route(s) for s in old]
    executor=P/'scripts/sensitivity_executor_v4_5.py'
    ledger=P/'manifests/global_SSD_resource_reservation_v4_3.json'
    if sha(ledger)!='6bd5c7cec9d9880ae630e5bfd431b548bd4cc2f3837d3ce12935fda326a9c590':raise RuntimeError('GLOBAL_RESERVATION_LEDGER_CHANGED')
    files=[prior,executor,Path(__file__),ledger]
    files.extend(P/'scripts'/n for n in ['terminal_commit_common_v2.py','native_stage_completion_v4_3.py'])
    files.extend(P/'reviews'/n for n in ['terminal_commit_common_review_seal_v2.json','genomicsem_native_stage_completion_review_v4_4.md','genomicsem_native_stage_completion_review_v4_4.json','independent_core_numerical_adjudication_v4.sha256','independent_whole_extension_v4.sha256','independent_whole_validation_v4.sha256','independent_sensitivity_preflight_v4_3_1.sha256','sensitivity_final_gate_correction_v4_3_1.sha256'])
    for file in files:plan['dependencies_sha256'][str(file)]=sha(file)
    for file,digest in plan['dependencies_sha256'].items():
        if sha(file)!=digest:raise RuntimeError('SENSITIVITY_FROZEN_DEPENDENCY_CHANGED: '+file)
    plan.update(schema='frozen_estimator_sensitivity_operational_plan_v4_5',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='PREPARED_NO_AUDITS_OR_FITS_LAUNCHED',executor=str(executor),preserved_v4_3_1_plan=str(prior),preserved_v4_3_1_plan_sha256=prior_sha,global_resource_ledger_path=str(ledger),global_resource_ledger_sha256=sha(ledger),reservation_arithmetic=json.loads(ledger.read_text()),terminal_contract='Private immutable owned PENDING marker,deferred catchable signals,read-only identity callback; hash-bound worker outputs and journals,audits,intersection proof,and native captures rechecked after result persistence before final commit.',execution_route='Single full --execute invocation:52 sequential stock merge audits,124 paired identities,then26 unchanged fits. No separate merge-only terminal epoch.')
    plan['guard']['global_reservation_bytes']=300*(1<<30)
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_SENSITIVITY_FLOORS_NOT_MET')
    for folder in ['fits_v4_5','proofs/intersections_v4_5','receipts_v4','logs_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    path=OUT/'sensitivity_operational_plan_v4_5.json';write_new(path,plan)
    print(json.dumps({'plan':str(path),'sha256':sha(path),'audit_commands':52,'paired_identities':124,'fit_commands':26,'estimates':62,'workers_launched':0}))


if __name__=='__main__':main()
