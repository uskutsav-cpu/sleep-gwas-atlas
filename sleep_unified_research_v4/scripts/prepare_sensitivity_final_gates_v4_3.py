#!/usr/bin/env python3
"""Preserve v4_2 review/fault witness; bind additive final-resource and baseline-identity corrections."""
import datetime
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    old=SSD/'sensitivity_operational_plan_v4_2.json';digest='dd35eb5004349d5787446fa1fc5237c4d3fdd972f944e4d0eb0acb59153db3ea'
    if sha(old)!=digest:raise RuntimeError('SEALED_V4_PLAN_CHANGED')
    plan=json.loads(old.read_text())
    for path,value in plan['dependencies_sha256'].items():
        if sha(path)!=value:raise RuntimeError('SEALED_V4_DEPENDENCY_CHANGED: '+path)
    for path in [old,P/'scripts/sensitivity_executor_v4_3.py',P/'scripts/verify_sensitivity_corrections_v4_3.py',Path(__file__),SSD/'proofs/sensitivity_correction_controls_v4_2.json',P/'reviews/independent_sensitivity_preflight_v4_2.md',P/'reviews/independent_sensitivity_preflight_v4_2.json',P/'reviews/independent_sensitivity_post_hash_fault_v4_2.py',P/'reviews/independent_sensitivity_post_hash_fault_receipt_v4_2.json']:
        plan['dependencies_sha256'][str(path)]=sha(path)
    plan['preserved_v4_2_plan_sha256']=digest
    plan['executor']=str(P/'scripts/sensitivity_executor_v4_3.py')
    for job in plan['jobs']+plan['audit_jobs']:
        prior=job['out_prefix'];new=prior.replace('/fits_v4/','/fits_v4_3/').replace('/intersections_v4/','/intersections_v4_3/')
        job['out_prefix']=new;job['ldsc_args'][job['ldsc_args'].index('--out')+1]=new;Path(new).parent.mkdir(parents=True,exist_ok=True)
    plan['control_oracle_addendum']='Additive v4_3: post-output/journal/within-fit/family identity-I/O final resource gates; exact monitor stage/plan and per-job full-precision output binding. All prior source, plans, reviews and fault witnesses retained. No scientific matrix/guard/tolerance change.'
    plan['schema']='frozen_estimator_sensitivity_operational_plan_v4_3'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    out=SSD/'sensitivity_operational_plan_v4_3.json'
    with out.open('x') as f:json.dump(plan,f,indent=2);f.write('\n')
    print(json.dumps(dict(plan=str(out),sha256=sha(out),executor=plan['executor'],real_workers_launched=0),indent=2))


if __name__=='__main__':main()
