#!/usr/bin/env python3
"""Preserve v4 control oracle; bind unique-PID controls without changing executor."""
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
    old=SSD/'sensitivity_operational_plan_v4.json';digest='0b8ea8c6d45e59e47e8cb71ac5668bcd55e9d30960e42b5e55f7c40ed7555000'
    if sha(old)!=digest:raise RuntimeError('SEALED_V4_PLAN_CHANGED')
    plan=json.loads(old.read_text())
    for path,value in plan['dependencies_sha256'].items():
        if sha(path)!=value:raise RuntimeError('SEALED_V4_DEPENDENCY_CHANGED: '+path)
    for path in [old,P/'scripts/verify_sensitivity_corrections_v4_2.py',Path(__file__),SSD/'proofs/sensitivity_v4_control_oracle_failure.json']:
        plan['dependencies_sha256'][str(path)]=sha(path)
    plan['preserved_v4_control_plan_sha256']=digest
    plan['control_oracle_addendum']='Unique PID per mocked software worker; preserved fixed-PID recovery receipt collision. Active deferred-signal executor_v4 unchanged.'
    plan['schema']='frozen_estimator_sensitivity_operational_plan_v4_2'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    out=SSD/'sensitivity_operational_plan_v4_2.json'
    with out.open('x') as f:json.dump(plan,f,indent=2);f.write('\n')
    print(json.dumps(dict(plan=str(out),sha256=sha(out),executor=plan['executor'],real_workers_launched=0),indent=2))


if __name__=='__main__':main()
