#!/usr/bin/env python3
"""Narrow additive metadata review; no body reads, mutex or worker.
Prior simulated operational controls are inherited only by whole-code identity.
"""
import ast
import hashlib
import json
from pathlib import Path
import time

P=Path(__file__).resolve().parents[1]
OUT=P/'reviews/independent_finngen_operational_prelaunch_v4_3_2.json'
PIN='9d6bc2046194d030b3c87bd8b3b3708d90d942cc52c213e8b7919aac8ed4ca41'
B={};checks=[]

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(65536),b''):h.update(x)
    return h.hexdigest()

def bind(p,h=None):
    actual=sha(p);assert h is None or actual==h,(p,actual,h);B[str(p)]=actual;return actual

def check(name,ok,**detail):
    checks.append(dict(name=name,pass_control=bool(ok),**detail));assert ok,checks[-1]

def main():
    t=time.monotonic();assert not OUT.exists()
    pf=P/'manifests/finngen_preprocessing_plan_v4_3_2.json';bind(pf,PIN);plan=json.loads(pf.read_text())
    oldf=P/'manifests/finngen_preprocessing_plan_v4_3_1.json';bind(oldf,'8e6ffda56823fdbf4f1a7b64d854072383f2732d69c9a02dd465fcb6aee6ed06');old=json.loads(oldf.read_text())
    sealfile=P/'reviews/independent_finngen_operational_prelaunch_seal_v4_3_1.json';bind(sealfile,'560c8955d3117b13da171147dabf5fc41fadb39580bc6b7cef0488f1e996df9e');seal=json.loads(sealfile.read_text())
    for p,h in seal['file_sha256'].items():bind(p,h)
    prior=json.loads((P/'reviews/independent_finngen_operational_prelaunch_v4_3_1.json').read_text())
    check('blocked_predecessor_preserved',prior['verdict']=='PRELAUNCH_BLOCKED_UNBOUND_CHILD_HELPER'
          and prior['control_count']==51 and prior['all_control_expectations_pass'])
    helper=P/'scripts/canonical_calibration_common_v4_3.py';bind(helper,'8b9e0a4d4b7a9e9c8c487b8c5991f5fa232c4171429951183b8dc5d84aff80eb')
    check('actual_child_helper_binding_corrected',plan['dependencies_sha256'].get(str(helper))==B[str(helper)]
          and str(helper) not in old['dependencies_sha256'])
    oldcode=P/'scripts/58_run_finngen_feasibility_stage_v3_1.py';code=P/'scripts/58_run_finngen_feasibility_stage_v3_2.py'
    text=code.read_text().replace('pipeline_replay_v3_2','pipeline_replay_v3_1')
    check('whole_controller_text_and_AST_unchanged_except_private_namespace',text==oldcode.read_text()
          and ast.dump(ast.parse(text),include_attributes=False)==ast.dump(ast.parse(oldcode.read_text()),include_attributes=False))
    def normalize(x):
        if isinstance(x,dict):return {k:normalize(v) for k,v in x.items()}
        if isinstance(x,list):return [normalize(v) for v in x]
        if isinstance(x,str):
            for a,b in [('pipeline_replay_v3_2','pipeline_replay_v3_1'),('_v4_3_2','_v4_3_1'),('_v3_2','_v3_1')]:x=x.replace(a,b)
        return x
    exclude=['prepared_utc','dependencies_sha256','resource_preflight']
    check('all_scientific_operational_fields_and_commands_preserved',normalize({k:v for k,v in plan.items() if k not in exclude})==
          {k:v for k,v in old.items() if k not in exclude})
    check('all_predecessor_dependencies_retained_unchanged',all(plan['dependencies_sha256'].get(p)==h for p,h in old['dependencies_sha256'].items()))
    check('exact_baseline3_1_not_accidental3_2',plan['baseline_gate_plan'].endswith('/sensitivity_operational_plan_v4_3_1.json')
          and plan['baseline_gate_plan_sha256']=='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4')
    unread=[]
    for p,h in plan['dependencies_sha256'].items():
        f=Path(p)
        if f.suffix in ['.py','.md','.json','.sha256'] or (f.suffix=='.tsv' and '/config/' in p) or (f.name=='python' and '/.ldsc-env/bin/' in p):bind(f,h)
        else:unread.append(dict(path=p,expected_sha256=h))
    queue=[code,P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py',P/'scripts/sensitivity_executor_v4_4.py',
           P/'scripts/30_prepare_and_run_ssd_native_campaign.py'];seen=set();edges=[]
    while queue:
        f=queue.pop()
        if f in seen:continue
        seen.add(f)
        for node in ast.walk(ast.parse(f.read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [a.name for a in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                if name:
                    target=P/'scripts'/(name.split('.')[0]+'.py')
                    if target.exists():
                        edges.append(dict(caller=str(f),dependency=str(target),bound=str(target) in plan['dependencies_sha256']));queue.append(target)
    check('complete_inspected_nine_file_local_import_closure',len(seen)==9 and all(e['bound'] for e in edges),edges=edges)
    check('only_preprocessing_no_eligibility_or_replication_admission',plan['stage']=='preprocessing' and plan['new_rg_commands']==0
          and len(plan['jobs'])==1 and not plan['scientific_source_admitted'] and not plan['independent_replication_established'])
    check('successor_admission_terminal_and_worker_absent',all(not Path(plan[k]).exists() for k in ['admission','stage_receipt','terminal_pending','terminal_seal'])
          and not Path(plan['jobs'][0]['worker_receipt']).exists())
    for p,h in B.items():assert sha(p)==h,('mutation',p)
    result=dict(schema='independent_finngen_operational_successor_binding_v4_3_2',verdict='QUALIFIED_NARROW_OPERATIONAL_PRELAUNCH_PASS',
                plan_path=str(pf),exact_plan_sha256=PIN,checks=checks,check_count=len(checks),all_checks_pass=True,
                code_metadata_before_after_sha256=B,unread_real_body_dependencies=unread,prior_blocked_review_seal_sha256=B[str(sealfile)],
                inherited_controller_controls=28,inherited_control_expectations=51,operational_controls_rerun=False,
                corrected_material_blocker='Exact actual child common3 SHA is now bound; prior blocked3_1 is retained.',
                actual_workers=0,actual_mutex_operations=0,production_body_reads=0,real_decompressions=0,fits=0,
                execution_admission_granted=False,scientific_source_admission=False,replication_established=False,
                qualification='Private immutable marker/deferred-signal contract; absent INFO/per-variant N; assumed N_eff; clinical/rights/cohort gates unresolved. Original exact seven-review admission/full190/runtime resource gates remain required.',
                elapsed_seconds=time.monotonic()-t)
    with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(dict(verdict=result['verdict'],checks=len(checks),receipt_sha256=sha(OUT),elapsed=result['elapsed_seconds'])))

if __name__=='__main__':main()
