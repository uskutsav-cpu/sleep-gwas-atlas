#!/usr/bin/env python3
"""Two narrow corrected-consumer controls and prospective metadata closure."""
import ast
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import resource
import sys
import time

R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_report_fixture',R/'independent_finngen_h2_report_commit_controls_v4_3_5.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
C,F=M.C,M.F;F.OUT=R/'independent_finngen_h2_final_consumer_controls_v4_3_6'
sys.path.insert(0,str(F.S))


def run():
    started=time.monotonic();F.OUT.mkdir(exist_ok=False)
    inputs={str(F.S/n):F.sha(F.S/n) for n in ['60_verify_finngen_h2_diagnostic_v5.py',
        '58_run_finngen_feasibility_stage_v3_6.py','59_prepare_finngen_feasibility_stage_v3_6.py']}
    rows=[]
    for label in ['healthy','prior_input_drift_during_consumption']:
        f=C.fixture(label);d=f['dir'];plan=f['plan'];path=d/'exact_3_6_plan.json'
        stage=d/'exact_3_6_stage.json';pending=d/'exact_3_6_PENDING.json';seal=d/'exact_3_6_seal.json'
        plan.update(executor_path=str(F.S/'58_run_finngen_feasibility_stage_v3_6.py'),stage_receipt=str(stage),
                    terminal_pending=str(pending),terminal_seal=str(seal))
        argv=plan['jobs'][0]['command_template'];argv[argv.index('--plan')+1]=str(path)
        F.write(path,plan);ph=F.sha(path)
        cap=json.loads(f['capture'].read_text());cap['plan_sha256']=ph;F.write(f['capture'],cap)
        worker=json.loads(f['worker'].read_text());worker['plan_sha256']=ph
        worker['command']=[ph if x=='{PLAN_SHA256}' else x for x in argv]
        worker['output_sha256'][str(f['capture'])]=F.sha(f['capture']);F.write(f['worker'],worker)
        old=json.loads(f['stage'].read_text());frozen={p:F.sha(p) for p in old['result_sha256']}
        F.write(stage,dict(status='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED',plan_sha256=ph,
                         admission_sha256='a'*64,owned_cleanup_verified=True,result_sha256=frozen))
        F.commit(pending,seal,dict(plan_sha256=ph,admission_sha256='a'*64,executor_sha256=F.sha(plan['executor_path'])),stage)
        v=F.module('_finn_v5_'+label,F.S/'60_verify_finngen_h2_diagnostic_v5.py');out=d/'diagnostic.json'
        with contextlib.redirect_stdout(io.StringIO()):v.verify(path,ph,out)
        report_pending,report_seal=v.report_terminal_paths(out)
        if label=='prior_input_drift_during_consumption':
            old_require=v.require_committed;changed=[False]
            def consume(a,b,binding,pair):
                value=old_require(a,b,binding,pair)
                if str(b)==str(report_seal) and not changed[0]:
                    f['derivative'].write_text('PRIVATE DRIFT DURING REPORT CONSUMPTION\n');changed[0]=True
                return value
            v.require_committed=consume
        accepted=False;error=None;proof=None
        try:proof=v.require_verified_report(path,ph,out);accepted=True
        except BaseException as e:error=type(e).__name__+': '+str(e)
        rows.append(dict(label=label,consumer_accepted=accepted,expected_acceptance=label=='healthy',
            passed=accepted==(label=='healthy'),error=error,report_producer_commit_pass=True,
            report_PENDING_absent=not report_pending.exists(),consumer_proof=proof,fixture_only=True))
    old=(F.S/'60_verify_finngen_h2_diagnostic_v4.py').read_text()
    new=(F.S/'60_verify_finngen_h2_diagnostic_v5.py').read_text()
    exact=old.replace("report_seal_hash=require_committed(pending,seal,binding,pair)\n",
        "report_seal_hash=require_committed(pending,seal,binding,pair)\n    producer.prior_preprocessing_gate(plan)\n").replace(
        "schema='single_observed_h2_arithmetic_stock_closure_v4'","schema='single_observed_h2_arithmetic_stock_closure_v5'")
    closure=dict(verifier_only_final_prior_repeat_and_schema_change=new==exact,
        runner_only_namespace_change=(F.S/'58_run_finngen_feasibility_stage_v3_6.py').read_text()==(F.S/'58_run_finngen_feasibility_stage_v3_5.py').read_text().replace('observed_h2_v3_5','observed_h2_v3_6'))
    a=ast.parse((F.S/'59_prepare_finngen_feasibility_stage_v3_5.py').read_text());b=ast.parse((F.S/'59_prepare_finngen_feasibility_stage_v3_6.py').read_text())
    def assignments(tree,name):
        return [ast.dump(n.value,include_attributes=False) for n in ast.walk(tree) if isinstance(n,ast.Assign)
                and any(isinstance(t,ast.Name) and t.id==name for t in n.targets)]
    for name in ['derivative','preprocessing_receipt','prior_plan','prior_stage','identity','expected']:
        closure[name+'_AST_unchanged']=assignments(a,name)==assignments(b,name)
    def h2job(tree):
        return [ast.dump(n,include_attributes=False) for n in ast.walk(tree) if isinstance(n,ast.Expr)
            and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute)
            and isinstance(n.value.func.value,ast.Name) and n.value.func.value.id=='job' and n.value.func.attr=='update']
    closure['h2_science_job_update_AST_unchanged']=h2job(a)==h2job(b)
    checks=[]
    for name,expected in [
        ('independent_finngen_h2_candidate_findings_seal_v4_3_3.json','a7e236c85b86d58c40b8612f37a4fd1d32790d7729854ccbbd968873a05a4c8f'),
        ('independent_finngen_h2_corrected_candidate_findings_seal_v4_3_4.json','717aaff9ee2cb7e646dcf4ab492e222388ff2149958efc16d4b6689327f32bc3'),
        ('independent_finngen_h2_report_commit_findings_seal_v4_3_5.json','8542d2fd6d2a3228941317a4c7920e2d1adc212fe23ca2bbce3e1dcad3126c35')]:
        p=R/name;s=json.loads(p.read_text());bad=[k for k,h in s['artifact_sha256'].items() if F.sha(k)!=h]
        assert F.sha(p)==expected and not bad
        checks.append(dict(seal=str(p),seal_sha256=expected,regular_artifacts=len(s['artifact_sha256']),all_current_hashes_match=True))
    receipt=dict(schema='independent_finngen_v5_consumer_prospective3_6_controls',controls=rows,
        all_two_controls_pass=all(x['passed'] for x in rows),prospective_closure_checks=closure,
        all_prospective_closure_checks_pass=all(closure.values()),preserved_old_seals_rechecked=checks,
        reviewed_input_sha256=inputs,reviewed_inputs_unchanged=all(F.sha(k)==v for k,v in inputs.items()),
        private_regular_sha256={str(p):F.sha(p) for p in sorted(F.OUT.rglob('*')) if p.is_file() and not p.is_symlink()},
        no_actual_body_worker_fit_network_mutex=True,actual_plan_prelaunch_review_complete=False,
        elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    F.write(R/'independent_finngen_h2_final_consumer_controls_receipt_v4_3_6.json',receipt)
    print(json.dumps(dict(controls_pass=receipt['all_two_controls_pass'],prospective_closure_pass=receipt['all_prospective_closure_checks_pass'],old_artifacts_rechecked=sum(x['regular_artifacts'] for x in checks))))


if __name__=='__main__':run()
