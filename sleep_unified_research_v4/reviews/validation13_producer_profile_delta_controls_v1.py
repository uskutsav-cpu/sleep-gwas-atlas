"""Focused genuine-producer origin/resume and scientific-AST controls only."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'scripts'))
from extension_replay_common_v4 import sha,write_new
from validation_acquisition_producer_profile_v1 import completed_acquisition

path=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/validation_pipeline_replay_v4/validation_pipeline_plan_v4.json')
expected='f8ea94af81b315fca233c556499f45661e15e47b94b6a46b2df429fd19083899'
assert sha(path)==expected
plan=json.loads(path.read_text());contract=plan['acquisition_terminal']
raw=json.loads(Path(contract['producer_plan_path']).read_text())
spec=importlib.util.spec_from_file_location('_producer_profile_controls108',contract['producer_executor_path'])
producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
cases=[]
completed_acquisition(plan);cases.append(dict(case='ACTUAL11_V3_2_V4_TERMINAL13_CONSUMPTION',status='PASS'))


def reject(name,call):
    try:call()
    except (RuntimeError,AssertionError) as e:cases.append(dict(case=name,status='EXPECTED_REJECTION',reason=str(e)))
    else:raise RuntimeError('NEGATIVE_PROFILE_CONTROL_ACCEPTED: '+name)


member=raw['members'][0];pair=contract['source_receipt_pairs'][0]
wrong=copy.deepcopy(raw);wrong['receipt_origins'][member['source_id']]['reused_v3']=False
reject('WRONG_INHERITED_RECEIPT_ORIGIN_REJECTED',lambda:producer.source_receipt_gate(wrong,member,dict(path=pair['primary_path'],sha256=pair['sha256']),contract['binding']['plan_sha256'],verify_body=False))
member=raw['members'][11];pair=contract['source_receipt_pairs'][11]
wrong=copy.deepcopy(raw);wrong['resume_source12']['prefix_bytes']=0
reject('WRONG_SOURCE12_RESUME_OFFSET_REJECTED',lambda:producer.source_receipt_gate(wrong,member,dict(path=pair['primary_path'],sha256=pair['sha256']),contract['binding']['plan_sha256'],verify_body=False))
wrong=copy.deepcopy(plan);wrong['acquisition_terminal']['source_receipt_pairs']=wrong['acquisition_terminal']['source_receipt_pairs'][:-1]
reject('INCOMPLETE12_MEMBER_CONSUMPTION_REJECTED',lambda:completed_acquisition(wrong))
receipt=json.loads(Path(pair['primary_path']).read_text())
assert receipt['observed_headers']['status']==206 and receipt['resume_offset']==57671680
reject('RESUMED206_CANNOT_BE_CREDITED_AS_FRESH200',lambda:producer.assert_headers(member,receipt['observed_headers'],0))
functions=lambda file:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(file.read_text()).body if isinstance(n,ast.FunctionDef)}
old=functions(P/'scripts/96_replay_original_validation_collector_v2.py');new=functions(P/'scripts/validation_collector_producer_profile_v1.py')
assert all(old[k]==new[k] for k in old if k!='completed_acquisition')
spec=importlib.util.spec_from_file_location('_producer_profile_runner_commands',P/'scripts/run_original_validation_pipeline_v4.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
args=SimpleNamespace(plan=path,plan_sha=expected)
commands=[job for member in plan['members'] for job in runner.commands(plan,args,member)]
assert len(commands)==39 and [x[0] for x in commands]==['source','collector','compare']*13
design=json.loads((P/'statistical_validation/validation_raw_replay_design_v1.json').read_text())
assert [m['source_id'] for m in plan['members']]==[m['source_id'] for m in design['members']]
assert all(a['frozen_source_metadata']==b['frozen_source_metadata'] and a['original_source_identity']['sha256']==b['historical_source_body']['observed_sha256'] and a['expected_output_rows']==int(b['original_qc']['output_rows']) and a['effective_N_serialized_12g']==b['effective_N_serialized_12g'] for a,b in zip(plan['members'],design['members']))
assert plan['member_count']==13 and plan['new_estimator_fits']==0 and plan['new_network_transfers']==0
cases.append(dict(case='UNCHANGED96_SCIENTIFIC_ASTS_ORIGINAL13_MEMBERSHIP_N_OUTPUT_ROWS_39_COMMANDS',status='PASS'))
out=P/'reviews/validation13_producer_profile_delta_controls_receipt_v1.json'
write_new(out,dict(status='FOCUSED_AUTHOR_PRODUCER_COMPATIBILITY_CONTROLS_PASS_NO_ADMISSION',cases=cases,case_count=len(cases),
                  plan_path=str(path),plan_sha256=expected,checker_sha256=sha(__file__),
                  body_reads=0,runtime_censuses=0,workers_launched=0,production_mutex_operations=0,old_suite_repeated=False))
print(json.dumps(dict(receipt=str(out),sha256=sha(out),case_count=len(cases)),indent=2))
