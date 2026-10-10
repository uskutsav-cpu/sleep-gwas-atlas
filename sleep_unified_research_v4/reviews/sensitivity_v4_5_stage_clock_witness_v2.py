#!/usr/bin/env python3
"""Additive private metadata witness for the recorded but unchecked stage clock."""
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode=True
R=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4/reviews')
SOURCE=R/'sensitivity_v4_5_prelaunch_controls.py'
s=importlib.util.spec_from_file_location('_stage_clock_private_controls',SOURCE)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
assert m.digest(SOURCE)==json.loads((R/'sensitivity_v4_5_prelaunch_controls_receipt.json').read_text())['control_script_sha256']
m.FIX=R/'sensitivity_v4_5_stage_clock_witness_controls_v2';m.FIX.mkdir(exist_ok=False)
code=inspect.getsource(m.terminal_control)
old="        new(target,{'baseline_gate_receipt_sha256'"
new="""        clock=m.SSD/'receipts_v4/stage_clock_v4_5.json'
        new(clock,{'plan_sha256':ph,'initial_started_epoch':0,'fixture_only_not_scientific_evidence':True})
        roles['clock']=clock
        new(target,{'stage_clock_sha256':digest(clock),'baseline_gate_receipt_sha256'"""
assert code.count(old)==1
code=code.replace(old,new).replace("    expect=case=='success'", "    expect=case in ['success','mutate_clock']")
code=code.replace('    if expect:\n', "    if case=='success':\n")
# Execute only this additive control function in the isolated review module.
exec(compile(code,str(__file__)+':private_function','exec'),m.__dict__)
m.terminal_control('mutate_clock')
d=m.FIX/'terminal_mutate_clock';target=d/'SSD/sensitivities/receipts_v4/sensitivity_execution_receipt_v4_5.json'
clock=d/'SSD/sensitivities/receipts_v4/stage_clock_v4_5.json'
stage=json.loads(target.read_text())
assert stage['stage_clock_sha256']!=m.digest(clock)
assert m.RECORDS[-1]['consumer_accepts'] is True
receipt={'schema':'independent_sensitivity_v4_5_stage_clock_binding_witness',
    'executor_sha256':m.digest(m.CODE),'base_controls_sha256':m.digest(SOURCE),
    'witness_script_sha256':m.digest(__file__),'private_control_function_sha256':hashlib.sha256(code.encode()).hexdigest(),
    'observation':m.RECORDS[-1],'recorded_stage_clock_sha256':stage['stage_clock_sha256'],
    'actual_mutated_clock_sha256':m.digest(clock),'terminal_accepted_clock_mismatch':True,
    'fixture_file_sha256':{str(q):m.digest(q) for q in d.rglob('*') if q.is_file() and not q.is_symlink()},
    'GWAS_body_reads':0,'reference_body_reads':0,'actual_workers':0,'actual_fits':0,'execution_admission_granted':False}
out=R/'sensitivity_v4_5_stage_clock_witness_receipt_v2.json';m.new(out,receipt)
print(json.dumps({'clock_mismatch_consumer_accepted':True,'receipt_sha256':m.digest(out)},indent=2))
