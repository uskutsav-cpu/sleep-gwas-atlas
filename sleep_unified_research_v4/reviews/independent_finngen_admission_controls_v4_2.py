#!/usr/bin/env python3
"""Admission-only synthetic metadata controls. Does not call execute/worker."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
sys.path.insert(0,str(P/'scripts'))
SOURCE=P/'scripts/58_run_finngen_feasibility_stage_v2.py'
RECEIPT=P/'reviews/independent_finngen_admission_controls_receipt_v4_2.json'
TEMP_PARENT=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/reviewer_software_controls')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    started=time.monotonic(); before=sha(SOURCE)
    spec=importlib.util.spec_from_file_location('_independent_finn_admission_v2',SOURCE)
    candidate=importlib.util.module_from_spec(spec);spec.loader.exec_module(candidate)
    rejected=[]
    with tempfile.TemporaryDirectory(prefix='independent-finn-admission-',dir=TEMP_PARENT) as td:
        directory=Path(td)
        original=json.loads((P/'manifests/finngen_preprocessing_plan_v4_2.json').read_text())
        base=directory/'synthetic-baseline.json';base.write_text('{}\n')
        review_paths=[directory/n for n in ['science.md','science.json','operational.md','operational.json']]
        for path in review_paths:path.write_text('{}\n' if path.suffix=='.json' else 'Synthetic software fixture.\n')
        hashes={str(p):sha(p) for p in review_paths}
        admission=directory/'synthetic-admission.json'
        plan=dict(original,dependencies_sha256={},baseline_gate_plan=str(base),baseline_gate_plan_sha256=sha(base),
            admission=str(admission),required_independent_review_paths=list(hashes))
        pp=directory/'synthetic-plan.json';pp.write_text(json.dumps(plan));ph=sha(pp)
        valid=dict(execution_admitted=True,plan_sha256=ph,scope=plan['scope'],
            independent_binding_review_pass=True,resource_plan_review_pass=True,independent_review_sha256=hashes)
        cases=[('empty_review_map',dict(valid,independent_review_sha256={})),
            ('review_map_missing', {k:v for k,v in valid.items() if k!='independent_review_sha256'}),
            ('wrong_review_path_set',dict(valid,independent_review_sha256={str(directory/'wrong'):sha(base)})),
            ('null_review_map',dict(valid,independent_review_sha256=None)),
            ('nonstring_digest',dict(valid,independent_review_sha256={k:(None if i==0 else v) for i,(k,v) in enumerate(hashes.items())})),
            ('wrong_length_digest',dict(valid,independent_review_sha256={k:('abc' if i==0 else v) for i,(k,v) in enumerate(hashes.items())})),
            ('nonhex_64_character_digest',dict(valid,independent_review_sha256={k:('g'*64 if i==0 else v) for i,(k,v) in enumerate(hashes.items())})),
            ('mutated_hash',dict(valid,independent_review_sha256={k:('0'*64 if i==0 else v) for i,(k,v) in enumerate(hashes.items())})),
            ('not_admitted',dict(valid,execution_admitted=False)),
            ('wrong_plan_hash',dict(valid,plan_sha256='0'*64)),
            ('wrong_scope',dict(valid,scope='INDEPENDENT_REPLICATION')),
            ('binding_review_false',dict(valid,independent_binding_review_pass=False)),
            ('resource_review_false',dict(valid,resource_plan_review_pass=False))]
        for label,a in cases:
            admission.write_text(json.dumps(a))
            try:candidate.require_admission(pp,ph)
            except (RuntimeError,KeyError) as error:rejected.append(dict(control=label,error=str(error)))
            else:raise AssertionError('fault admission accepted: '+label)
        admission.write_text(json.dumps(valid))
        got=candidate.require_admission(pp,ph)
        assert got[0]==plan and got[1]==valid and got[2]==sha(admission)
        # Separate tiny CLI/import smoke test under the exact frozen interpreter.
        python=original['jobs'][0]['command_template'][0]
        help_result=subprocess.run([python,'-B',str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py'),'--help'],
            capture_output=True,text=True,timeout=5,check=True)
        assert '--inherited-heavy-lock-fd' in help_result.stdout
    assert sha(SOURCE)==before
    r=dict(status='PASS_ADMISSION_REJECTION_CONTROLS_AND_PINNED_INTERPRETER_IMPORT',
        executor_sha256=before,checker_sha256=sha(__file__),rejected_count=len(rejected),rejections=rejected,
        synthetic_complete_review_binding_accepted=True,pinned_interpreter_CLI_import_pass=True,
        controls_modify_synthetic_plan_metadata_only=True,real_study_admission_created=False,
        actual_preprocessing_or_GWAS_or_native_workers=0,tiny_CLI_helper_count=1,GWAS_body_reads=0,
        source_code_unchanged=True,elapsed_seconds=time.monotonic()-started)
    with RECEIPT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r))

if __name__=='__main__':main()
