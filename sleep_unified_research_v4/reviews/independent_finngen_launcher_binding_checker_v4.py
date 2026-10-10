#!/usr/bin/env python3
"""Small-buffer metadata/code binding review; never hash or open a GWAS body."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import shutil
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--expected-plan-sha256',required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    assert not args.receipt.exists()
    assert sha(args.plan)==args.expected_plan_sha256
    plan=json.loads(args.plan.read_text())
    assert plan['stage']=='preprocessing' and plan['new_rg_commands']==0
    source=Path(plan['source'])
    assert not source.is_symlink() and source.is_file()
    assert source.stat().st_size==plan['source_identity']['bytes']==809346932
    assert plan['source_identity']==dict(bytes=809346932,md5='f16c21acbf8b9ebfeb0f26a19f0ecfe3',sha256='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049')
    assert plan['source_generation']=='1777989563097164'
    assert plan['input_sha256']=={str(source):plan['source_identity']['sha256']}
    deps=plan['dependencies_sha256']
    assert str(source) not in deps
    # Fail closed if later plans unexpectedly add a body to the dependency set.
    sizes={path:Path(path).stat().st_size for path in deps}
    assert max(sizes.values()) < (32<<20) and sum(sizes.values()) < (128<<20)
    initial=shutil.disk_usage('/System/Volumes/Data').free
    assert initial > (512<<20)
    before={path:sha(path) for path in deps}
    assert before==deps
    base=Path(plan['baseline_gate_plan'])
    assert sha(base)==plan['baseline_gate_plan_sha256']=='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
    previous=json.loads(base.read_text())
    assert all(deps.get(k)==v for k,v in previous['dependencies_sha256'].items())
    assert plan['worker_code'] in deps and deps[plan['worker_code']]=='dad16c6de646ce4090d7d753cf921872bd85606e96ff54229a15748a6e4ceee6'
    job=plan['jobs'][0]
    assert len(plan['jobs'])==1
    interpreter=job['command_template'][0]
    interpreter_bound=interpreter in deps
    assert job['command_template'][2]==str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py')
    for name in ['56_preprocess_finngen_insomnia_feasibility.py','56_preprocess_finngen_insomnia_feasibility_v2.py',
                 '57_verify_finngen_row_controls.py','57_verify_finngen_row_controls_v2.py']:
        assert str(P/'scripts'/name) in deps
    after={path:sha(path) for path in deps}
    assert after==before==deps and sha(args.plan)==args.expected_plan_sha256
    final=shutil.disk_usage('/System/Volumes/Data').free
    assert final > (128<<20)
    usage=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert usage < (128<<20) # macOS reports bytes
    receipt=dict(status='PASS_BINDINGS_WITH_UNBOUND_INTERPRETER_FAULT' if not interpreter_bound else 'PASS_FROZEN_BINDINGS_ONLY',
        completed_utc=datetime.now(timezone.utc).isoformat(),plan_sha256=args.expected_plan_sha256,
        checker_sha256=sha(__file__),before_after_dependency_hashes_pass=True,dependency_count=len(deps),
        dependency_bytes=sum(sizes.values()),dependencies={k:dict(sha256=v,bytes=sizes[k]) for k,v in deps.items()},
        all_baseline_dependencies_preserved=True,interpreter=interpreter,interpreter_hash_bound=interpreter_bound,
        source_metadata_matches_sealed_identity=True,source_or_derivative_body_reads=0,
        source_body_hash_independently_rerun=False,source_generation_independently_requeried=False,
        acquisition_receipt_hash_and_metadata_only=True,admission_exists=Path(plan['admission']).exists(),
        stage_receipt_exists=Path(plan['stage_receipt']).exists(),worker_receipt_exists=Path(job['worker_receipt']).exists(),
        internal_initial_free_bytes=initial,internal_final_free_bytes=final,peak_RSS_bytes=usage,
        read_buffer_bytes=65536,elapsed_seconds=time.monotonic()-started,
        actual_preprocessing_or_native_workers=0,all190_completion_or_scientific_admission_claimed=False)
    with args.receipt.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='dependencies'}))

if __name__=='__main__':main()
