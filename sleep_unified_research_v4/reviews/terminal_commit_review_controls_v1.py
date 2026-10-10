#!/usr/bin/env python3
"""Repeat only the two v1 rejection witnesses in fresh private metadata fixtures."""
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
P = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
HELPER = P/'scripts/terminal_commit_common_v1.py'
EXPECTED = '4ed063eb62f5768455cf38957f38ff6b1f32f2d9649757d3c8544f63d67e7e90'
ROOT = P/'reviews/terminal_commit_controls_v1_recheck'

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def module():
    spec = importlib.util.spec_from_file_location('_terminal_v1_private_recheck', HELPER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def main():
    assert digest(HELPER) == EXPECTED
    ROOT.mkdir(exist_ok=False)
    observations = []
    for case in ['broken_addendum', 'empty_stage_binding']:
        d = ROOT/case; d.mkdir()
        m = module(); pending=d/'PENDING'; seal=d/'seal.json'; result=d/'receipt.json'
        result.write_text('{"status":"PROVISIONAL_SUCCESS"}\n')
        binding = {'stage':'metadata_fixture','plan_sha256':'a'*64} if case == 'broken_addendum' else {}
        if case == 'broken_addendum':
            Path(str(result)+'.failure.json').symlink_to(d/'absent-target')
        t=m.TerminalCommit(pending,seal,binding)
        receipts={str(result):digest(result)}
        committed=t.commit(receipts,identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None)
        consumer_seal=m.require_committed(pending,seal,binding,receipts)
        assert committed and consumer_seal == digest(seal) and not pending.exists()
        observations.append({'case':case,'commit_success':committed,'consumer_accepts':True,
                             'false_admission_witness':True,'seal_sha256':consumer_seal})
    # Preserve and independently recheck the original witnesses without writing to them.
    original = P/'reviews/terminal_commit_controls_v1'
    for case in ['broken_addendum','empty_stage_binding']:
        d=original/case; pending=d/'PENDING'; seal=d/'seal.json'; result=d/'receipt.json'
        s=json.loads(seal.read_text())
        accepted=module().require_committed(pending,seal,s['binding'],s['result_receipt_sha256'])
        assert accepted == digest(seal)
    files={}; links={}
    for tree in [original,ROOT]:
        for f in sorted(tree.rglob('*')):
            if f.is_symlink(): links[str(f)]=os.readlink(f)
            elif f.is_file(): files[str(f)]=digest(f)
    assert digest(HELPER) == EXPECTED
    receipt={'schema':'independent_terminal_commit_v1_rejection_controls',
             'helper_sha256':EXPECTED,'control_script_sha256':digest(__file__),
             'observations':observations,'original_witnesses_consumer_rechecked':2,
             'verdict':'REJECT_FALSE_ADMISSION','fixture_regular_file_sha256':files,
             'fixture_symlink_literal_targets':links,'GWAS_body_reads':0,
             'reference_body_reads':0,'decompressions':0,'workers':0,'fits':0,
             'execution_admission_granted':False}
    out=P/'reviews/terminal_commit_controls_receipt_v1.json'
    with out.open('x') as f:
        json.dump(receipt,f,indent=2,sort_keys=True); f.write('\n'); f.flush(); os.fsync(f.fileno())
    print(json.dumps({'receipt_sha256':digest(out),'observations':observations},indent=2))

if __name__ == '__main__': main()
