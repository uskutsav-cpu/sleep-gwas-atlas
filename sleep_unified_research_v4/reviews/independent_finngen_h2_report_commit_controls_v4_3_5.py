#!/usr/bin/env python3
"""Narrow v4 report-pair producer/consumer controls on private metadata."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import resource
import sys
import time

R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_corrected_fixture_builder',R/'independent_finngen_h2_corrected_candidate_controls_v4_3_4.py')
C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)
F=C.F;F.OUT=R/'independent_finngen_h2_report_commit_controls_v4_3_5'
sys.path.insert(0,str(F.S))


def run():
    started=time.monotonic();F.OUT.mkdir(exist_ok=False);rows=[]
    path=F.S/'60_verify_finngen_h2_diagnostic_v4.py';fixed=F.sha(path)
    labels=['healthy','postpersist_input_symlink','postpersist_output_mutation','postpersist_second_write_failure',
        'current_prior_input_drift','current_report_bytes_drift','current_report_pending','current_report_failure',
        'current_native_pending','current_stock_symlink','prior_input_drift_during_consumption']
    for label in labels:
        f=C.fixture(label);out=f['dir']/'diagnostic.json';table=out.with_suffix('.tsv')
        v=F.module('_finn_report_'+label,path);old_regular=v.regular_hashes
        if label in ['postpersist_input_symlink','postpersist_output_mutation']:
            def gate(mapping):
                if str(out) in mapping:
                    if label=='postpersist_input_symlink':C.replace_symlink(f['log'])
                    else:out.write_text(out.read_text()+' ')
                return old_regular(mapping)
            v.regular_hashes=gate
        if label=='postpersist_second_write_failure':table.write_text('PRIVATE PRESERVED COLLISION\n')
        verify_accepted=False;verify_error=None
        try:
            with contextlib.redirect_stdout(io.StringIO()):v.verify(f['plan_path'],f['hash'],out)
            verify_accepted=True
        except BaseException as e:verify_error=type(e).__name__+': '+str(e)
        v.regular_hashes=old_regular
        pending,seal=v.report_terminal_paths(out)
        if label=='current_prior_input_drift':f['derivative'].write_text('DRIFT\n')
        if label=='current_report_bytes_drift':out.write_text(out.read_text()+' ')
        if label=='current_report_pending':pending.write_text('VETO\n')
        if label=='current_report_failure':Path(str(out)+'.failure.json').write_text('VETO\n')
        if label=='current_native_pending':f['pending'].write_text('VETO\n')
        if label=='current_stock_symlink':C.replace_symlink(f['log'])
        if label=='prior_input_drift_during_consumption':
            old_require=v.require_committed;changed=[False]
            def consume(a,b,binding,pair):
                value=old_require(a,b,binding,pair)
                if str(b)==str(seal) and not changed[0]:
                    f['derivative'].write_text('DRIFT DURING REPORT CONSUMPTION\n');changed[0]=True
                return value
            v.require_committed=consume
        accepted=False;error=None;proof=None
        try:proof=v.require_verified_report(f['plan_path'],f['hash'],out);accepted=True
        except BaseException as e:error=type(e).__name__+': '+str(e)
        expected=label=='healthy'
        rows.append(dict(label=label,verify_accepted=verify_accepted,verify_error=verify_error,
            consumer_accepted=accepted,consumer_error=error,expected_consumer_acceptance=expected,
            consumer_control_pass=accepted==expected,report_PENDING_present=pending.exists() or pending.is_symlink(),
            report_seal_present=seal.exists(),consumer_proof=proof,fixture_only=True))
    receipt=dict(schema='independent_finngen_v4_report_commit_controls_v4_3_5',controls=rows,
        reviewed_verifier_sha256=fixed,reviewed_verifier_unchanged=F.sha(path)==fixed,
        all_three_old_provisional_report_failures_vetoed=all(not row['consumer_accepted'] and row['report_PENDING_present'] for row in rows if row['label'].startswith('postpersist_')),
        healthy_committed_report_consumer_pass=rows[0]['consumer_accepted'],
        unexpected_consumer_acceptances=[row['label'] for row in rows if row['consumer_accepted'] and not row['expected_consumer_acceptance']],
        private_regular_file_sha256={str(p):F.sha(p) for p in sorted(F.OUT.rglob('*')) if p.is_file() and not p.is_symlink()},
        no_actual_body_worker_fit_network_mutex=True,
        elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    F.write(R/'independent_finngen_h2_report_commit_controls_receipt_v4_3_5.json',receipt)
    print(json.dumps(dict(controls=len(rows),healthy_consumer=rows[0]['consumer_accepted'],
        three_old_failures_vetoed=receipt['all_three_old_provisional_report_failures_vetoed'],
        unexpected_consumer_acceptances=receipt['unexpected_consumer_acceptances'])))


if __name__=='__main__':run()
