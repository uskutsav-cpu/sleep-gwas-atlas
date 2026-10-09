#!/usr/bin/env python3
"""Replay the original bounded-memory prefilter, separately from LDSC fits."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
SSD_ROOT = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS')
PYTHON = SSD_ROOT/'FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python'
ALLOWLIST = SSD_ROOT/'FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr/w_hm3.snplist'
EXPECTED_ALLOWLIST = 'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed'
DEST = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/lipid_prefilter_reproduction_v2')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024**2), b''):
            h.update(block)
    return h.hexdigest()

def save(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')

def main(trait):
    acquisition_path = PACKAGE/'logs'/(trait+'_ranged_acquisition_receipt_v2.json')
    source = json.loads(acquisition_path.read_text())
    if source['trait_id'] != trait or source['status'] != 'EXACT_HISTORICAL_SOURCE_REACQUIRED':
        raise SystemExit('EXACT_SOURCE_ACQUISITION_GATE_FAILED')
    raw = Path(source['path'])
    qc = SSD_ROOT/'Sep1-local-dependencies/data/harmonized'/(trait+'.qc.txt')
    historical = dict(line.split('\t',1) for line in qc.read_text().splitlines() if '\t' in line)
    if source['expected_sha256'] != historical['prefilter_source_sha256'] or sha(raw) != source['expected_sha256']:
        raise SystemExit('RAW_SOURCE_HASH_GATE_FAILED')
    if sha(ALLOWLIST) != EXPECTED_ALLOWLIST or historical['prefilter_allowlist_sha256'] != EXPECTED_ALLOWLIST:
        raise SystemExit('ORIGINAL_ALLOWLIST_HASH_GATE_FAILED')
    script = ROOT/'scripts/21_prefilter_hm3.py'
    current_blob = subprocess.check_output(['git','rev-parse','HEAD:scripts/21_prefilter_hm3.py'],cwd=ROOT,text=True).strip()
    historical_blob = subprocess.check_output(['git','rev-parse','659d01cf:scripts/21_prefilter_hm3.py'],cwd=ROOT,text=True).strip()
    disk_blob = subprocess.check_output(['git','hash-object',str(script)],cwd=ROOT,text=True).strip()
    if current_blob != historical_blob or disk_blob != current_blob:
        raise SystemExit('ORIGINAL_NATIVE_SCRIPT_IDENTITY_GATE_FAILED')
    outdir = DEST/trait
    if outdir.exists():
        raise SystemExit('PRIOR_PREFILTER_DESTINATION_PRESERVED')
    outdir.mkdir(parents=True)
    (outdir/'tmp').mkdir()
    internal_free = shutil.disk_usage(ROOT).free
    if internal_free < 512*1024**2 or shutil.disk_usage(outdir).free < 5*1024**3:
        raise SystemExit('LIGHTWEIGHT_PREFILTER_RESOURCE_GATE_FAILED')
    output = outdir/(trait+'.hm3.tsv.gz')
    provenance = outdir/(trait+'.hm3.provenance.json')
    command = [str(PYTHON),str(script),'--input',str(raw),'--output',str(output),
               '--provenance',str(provenance),'--snp-column','rsID','--allowlist',str(ALLOWLIST),
               '--expected-input-bytes',str(source['expected_bytes']),
               '--expected-input-sha256',source['expected_sha256']]
    bound = {'raw_source': raw, 'allowlist': ALLOWLIST, 'native_script': script,
             'native_python_binary': PYTHON.resolve(), 'historical_qc': qc, 'acquisition_receipt': acquisition_path}
    before = {key:sha(path) for key,path in bound.items()}
    plan = {'frozen_utc': datetime.now(timezone.utc).isoformat(), 'trait_id': trait, 'command': command,
            'inputs_sha256_before': before, 'native_script_git_blob': current_blob,
            'native_script_original_main_git_blob': historical_blob, 'wrapper_sha256': sha(__file__),
            'python_version': subprocess.check_output([str(PYTHON),'--version'],text=True).strip(),
            'workers': 1, 'expected_maximum_retained_output_bytes': 128*1024**2,
            'observed_worker_RSS_stop_bytes': 768*1024**2, 'maximum_seconds': 3600,
            'internal_free_before': internal_free, 'internal_minimum_before': 512*1024**2,
            'internal_emergency_stop_bytes': 128*1024**2, 'TMPDIR': str(outdir/'tmp'),
            'scope': 'bounded native HapMap3 rsID prefilter only; not harmonization, munging or LDSC',
            'full_native_LDSC_3GiB_internal_guard_unchanged': True,
            'original_SSD_and_v1_outputs_modified': False}
    plan_path = PACKAGE/'manifests'/(trait+'_prefilter_execution_plan_v2.json')
    save(plan_path, plan)
    print(json.dumps({'trait':trait,'plan_frozen':str(plan_path),'scope':plan['scope']}),flush=True)
    environment = os.environ.copy()
    environment.update(PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(outdir/'tmp'),
                       OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    log_path = PACKAGE/'logs'/(trait+'_native_prefilter_v2.log')
    began = time.time()
    peak = 0
    stopped = None
    termination_escalated = False
    with log_path.open('x') as log:
        process = subprocess.Popen(command,cwd=ROOT,env=environment,stdout=log,stderr=subprocess.STDOUT)
        worker_output_paths = [output, provenance,
                               Path(str(output)+'.tmp.'+str(process.pid)),
                               Path(str(provenance)+'.tmp.'+str(process.pid))]
        def retained_size():
            n = 0
            for path in worker_output_paths:
                try:
                    n += path.stat().st_size
                except FileNotFoundError:
                    # The native worker atomically renames temporary output.
                    # A completed run is checked again after the worker exits.
                    pass
            return n
        while process.poll() is None:
            rss = subprocess.run(['ps','-p',str(process.pid),'-o','rss='],capture_output=True,text=True).stdout.strip()
            if rss:
                peak = max(peak,int(rss)*1024)
            if peak > plan['observed_worker_RSS_stop_bytes']:
                stopped = 'WORKER_RSS_LIMIT'
            elif shutil.disk_usage(ROOT).free < plan['internal_emergency_stop_bytes']:
                stopped = 'INTERNAL_EMERGENCY_FLOOR'
            elif time.time()-began > plan['maximum_seconds']:
                stopped = 'TIME_LIMIT'
            elif retained_size() > plan['expected_maximum_retained_output_bytes']:
                stopped = 'OUTPUT_SIZE_LIMIT'
            if stopped:
                process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    termination_escalated = True
                    process.kill()
                    process.wait()
                break
            time.sleep(2)
    after = {key:sha(path) for key,path in bound.items()}
    if retained_size() > plan['expected_maximum_retained_output_bytes']:
        stopped = 'OUTPUT_SIZE_LIMIT_AFTER_EXIT'
    elif shutil.disk_usage(ROOT).free < plan['internal_emergency_stop_bytes']:
        stopped = 'INTERNAL_EMERGENCY_FLOOR_AFTER_EXIT'
    elif time.time()-began > plan['maximum_seconds']:
        stopped = 'TIME_LIMIT_AFTER_EXIT'
    numerical = json.loads(provenance.read_text()) if provenance.exists() else None
    success = process.returncode == 0 and stopped is None and before == after and numerical is not None
    comparisons = {} if numerical is None else {
        'source_rows': [int(historical['prefilter_source_rows']),numerical['source_rows']],
        'retained_rows': [int(historical['prefilter_retained_rows']),numerical['retained_rows']],
        'gzip_SHA256': [historical['infile_sha256'],sha(output)],
        'allowlist_SHA256': [historical['prefilter_allowlist_sha256'],numerical['allowlist_sha256']]}
    agrees = success and all(a == b for a,b in comparisons.values())
    result = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'trait_id': trait,
              'elapsed_seconds': time.time()-began, 'exit_code': process.returncode, 'stop_reason': stopped,
              'owned_worker_termination_escalated_to_kill': termination_escalated,
              'worker_peak_observed_RSS_bytes': peak, 'inputs_sha256_before': before, 'inputs_sha256_after': after,
              'worker_output_bytes_final_and_temporary': retained_size(),
              'worker_output_paths': [str(p) for p in worker_output_paths if p.exists()],
              'input_hashes_unchanged': before == after, 'resource_plan_sha256': sha(plan_path),
              'log_sha256': sha(log_path), 'comparison': comparisons,
              'native_prefilter_provenance_path': str(provenance),
              'native_prefilter_provenance_sha256': sha(provenance) if provenance.exists() else None,
              'output_path': str(output), 'output_sha256': sha(output) if output.exists() else None,
              'status': 'NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS' if agrees else 'EXECUTION_FAILED_OR_PREFILTER_DIFFERENCE_PRESERVED',
              'original_SSD_and_v1_outputs_modified': False, 'complete_GWAS_QC_or_estimator_reproduction': False}
    save(PACKAGE/'logs'/(trait+'_native_prefilter_receipt_v2.json'),result)
    print(json.dumps(result),flush=True)
    if not agrees:
        raise SystemExit('PREFILTER_EXECUTION_OR_HISTORICAL_COMPARISON_FAILED')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--trait',required=True,choices=['hdl','ldl','triglycerides'])
    main(parser.parse_args().trait)
