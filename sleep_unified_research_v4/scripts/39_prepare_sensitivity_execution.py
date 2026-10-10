#!/usr/bin/env python3
"""Freeze operational hashes/commands only; does not import LDSC or execute workers."""
import datetime
import hashlib
import json
import shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
MANIFEST=P/'manifests/frozen_estimator_sensitivity_members_v4.json'
MEMBER_SHA='3ce7b50bb4e1289ef5a69d8e46dcef759963703ef7b317e156bdee0d4c40e5e1'


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    if sha(MANIFEST)!=MEMBER_SHA:raise ValueError('scientific member freeze differs')
    m=json.loads(MANIFEST.read_text())
    baseline_path=P/'manifests/ssd_native_execution_plan_v4_3.json'
    baseline=json.loads(baseline_path.read_text())
    proof_path=SSD/'proofs/total_N_derivatives_v1.json'
    proof=json.loads(proof_path.read_text())
    if proof['manifest_sha256']!=MEMBER_SHA or proof['status']!='RESULT_FREE_ONLY_FINITE_N_CHANGED':raise ValueError('source proof not admitted')
    inputs={}
    for job in m['jobs']:
        for source in job['inputs']:
            expected=m['source_hashes'].get(source)
            if expected is None:
                match=next(r for r in proof['results'] if r['derivative']==source)
                if not match['original_and_derivative_invariant_match'] or not match['gzip']['independent_second_compression_matches']:raise ValueError('derivative proof missing')
                expected=match['derivative_sha256']
            if sha(source)!=expected:raise ValueError('source changed: '+source)
            inputs[source]=expected
    for r in proof['results']:
        if sha(r['source'])!=r['source_sha256_before'] or r['source_sha256_before']!=r['source_sha256_after']:raise ValueError('original source changed')
    for directory in ['tmp','cache','logs','receipts','fits','proofs/intersections']:(SSD/directory).mkdir(parents=True,exist_ok=True)
    code=ROOT.parent/'ldsc-code'
    python=next(p for p in baseline['dependencies_sha256'] if p.endswith('/.ldsc-env/bin/python'))
    ref=Path(next(p for p in baseline['dependencies_sha256'] if p.endswith('/1.l2.ldscore.gz'))).parent
    dependencies=dict(baseline['dependencies_sha256'])
    dependencies.update(m['files'])
    dependencies[str(MANIFEST)]=MEMBER_SHA
    dependencies[str(proof_path)]=sha(proof_path)
    dependencies[str(baseline_path)]=sha(baseline_path)
    names=['36_materialize_total_n.py','37_stock_ldsc_intersection_audit.py','38_sensitivity_ldsc_capture.py',
           '39_prepare_sensitivity_execution.py','40_run_estimator_sensitivities.py','sensitivity_capture_common.py',
           'verify_sensitivity_operational_plan.py','30_prepare_and_run_ssd_native_campaign.py']
    for name in names:dependencies[str(P/'scripts'/name)]=sha(P/'scripts'/name)
    for path,expected in dependencies.items():
        if sha(path)!=expected:raise ValueError('frozen dependency hash differs: '+path)
    original_jobs={j['job_id']:j for j in baseline['jobs']}
    jobs=[];audits=[]
    def argv(kind,sources,options,out):
        return ['--'+kind,','.join(sources),'--ref-ld-chr',str(ref)+'/',
                '--w-ld-chr',str(ref)+'/','--n-blocks','200','--print-delete-vals','--out',str(out)]+options
    for job in m['jobs']:
        j=dict(job);j['out_prefix']=str(SSD/'fits'/job['job_id'])
        j['ldsc_args']=argv(j['kind'],j['inputs'],j['options'],j['out_prefix'])
        if j['kind']=='rg':
            sleep=Path(j['inputs'][0]).name.split('.sumstats.gz')[0]
            baseline_id='core_rg_'+sleep
            base_inputs=[next(p for p in m['source_hashes'] if Path(p).name==Path(x).name) for x in j['inputs']]
            base_options=[]
        else:
            trait=Path(j['inputs'][0]).name.split('.sumstats.gz')[0]
            baseline_id='core_h2_'+trait
            base_inputs=original_jobs[baseline_id]['inputs'];base_options=original_jobs[baseline_id]['options']
        if baseline_id not in original_jobs:raise ValueError('missing historical baseline job')
        j['baseline_job_id']=baseline_id
        jobs.append(j)
        for arm,sources,options in [('baseline',base_inputs,base_options),('sensitivity',j['inputs'],j['options'])]:
            aid=arm+'__'+j['job_id'];out=str(SSD/'proofs/intersections'/aid)
            for source in sources:
                if source not in inputs:
                    expected=m['source_hashes'][source]
                    if sha(source)!=expected:raise ValueError('baseline audit source changed')
                    inputs[source]=expected
            audits.append(dict(audit_id=aid,sensitivity_job_id=j['job_id'],baseline_job_id=baseline_id,arm=arm,
                               inputs=sources,expected_identities=j['estimates'],kind=j['kind'],out_prefix=out,
                               ldsc_args=argv(j['kind'],sources,options,out),
                               historical_fit_instrumented=False,scope='exact pinned per-pair read/merge/allele replay; baseline input subset grouped for bounded reads'))
    if len(jobs)!=26 or sum(j['estimates'] for j in jobs)!=62 or len(audits)!=52:raise ValueError('frozen cardinalities differ')
    snapshot=dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,SSD_free_bytes=shutil.disk_usage(SSD).free)
    if snapshot['internal_free_bytes']<3*(1<<30) or snapshot['SSD_free_bytes']<5*(1<<30):raise RuntimeError('original resource guard fails')
    plan=dict(schema='frozen_estimator_sensitivity_operational_plan_v1',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              scientific_member_manifest=str(MANIFEST),scientific_member_manifest_sha256=MEMBER_SHA,
              status='PREPARED_ONLY_NO_AUDITS_OR_FITS_LAUNCHED',ldsc_dir=str(code),python=python,reference_prefix=str(ref)+'/',
              environment=baseline['environment'],input_sha256=inputs,dependencies_sha256=dependencies,
              derivative_proof=str(proof_path),derivative_proof_sha256=sha(proof_path),
              baseline_execution_plan=str(baseline_path),baseline_execution_plan_sha256=sha(baseline_path),
              baseline_dependency_sha256=baseline['dependencies_sha256'],
              baseline_input_sha256={x['path']:x['expected_sha256'] for x in baseline['inputs_verified']},
              baseline_support_package=baseline['ssd_support_package'],original_190_jobs=baseline['jobs'],
              jobs=jobs,audit_jobs=audits,native_command_count=26,new_fitted_estimates=62,
              stock_merge_only_command_count=52,stock_merge_only_intersection_identities=124,
              guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3*(1<<30),SSD_floor_bytes=5*(1<<30),
                         observed_aggregate_worker_RSS_limit_bytes=2*(1<<30),new_output_limit_bytes=4*(1<<30),deadline_seconds=36*3600,poll_seconds=2,
                         shared_heavy_worker_lock=str(SSD.parent/'native_heavy_worker.lock'),
                         shared_lock_operation='fcntl.flock LOCK_EX|LOCK_NB held across all audits/fits and owned-group teardown'),
              execution_order='Verify completed historical 190 jobs and baseline controls; sequential merge-only audits; seal 124 paired intersection identities; then 26 sensitivity commands.',
              resource_preflight=snapshot,estimator_calls_in_preparation=0,pinned_code_edits=0,
              no_biological_testing_family_or_QC_threshold_changes=True)
    out=SSD/'sensitivity_operational_plan_v1.json'
    with out.open('x') as f:json.dump(plan,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(plan=str(out),plan_sha256=sha(out),prepared_native_commands=26,prepared_merge_only_commands=52,new_estimates=62,launched_commands=0),indent=2))


if __name__=='__main__':main()
