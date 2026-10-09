#!/usr/bin/env python3
"""Sequential native reproduction in new outputs, with explicit resource gates."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/'sleep_unified_research_v1'
DATA=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged')
REF=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr')
PYTHON=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python')
RECOVERED=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs')

def tsv(path):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def execution_dependencies():
    code=ROOT.parent/'ldsc-code'
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=code,capture_output=True,text=True,check=True).stdout.strip()
    if revision!='6c673952cee74bd5c57aef1555a03b1c015399a0':raise RuntimeError('LDSC_REVISION_GATE_FAILED')
    names=['ldsc.py']+[str(p.relative_to(code)) for p in sorted((code/'ldscore').glob('*.py'))]
    if len(names)!=10:raise RuntimeError('LDSC_CODE_CARDINALITY_FAILED')
    for name in names:
        tracked=subprocess.run(['git','show','HEAD:'+name],cwd=code,capture_output=True,check=True).stdout
        if digest(code/name)!=hashlib.sha256(tracked).hexdigest():raise RuntimeError('LDSC_TRACKED_CODE_CHANGED')
    paths=[code/name for name in names]+[PACKAGE/'scripts/native_ldsc_capture.py',PYTHON]
    paths += [REF/f'{chrom}.l2.{suffix}' for chrom in range(1,23) for suffix in ['ldscore.gz','M_5_50']]
    return {str(path):digest(path) for path in paths}

def finite_estimates(kind,values,expected):
    def numeric(v):return isinstance(v,(int,float)) and math.isfinite(v)
    if len(values)!=expected:return False
    if kind=='h2':
        v=values[0]
        return all(numeric(v.get(k)) for k in ['tot','tot_se','intercept','intercept_se']) and v['tot_se']>0
    for v in values:
        if v.get('status')!='NATIVE_ESTIMATE_RETURNED':return False
        if not all(numeric(v.get(k)) for k in ['rg_ratio','rg_se','z','p']):return False
        if not(v['rg_se']>0 and 0<=v['p']<=1):return False
        for prefix in ['hsq1','hsq2','gencov']:
            if not all(numeric(v.get(prefix,{}).get(k)) for k in ['tot','tot_se','intercept','intercept_se']):return False
    return True

def jobs():
    core=tsv(ROOT/'config/analysis_panel.tsv'); sleep=[r for r in core if r['domain']=='sleep']; disease=[r for r in core if r['domain']!='sleep']
    assert len(core)==45 and len(sleep)==12 and len(disease)==33
    out=[]
    for r in core:
        options=[]
        if r['type']=='binary':
            if r['pop_prev_citation'] in ['', 'NA','UNRESOLVED']:raise RuntimeError('Missing liability convention citation')
            prevalence=round(float(r['ncase'])/(float(r['ncase'])+float(r['ncontrol'])),6)
            options=['--samp-prev',str(prevalence),'--pop-prev',r['pop_prev']]
        out.append(dict(job_id='core_h2_'+r['trait_id'],stage='core',kind='h2',inputs=[str(DATA/(r['trait_id']+'.sumstats.gz'))],options=options,estimates=1))
    for r in sleep:
        inputs=[str(DATA/(r['trait_id']+'.sumstats.gz'))]+[str(DATA/(d['trait_id']+'.sumstats.gz')) for d in disease]
        out.append(dict(job_id='core_rg_'+r['trait_id'],stage='core',kind='rg',inputs=inputs,options=[],estimates=33))
    ext=tsv(ROOT/'discovery_extension/config/candidate_traits.tsv'); assert len(ext)==100
    sources=[str(RECOVERED/'discovery_extension/data/munged'/(r['extension_trait_id']+'.sumstats.gz')) for r in ext]
    for r,s in zip(ext,sources):
        out.append(dict(job_id='extension_h2_'+r['extension_trait_id'],stage='extension',kind='h2',inputs=[s],options=[],estimates=1))
    for r in sleep:
        out.append(dict(job_id='extension_rg_'+r['trait_id'],stage='extension',kind='rg',inputs=[str(DATA/(r['trait_id']+'.sumstats.gz'))]+sources,options=[],estimates=100))
    source_ids=sorted({r['replication_source_id'] for r in tsv(ROOT/'discovery_extension/config/replication_manifest.tsv') if r['source_curation_status']=='COMPLETE_BEFORE_RESULTS'})
    assert len(source_ids)==13
    for source in source_ids:
        out.append(dict(job_id='validation_h2_'+source,stage='validation',kind='h2',inputs=[str(RECOVERED/'discovery_extension/data/replication/munged'/(source+'.sumstats.gz'))],options=[],estimates=1))
    rep=tsv(ROOT/'discovery_extension/results/replication/replication_rg_jobs.tsv')
    assert sum(int(r['source_count']) for r in rep)==41
    for r in rep:
        inputs=[str(DATA/(r['sleep_trait']+'.sumstats.gz'))]+[str(RECOVERED/'discovery_extension/data/replication/munged'/(s+'.sumstats.gz')) for s in r['source_ids'].split(',')]
        out.append(dict(job_id='validation_rg_'+r['sleep_trait'],stage='validation',kind='rg',inputs=inputs,options=[],estimates=len(inputs)-1))
    return out

def preflight():
    return {'recorded_utc':datetime.now(timezone.utc).isoformat(),
            'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,
            'ssd_free_bytes':shutil.disk_usage(RECOVERED.parent if RECOVERED.parent.exists() else '/Volumes/Extreme SSD').free,
            'internal_minimum_bytes':3*1024**3,'ssd_minimum_bytes':5*1024**3,
            'memory':subprocess.run(['sysctl','hw.memsize','vm.swapusage'],capture_output=True,text=True).stdout}

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['core','extension','validation','all'],default='core')
    p.add_argument('--execute',action='store_true');p.add_argument('--limit',type=int)
    args=p.parse_args();family=jobs();selected=[j for j in family if args.stage=='all' or j['stage']==args.stage]
    if args.limit:selected=selected[:args.limit]
    plan={'protocol':str(PACKAGE/'FROZEN_NEW_ANALYSIS_PROTOCOL.md'),
          'protocol_sha256':digest(PACKAGE/'FROZEN_NEW_ANALYSIS_PROTOCOL.md'),
          'ldsc_commit':'6c673952cee74bd5c57aef1555a03b1c015399a0','jobs':family,
          'numerical_method_changes':False,'worker_count':1,'blas_threads':1,
          'expected_rg_estimates':{'core':396,'extension':1200,'validation':41}}
    plan_path=PACKAGE/'manifests/native_reproduction_jobs_v1.json'
    if plan_path.exists() and json.loads(plan_path.read_text())!=plan:raise SystemExit('FROZEN_PLAN_CHANGED; require a new version')
    if not plan_path.exists():plan_path.write_text(json.dumps(plan,indent=2)+'\n')
    resource=preflight(); resource['resource_gate_pass']=resource['internal_free_bytes']>=resource['internal_minimum_bytes'] and resource['ssd_free_bytes']>=resource['ssd_minimum_bytes']
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    (PACKAGE/'logs'/f'native_resource_preflight_{stamp}.json').write_text(json.dumps(resource,indent=2)+'\n')
    print(json.dumps({'stage':args.stage,'jobs_selected':len(selected),'execute':args.execute,**resource},indent=2),flush=True)
    if not args.execute:return
    if not resource['resource_gate_pass']:raise SystemExit('RESOURCE_GATE_FAILED; no full native job launched')
    core_hashes=tsv(PACKAGE/'tables/native_input_hash_checks.tsv')
    admitted={r['path']:r['actual_sha256'] for r in core_hashes if r['kind']=='core_munged' and r['actual_sha256']}
    recovery_table=PACKAGE/'tables/archived_extension_recovery.tsv'
    if recovery_table.exists():
        admitted.update({r['path']:r['actual_sha256'] for r in tsv(recovery_table) if r['status']=='EXACT_RECEIPT_HASH_RECOVERED'})
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    for index,j in enumerate(selected,1):
        resource=preflight()
        if resource['internal_free_bytes']<resource['internal_minimum_bytes'] or resource['ssd_free_bytes']<resource['ssd_minimum_bytes']:raise SystemExit('RESOURCE_FLOOR_REACHED_BETWEEN_JOBS')
        for source in j['inputs']:
            path=Path(source)
            if not path.is_file() or source not in admitted or digest(path)!=admitted[source]:raise SystemExit(f'INPUT_NOT_HASH_ADMITTED {source}')
        output=PACKAGE/'native'/f"{j['stage']}_reproduction_v1"/j['job_id'];output.parent.mkdir(parents=True,exist_ok=True)
        receipt_path=Path(str(output)+'.execution_receipt.json')
        command=[str(PYTHON),'-u',str(PACKAGE/'scripts/native_ldsc_capture.py'),'--ldsc-dir',str(ROOT.parent/'ldsc-code'),
                 '--'+j['kind'],','.join(j['inputs']),'--ref-ld-chr',str(REF)+'/',
                 '--w-ld-chr',str(REF)+'/','--print-delete-vals','--out',str(output)]+j['options']
        dependencies=execution_dependencies()
        if receipt_path.exists():
            previous=json.loads(receipt_path.read_text())
            full=Path(str(output)+'.full_precision.json')
            sealed_outputs=previous.get('all_output_sha256',{})
            if (previous.get('returncode')==0 and previous.get('scientific_cardinality_gate_pass') is True
                and previous.get('execution_identity_gate_pass') is True and previous.get('command')==command
                and previous.get('dependency_sha256_before')==dependencies and sealed_outputs
                and all(Path(p).is_file() and digest(Path(p))==h for p,h in sealed_outputs.items())
                and full.is_file() and previous.get('output_sha256')==digest(full)):
                print(f'CHECKPOINT_VERIFIED {index}/{len(selected)} {j["job_id"]}',flush=True);continue
            raise SystemExit('PRIOR_EXECUTION_REQUIRES_EXPLICIT_AUDIT; no overwrite')
        if Path(str(output)+'.log').exists():raise SystemExit('UNSEALED_PREVIOUS_RUN_PRESERVED')
        print(f'NATIVE_START {index}/{len(selected)} {j["job_id"]}',flush=True)
        started=time.time()
        started_utc=datetime.now(timezone.utc).isoformat()
        with Path(str(output)+'.stdout.log').open('x') as log:
            run=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        full=Path(str(output)+'.full_precision.json')
        valid=False
        if run.returncode==0 and full.exists():
            values=json.loads(full.read_text())['estimates']
            valid=finite_estimates(j['kind'],values,j['estimates'])
        inputs_after={s:digest(Path(s)) for s in j['inputs']}
        dependencies_after=execution_dependencies()
        identity=inputs_after=={s:admitted[s] for s in j['inputs']} and dependencies_after==dependencies
        output_hashes={str(p):digest(p) for p in sorted(output.parent.glob(output.name+'*')) if p.is_file()}
        receipt={'job':j,'command':command,'returncode':run.returncode,'scientific_cardinality_gate_pass':valid,
                 'elapsed_seconds':time.time()-started,'started_utc':started_utc,'completed_utc':datetime.now(timezone.utc).isoformat(),
                 'input_sha256':{s:admitted[s] for s in j['inputs']},
                 'input_sha256_after':inputs_after,'dependency_sha256_before':dependencies,
                 'dependency_sha256_after':dependencies_after,'execution_identity_gate_pass':identity,
                 'all_output_sha256':output_hashes,
                 'output_sha256':digest(full) if full.exists() else None,
                 'original_outputs_modified':False,'resource_preflight':resource}
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
        if not(valid and identity):raise SystemExit('NATIVE_JOB_FAILURE_PRESERVED; investigate before any retry')
        print(f'NATIVE_COMPLETE {j["job_id"]} elapsed_seconds={receipt["elapsed_seconds"]:.1f}',flush=True)

if __name__=='__main__':main()
