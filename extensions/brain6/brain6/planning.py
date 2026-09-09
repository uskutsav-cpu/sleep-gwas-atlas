"""Expand already-prepared genome-wide pairs into resumable native PLACO jobs."""
from __future__ import annotations
from pathlib import Path
from .artifacts import verify_artifact
from .io import read_json,write_json,write_tsv,file_record,require,sha256,check_hash
from .pairs import verify_lock
from .pipeline import create_job
from .executor import run_job


def placo_jobs(pair_dir,chunk_dir,pair_lock,pair_id,source_file,source_hash,out,*,reviewed=False):
    """Requires real prepared inputs; no fake data binding or automatic pair choice.

    Creates one global-parameter job now, followed by deterministic chunk job
    materialization via execute_placo_plan after parameters exist. Does not run
    a single statistical calculation until the caller explicitly executes.
    """
    require(reviewed,'Review the new-family PLACO settings before planning')
    require(pair_id!='insomnia__adhd','Existing Pair B is read-only; use import-legacy rather than rerunning')
    lock=verify_lock(pair_lock)
    require(pair_id in {p['pair_id'] for p in lock['pairs']},'Pair not in selected family')
    pair_dir,chunk_dir=Path(pair_dir).resolve(),Path(chunk_dir).resolve()
    pr,cr=verify_artifact(pair_dir),verify_artifact(chunk_dir)
    require(pr['synthetic']==cr['synthetic']==lock['synthetic'],'Synthetic/real mismatch')
    require(any(r['sha256']==sha256(pair_dir/'pair.tsv.gz') for r in cr['inputs']),
            'Chunks were not created from this paired genome-wide input')
    method=lock['native_policy']['pleiotropy']
    require(method=='PLACO_PLUS','This automated new-family planner supports PLACO_PLUS; original PLACO requires a separate uncorrelated-trait review')
    check_hash(source_file,source_hash)
    out=Path(out).resolve();require(not out.exists(),'Choose a new immutable plan directory');out.mkdir(parents=True)
    common={'method':method,'placo_source':str(Path(source_file).resolve()),'extreme_z2':80,
            'p_threshold':1e-4,'absolute_tolerance':1e-13,'seed':20260908}
    settings={**common,'mode':'estimate','pair_file':str(pair_dir/'pair.tsv.gz'),'min_variants':1000000}
    sp=out/'estimate.settings.json';write_json(sp,settings)
    jp=create_job(sp,method='placo',job_id='parameters',inputs={
        'pair':pair_dir/'pair.tsv.gz','placo_source':source_file,'pair_lock':pair_lock},
        outputs={'parameters':'parameters.json','semantic_status':'status.json'},root=out/'jobs',
        reviewed=True,synthetic=lock['synthetic'])
    chunks=read_json(chunk_dir/'chunks.json')
    plan={'schema_version':1,'synthetic':lock['synthetic'],'reviewed':True,'pair_id':pair_id,
          'pair_lock':file_record(pair_lock),'pair':file_record(pair_dir/'pair.tsv.gz'),
          'source':file_record(source_file),'chunks_manifest':file_record(chunk_dir/'chunks.json'),
          'chunk_root':str(chunk_dir),'common':common,'estimate_job':str(jp),'n_pairs':lock['n_pairs'],
          'expected_rows':chunks['total_rows'],'execution_root':str(out/'native_results')}
    from .io import json_hash
    plan['plan_sha256']=json_hash(plan)
    write_json(out/'plan.json',plan)
    return out/'plan.json'


def execute_placo_plan(path,*,allow_synthetic=False):
    plan=read_json(path)
    from .io import json_hash
    digest=plan.pop('plan_sha256',None)
    require(digest==json_hash(plan),'PLACO plan changed after freezing')
    require(not plan['synthetic'] or allow_synthetic,'Synthetic execution requires explicit opt-in')
    for key in ['pair_lock','pair','source','chunks_manifest']:
        check_hash(plan[key]['path'],plan[key]['sha256'])
    lock=verify_lock(plan['pair_lock']['path'])
    require(plan['n_pairs']==lock['n_pairs'],'Family size drift')
    require(plan['pair_id']!='insomnia__adhd','Existing Pair B may not be rerun')
    verify_artifact(plan['chunk_root'])
    root=Path(plan['execution_root']);base=Path(path).resolve().parent
    run_job(plan['estimate_job'],root,allow_synthetic=allow_synthetic)
    parameters=root/'parameters/parameters.json'
    chunks=read_json(plan['chunks_manifest']['path'])['chunks']
    results=[]
    for i,c in enumerate(chunks):
        name=f'chunk_{i:06d}';chunk=Path(plan['chunk_root'])/c['path']
        settings={**plan['common'],'mode':'chunk','chunk_file':str(chunk),'parameters':str(parameters)}
        settings_path=base/'jobs'/(name+'.settings.json')
        if settings_path.exists():require(read_json(settings_path)==settings,'Chunk settings drift')
        else:write_json(settings_path,settings)
        jp=create_job(settings_path,method='placo',job_id=name,
           inputs={'chunk':chunk,'parameters':parameters,'placo_source':plan['source']['path']},
           outputs={'chunk':'chunk.tsv','errors':'numerical_errors.json','semantic_status':'status.json'},
           root=base/'jobs',reviewed=True,synthetic=plan['synthetic'])
        run_job(jp,root,allow_synthetic=allow_synthetic)
        artifact=root/name/'chunk.tsv'
        results.append({'path':str(artifact),'sha256':sha256(artifact)})
    manifest=base/'completed_chunks.tsv'
    if manifest.exists():
        from .io import read_tsv
        require(list(read_tsv(manifest))==results,'Completed chunk manifest drift')
    else:write_tsv(manifest,['path','sha256'],results)
    from .pleiotropy import collate_placo
    if (root/'collated').exists():return verify_artifact(root/'collated')
    result=collate_placo(manifest,root,'collated',expected_rows=plan['expected_rows'],
                         n_pairs=plan['n_pairs'],synthetic=plan['synthetic'])
    return verify_artifact(result)
