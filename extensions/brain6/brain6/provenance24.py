"""Verify the transitive artifact DAG, not only the top-level success marker."""
from __future__ import annotations
from pathlib import Path
from .artifacts import verify_artifact, transaction
from .io import read_json, require, check_hash, write_json, write_tsv


def verify_graph(roots, *, synthetic=False, maximum_artifacts=10000):
    todo=[(Path(r['path']).resolve(),r['fingerprint']) for r in roots]
    seen={};files={};records=[]
    while todo:
        path,expected=todo.pop()
        if path in seen:
            require(expected is None or seen[path]==expected,'Same artifact bound to conflicting fingerprints');continue
        require(len(seen)<maximum_artifacts,'Artifact graph exceeds declared traversal bound')
        rr=verify_artifact(path,expected)
        require(rr.get('synthetic')==synthetic,'Synthetic/empirical provenance contamination')
        seen[path]=rr['fingerprint']
        for x in rr['inputs']:
            p=Path(x['path']).resolve();stat=p.stat();cachekey=(str(p),stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
            require(stat.st_size==x['bytes'],'Input size drift in artifact ancestry')
            if cachekey in files:require(files[cachekey]==x['sha256'],'Conflicting input identities')
            else:check_hash(p,x['sha256']);files[cachekey]=x['sha256']
            if p.name=='receipt.json':
                child=read_json(p)
                if isinstance(child,dict) and 'fingerprint' in child and 'outputs' in child:
                    todo.append((p.parent,child['fingerprint']))
        records.append({'path':str(path),'fingerprint':rr['fingerprint'],'stage':rr['stage'],
            'scientific_status':rr.get('scientific_status','NOT_ASSIGNED'),'inputs':len(rr['inputs']),'outputs':len(rr['outputs'])})
    return sorted(records,key=lambda r:r['path'])


def audit(manifest,root,name):
    c=read_json(manifest);require(c.get('artifacts'),'Empty artifact graph')
    records=verify_graph(c['artifacts'],synthetic=c.get('synthetic',False))
    inputs=[manifest]+[Path(r['path'])/'receipt.json' for r in records]
    with transaction(root,name,stage='artifact_provenance_audit',inputs=inputs,parameters=c,synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'artifacts.tsv',list(records[0]),records)
        write_json(work/'status.json',{'status':'PASS','artifacts_verified':len(records),
            'interpretation':'Integrity of artifacts and their inputs, not validity of every scientific model'})
        meta['scientific_status']='PASS'
    return Path(root)/name
