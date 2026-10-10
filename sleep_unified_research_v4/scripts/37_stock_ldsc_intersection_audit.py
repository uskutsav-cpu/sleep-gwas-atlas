#!/usr/bin/env python3
"""Replay pinned reads/merges/allele filters with all regression constructors forbidden."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import runpy
import sys


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha',required=True);p.add_argument('--audit-id',required=True);args=p.parse_args()
    if sha(args.plan)!=args.plan_sha:raise ValueError('operational plan hash mismatch')
    plan=json.loads(args.plan.read_text());job=next(j for j in plan['audit_jobs'] if j['audit_id']==args.audit_id)
    for path,expected in plan['dependencies_sha256'].items():
        if sha(path)!=expected:raise ValueError('dependency changed: '+path)
    before={x:sha(x) for x in job['inputs']}
    if any(h!=plan['input_sha256'][x] for x,h in before.items()):raise ValueError('input identity differs')
    out=Path(job['out_prefix']);out.parent.mkdir(parents=True,exist_ok=True)
    if list(out.parent.glob(out.name+'*')):raise FileExistsError('prior audit preserved; no overwrite')
    code=Path(plan['ldsc_dir']);sys.path.insert(0,str(code))
    # Non-main run exposes the exact stock parser but does not enter its estimator dispatch.
    stock=runpy.run_path(str(code/'ldsc.py'),run_name='_stock_parser_only')
    ss=stock['sumstats'];np=stock['np'];pd=stock['pd']
    import scipy
    libraries=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
    if not sys.version.startswith('3.9.23') or libraries!={k:plan['environment'][k] for k in libraries}:raise RuntimeError('pinned environment differs')
    from sensitivity_capture_common import IntersectionCapture,audit_without_fits
    forbidden_calls=[]
    def forbidden(*a,**k):
        forbidden_calls.append(True)
        raise RuntimeError('ESTIMATOR_CONSTRUCTOR_FORBIDDEN_IN_RESULT_FREE_AUDIT')
    for name in ['Hsq','RG','Gencov','LD_Score_Regression']:setattr(ss.reg,name,forbidden)
    class Logger:
        def __init__(self,path):self.f=path.open('x')
        def log(self,msg):self.f.write(str(msg)+'\n');self.f.flush()
    log=Logger(Path(str(out)+'.merge_only.log'))
    capture=IntersectionCapture(ss,np,pd);capture.install()
    parsed=stock['parser'].parse_args(job['ldsc_args'])
    try:identities=audit_without_fits(ss,np,parsed,log,capture)
    finally:log.f.close()
    after={x:sha(x) for x in job['inputs']}
    if before!=after or sha(args.plan)!=args.plan_sha:raise ValueError('input or plan changed during audit')
    if len(identities)!=job['expected_identities'] or forbidden_calls:raise RuntimeError('audit cardinality or forbidden estimator call')
    record=dict(schema='stock_final_intersection_audit_v1',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                audit_id=args.audit_id,arm=job['arm'],plan_sha256=args.plan_sha,arguments=job['ldsc_args'],
                python=sys.version,libraries=libraries,
                source_sha256_before=before,source_sha256_after=after,final_intersections=identities,
                estimator_calls=0,estimator_constructors_forbidden=True,stock_parser_main_dispatch_entered=False,
                dependency_sha256=plan['dependencies_sha256'],
                scope='estimator-free replay of pinned per-pair read/merge/allele and final Z-filter code; not instrumentation of an earlier fit')
    path=Path(str(out)+'.intersection.json')
    with path.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(audit_id=args.audit_id,identities=len(identities),receipt=str(path),sha256=sha(path),estimator_calls=0)),flush=True)


if __name__=='__main__':main()
