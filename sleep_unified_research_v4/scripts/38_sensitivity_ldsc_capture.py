#!/usr/bin/env python3
"""Observational return-value/intersection capture around unmodified pinned LDSC."""
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
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--plan-sha',required=True);parser.add_argument('--job-id',required=True)
    known,forwarded=parser.parse_known_args()
    if sha(known.plan)!=known.plan_sha:raise ValueError('operational plan changed')
    plan=json.loads(known.plan.read_text());job=next(j for j in plan['jobs'] if j['job_id']==known.job_id)
    if forwarded!=job['ldsc_args']:raise ValueError('LDSC command differs from frozen operational job')
    for path,expected in plan['dependencies_sha256'].items():
        if sha(path)!=expected:raise ValueError('dependency changed: '+path)
    before={p:sha(p) for p in job['inputs']}
    if any(h!=plan['input_sha256'][p] for p,h in before.items()):raise ValueError('source identity differs')
    code=Path(plan['ldsc_dir']);sys.path.insert(0,str(code))
    import ldscore.sumstats as ss
    import numpy as np
    import pandas as pd
    import scipy
    libraries=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
    if not sys.version.startswith('3.9.23') or libraries!={k:plan['environment'][k] for k in libraries}:raise RuntimeError('pinned environment differs')
    from sensitivity_capture_common import IntersectionCapture,clean,scalar_attrs
    capture=IntersectionCapture(ss,np,pd);capture.install()
    original_rg,original_h2,original_pair=ss.estimate_rg,ss.estimate_h2,ss._rg
    pending=[]
    def pair(frame,args,log,M,ref,w,i):
        identity=capture.pair_frame_identity(frame,args)
        result=original_pair(frame,args,log,M,ref,w,i)
        if result is not None and identity['two_step_masks'] is not None:
            for name in ['hsq1','hsq2','gencov']:
                actual=getattr(getattr(result,name),'twostep_filtered',None)
                expected=identity['final_ordered_SNP_count']-identity['two_step_masks'][name]['count']
                if actual!=expected:raise RuntimeError('stock fitted two-step mask disagrees with read-only captured mask')
        return result
    def write(args,rows,identities):
        after={p:sha(p) for p in job['inputs']}
        if before!=after or sha(known.plan)!=known.plan_sha:raise RuntimeError('input/plan changed while estimating')
        record=dict(schema='stock_sensitivity_full_precision_capture_v1',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    job_id=known.job_id,plan_sha256=known.plan_sha,arguments=vars(args),python=sys.version,
                    libraries=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
                    estimates=rows,final_intersections=identities,input_sha256_before=before,input_sha256_after=after,
                    instrumentation='read-only stock read/merge/filter observations and original estimator return capture; pinned files unchanged')
        pending.append((args.out,record))
    def rg(args,log):
        results=original_rg(args,log)
        paths=args.rg.split(',');rows=[]
        for path,obj in zip(paths[1:],results):
            row=dict(p1=paths[0],p2=path,status='ESTIMATION_FAILED' if obj is None else 'NATIVE_ESTIMATE_RETURNED')
            if obj is not None:
                row.update({k:clean(getattr(obj,k,None),np) for k in ['rg_ratio','rg','rg_jknife','rg_se','z','p','intercept_gencov']})
                row['rg']=row['rg_ratio']
                row['rg_alias_scope']='rg aliases the stock rg_ratio point estimate; rg_jknife remains separately preserved'
                for name in ['hsq1','hsq2','gencov']:row[name]=scalar_attrs(getattr(obj,name),np)
            rows.append(row)
        write(args,rows,capture.rows)
        return results
    def h2(args,log):
        result=original_h2(args,log)
        row=dict(input=args.h2,**scalar_attrs(result,np))
        mask=capture.h2_identity['two_step_hsq_mask']
        if mask is not None and result.twostep_filtered!=mask['total']-mask['count']:raise RuntimeError('h2 two-step mask mismatch')
        write(args,[row],[capture.h2_identity])
        return result
    ss._rg=pair;ss.estimate_rg=rg;ss.estimate_h2=h2
    # The stock main parses and fits using exactly the frozen forwarded arguments.
    sys.argv=[str(code/'ldsc.py')]+forwarded
    stock=runpy.run_path(str(code/'ldsc.py'),run_name='__main__')
    stock['log'].log_fh.flush();stock['log'].log_fh.close()
    if len(pending)!=1:raise RuntimeError('capture dispatch cardinality mismatch')
    output,record=pending[0]
    record['warnings_and_errors']=[dict(line=i,text=line) for i,line in enumerate(Path(output+'.log').read_text().splitlines(),1)
                                    if 'WARNING' in line or 'ERROR' in line or 'out of bounds' in line]
    record['stock_log_sha256']=sha(output+'.log')
    record['stock_delete_array_sha256']={str(p):sha(p) for p in sorted(Path(output).parent.glob(Path(output).name+'*'))
                                       if p.is_file() and str(p).endswith(('.delete','.part_delete'))}
    if {p:sha(p) for p in job['inputs']}!=before or sha(known.plan)!=known.plan_sha:raise RuntimeError('post-main identity changed')
    with Path(output+'.full_precision.json').open('x') as f:json.dump(clean(record,np),f,indent=2,allow_nan=False,default=str);f.write('\n')


if __name__=='__main__':main()
