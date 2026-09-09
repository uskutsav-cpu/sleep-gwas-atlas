"""Evaluate a frozen specification grid without choosing whichever result looks best."""
from __future__ import annotations
import collections
import math
from pathlib import Path
from .io import read_json, read_tsv, write_json, write_tsv, require, check_hash
from .artifacts import transaction


def summarize(manifest, root, name):
    c=read_json(manifest)
    require(c.get('reviewed') is True and c.get('grid_defined_before_results') is True,'Freeze specification grid before sensitivity results')
    grid=list(read_tsv(c['grid'],['unit_id','specification_id','primary','metric','threshold','direction']))
    raw=list(read_tsv(c['results'],['unit_id','specification_id','status','value']))
    keys=[(r['unit_id'],r['specification_id']) for r in grid]
    require(keys and len(keys)==len(set(keys)),'Empty/duplicate sensitivity grid')
    expected=set(keys);got={};groups=collections.defaultdict(list)
    for r in grid:
        require(r['primary'] in {'0','1'} and r['direction'] in {'greater_equal','less_equal','same_sign'},'Unknown decision rule')
        threshold=float(r['threshold']);require(math.isfinite(threshold),'Invalid sensitivity threshold')
        groups[r['unit_id']].append(r)
    for r in raw:
        k=r['unit_id'],r['specification_id'];require(k in expected and k not in got,'Unexpected/duplicate sensitivity result')
        require(r['status'] in {'PASS','NO_SIGNAL','FAILED','NOT_RUN','UNDERPOWERED'},'Unknown result state')
        if r['status']=='PASS':require(math.isfinite(float(r['value'])),'Nonfinite sensitivity metric')
        got[k]=r
    decisions=[];scope=[]
    for unit,rs in sorted(groups.items()):
        primary=[r for r in rs if r['primary']=='1'];require(len(primary)==1,'Exactly one prespecified primary specification per unit')
        require(len({r['metric'] for r in rs})==len({r['direction'] for r in rs})==len({r['threshold'] for r in rs})==1,'Inconsistent comparison metric or threshold')
        base=got.get((unit,primary[0]['specification_id']))
        base_value=float(base['value']) if base and base['status']=='PASS' else None
        values=[];flags=[];missing=0
        for r in rs:
            result=got.get((unit,r['specification_id']),{'status':'NOT_RUN','value':'NA'})
            value=float(result['value']) if result['status']=='PASS' else None
            decision=None
            if value is not None:
                values.append(value);threshold=float(r['threshold'])
                decision=value>=threshold if r['direction']=='greater_equal' else value<=threshold if r['direction']=='less_equal' else value*base_value>0 if base_value is not None else None
                if decision is not None:flags.append(decision)
            else:missing+=1
            decisions.append({'unit_id':unit,'specification_id':r['specification_id'],'primary':r['primary'],
              'metric':r['metric'],'status':result['status'],'value':value if value is not None else 'NA',
              'criterion_met':decision if decision is not None else 'NA'})
        state='MISSING_OR_FAILED_SPECIFICATION' if missing else 'NO_PRIMARY_METRIC' if base_value is None else 'SENSITIVE_TO_SPECIFICATION' if len(set(flags))>1 else 'CONSISTENT_WITHIN_TESTED_GRID'
        scope.append({'unit_id':unit,'status':state,'planned_specifications':len(rs),'missing_or_failed':missing,
          'metric_min':min(values) if values else 'NA','metric_max':max(values) if values else 'NA',
          'criterion_met_in_all':all(flags) if len(flags)==len(rs) else 'NA'})
    with transaction(root,name,stage='robustness_family',inputs=[manifest,c['grid'],c['results']],parameters=c,synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'specifications.tsv',list(decisions[0]),decisions)
        write_tsv(work/'robustness.tsv',list(scope[0]),scope)
        meta['scientific_status']='PASS' if all(r['status']=='CONSISTENT_WITHIN_TESTED_GRID' for r in scope) else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'units':len(scope),'planned_specifications':len(grid),
            'note':'Consistency may include consistently negative results; neither small P values nor this check establish causality. No alternative replaced the primary analysis.'})
    return Path(root)/name
