"""Complete accounting of six-disorder deliverables, including missing work.

A coverage audit is not a scientific seal. PASS means a verified stage receipt,
not that the result is novel, causal, or ready for publication.
"""
from __future__ import annotations
from pathlib import Path
from .artifacts import verify_artifact, transaction
from .io import read_json, require, write_json, write_tsv
from .pairs import BRAIN, verify_lock

STAGES=('dense_inputs','replication','local_rg','pleiotropy','clumping',
        'finemapping','colocalization','molecular_qtl','cell_enrichment',
        'mr_forward','mr_reverse','cross_disorder')
# Required explicit labels; no report renders an absent step as zero discoveries.
UNBOUND={'NOT_RUN','BLOCKED_BY_DATA','BLOCKED_BY_SOFTWARE','BLOCKED_BY_COMPUTE',
         'FAILED_QC','NOT_APPLICABLE_NO_SIGNAL','REVIEW_REQUIRED'}
NATIVE_METHODS={'local_rg':{'lava','hdl_l','supergnova'},'pleiotropy':{'placo'},
                'clumping':{'clump'},'finemapping':{'susie'},'colocalization':{'coloc'},
                'mr_forward':{'mr'},'mr_reverse':{'mr'}}


def template(pair_lock,out):
    lock=verify_lock(pair_lock)
    rows=[]
    for pair in lock['pairs']:
        for stage in STAGES:
            rows.append({'pair_id':pair['pair_id'],'stage':stage,'path':None,
                         'fingerprint':None,'status':'NOT_RUN','reason':'Unbound stage'})
    value={'schema_version':1,'pair_lock':str(Path(pair_lock).resolve()),'reviewed':False,
           'requirements':rows,'global_stages':[
             {'stage':'genomic_sem','path':None,'fingerprint':None,'status':'NOT_RUN','reason':'Unbound multivariate S/V and model'},
             {'stage':'conjfdr','path':None,'fingerprint':None,'status':'NOT_RUN','reason':'MATLAB and reviewed official workflow required'}],
           'note':'Explicitly bind one accepted stage-summary artifact or native unit per slot. A single passing locus does not establish completion of all loci.'}
    write_json(out,value)
    return value


def audit(manifest,root,name):
    c=read_json(manifest)
    require(c.get('schema_version')==1,'Unknown campaign schema')
    lock=verify_lock(c['pair_lock'])
    expected={(r['pair_id'],s) for r in lock['pairs'] for s in STAGES}
    supplied={}
    for r in c['requirements']:
        key=r['pair_id'],r['stage']
        require(key in expected and key not in supplied,'Unexpected/duplicate campaign slot')
        supplied[key]=r
    require(set(supplied)==expected,'Incomplete campaign requirement manifest')
    rows=[];inputs=[manifest,c['pair_lock']]
    def one(pair,stage,b):
        row={'pair_id':pair,'stage':stage,'status':b.get('status','NOT_RUN'),
             'reason':b.get('reason',''),'artifact':b.get('path') or ''}
        if b.get('path'):
            require(c.get('reviewed') is True,'Bound campaign receipts require scope review')
            require(b.get('fingerprint'),'Pin every accepted campaign artifact fingerprint')
            r=verify_artifact(b['path'],b['fingerprint'])
            require(r['synthetic']==lock['synthetic'],'Synthetic campaign contamination')
            require(b.get('scope_complete') is True,'A passing single locus is not a completed stage: verify scope')
            if stage in NATIVE_METHODS and r['stage']=='native_job':
                require(r['parameters'].get('method') in NATIVE_METHODS[stage],'Native method/stage mismatch')
            inputs.append(Path(b['path'])/'receipt.json')
            inputs.extend(Path(b['path'])/o['path'] for o in r['outputs'])
            state=r.get('scientific_status')
            require(state in {'PASS','NO_SIGNAL','INSUFFICIENT_EVIDENCE','FAILED_QC_NOT_CONSUMED'},
                    'Stage must state its scientific status, not just process completion')
            row.update(status=state,reason='Verified output bytes and declared scope; inspect scientific interpretation')
        else:
            require(row['status'] in UNBOUND,'Absent artifact cannot be marked PASS')
            require(row['reason'],'Unbound stage requires reason')
            if row['status']=='NOT_APPLICABLE_NO_SIGNAL':
                require(b.get('upstream_no_signal_receipt'),'No-signal applicability needs upstream evidence')
                rr=verify_artifact(b['upstream_no_signal_receipt'])
                require(rr['synthetic']==lock['synthetic'] and rr.get('scientific_status')=='NO_SIGNAL',
                        'No-signal receipt is not valid')
                inputs.append(Path(b['upstream_no_signal_receipt'])/'receipt.json')
        return row
    for pair,stage in sorted(expected):
        rows.append(one(pair,stage,supplied[(pair,stage)]))
    global_rows=c.get('global_stages',[])
    require({r['stage'] for r in global_rows}=={'genomic_sem','conjfdr'} and len(global_rows)==2,
            'Keep Genomic SEM and conjunction-FDR visible')
    for r in global_rows:rows.append(one('GLOBAL',r['stage'],r))
    selected={r['disease_trait'] for r in lock['pairs']}
    for disease in BRAIN:
        if disease not in selected:
            rows.append({'pair_id':disease,'stage':'pair_selection','status':'NO_ELIGIBLE_PAIR',
                         'reason':'No global-screen connection selected; do not manufacture a positive result','artifact':''})
    complete_states={'PASS','NO_SIGNAL','NOT_APPLICABLE_NO_SIGNAL','NO_ELIGIBLE_PAIR'}
    resolved=sum(r['status'] in complete_states for r in rows)
    with transaction(root,name,stage='campaign_coverage',inputs=inputs,
                     parameters={'expected_requirements':len(rows)},synthetic=lock['synthetic']) as (work,meta):
        write_tsv(work/'coverage.tsv',['pair_id','stage','status','reason','artifact'],rows)
        # Never turn a coverage checklist into a scientific release certification.
        meta['scientific_status']='INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':'REVIEW_REQUIRED' if resolved==len(rows) else 'INCOMPLETE',
                  'requirements':len(rows),'resolved':resolved,'unresolved':len(rows)-resolved,
                  'synthetic':lock['synthetic'],'publication_ready':False})
        from html import escape
        with (work/'coverage.html').open('w') as f:
            f.write('<!doctype html><meta charset="utf-8"><title>Brain6 coverage</title><h1>Brain6 execution coverage</h1>')
            f.write('<p>'+('SYNTHETIC SOFTWARE TEST — NOT RESEARCH RESULTS' if lock['synthetic'] else 'Receipts audit, not scientific certification')+'</p>')
            f.write('<table border="1"><tr><th>Pair</th><th>Stage</th><th>Status</th><th>Reason</th></tr>')
            for r in rows:f.write('<tr>'+''.join('<td>'+escape(str(r[k]))+'</td>' for k in ('pair_id','stage','status','reason'))+'</tr>')
            f.write('</table>')
    return Path(root)/name
