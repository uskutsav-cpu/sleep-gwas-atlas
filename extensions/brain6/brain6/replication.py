"""Evaluate sourced replication records without equating a new release with independence."""
from __future__ import annotations
import math
from pathlib import Path
from .artifacts import transaction
from .io import read_json, read_tsv, require, write_tsv, write_json, file_record, check_hash
from .pairs import verify_lock
from .stats import bh


def evaluate(manifest_path, root, name):
    c=read_json(manifest_path)
    require(c.get('reviewed') is True, 'Replication phenotype/cohort audit requires review')
    lock=verify_lock(c['pair_lock'])
    require(0<c.get('fdr_max',.05)<1 and math.isfinite(c.get('minimum_h2_z',4)) and c.get('minimum_h2_z',4)>0,'Invalid replication policy')
    expected={p['pair_id']:p for p in lock['pairs']}
    bindings={}
    for b in c.get('results',[]):
        require(b['pair_id'] in expected and b['pair_id'] not in bindings,'Unexpected/duplicate replication pair')
        bindings[b['pair_id']]=b
    rows=[];inputs=[manifest_path,c['pair_lock']]
    for pair, discovery in sorted(expected.items()):
        row={'pair_id':pair,'status':'NO_VALID_REPLICATION','discovery_rg':discovery['rg'],
             'replication_rg':None,'se':None,'p':None,'independence_scope':'NONE',
             'phenotype_match':'UNRESOLVED','source':'','reason':'No verified replication source'}
        b=bindings.get(pair)
        if b is not None:
            row.update(source=b.get('source',''),reason=b.get('reason',''),
                       independence_scope=b.get('independence_scope','NONE'),
                       phenotype_match=b.get('phenotype_match','UNRESOLVED'))
            if b.get('status')=='NO_VALID_REPLICATION':
                require(row['reason'],'Explain unavailable replication')
            else:
                require(b.get('cohort_audit') and b.get('source'), 'Missing source or cohort-overlap audit')
                require(b.get('cohort_audit_sha256'), 'Pin the cohort-overlap audit document')
                rec=file_record(b['cohort_audit'])
                require(rec['sha256']==b['cohort_audit_sha256'],'Cohort audit changed')
                inputs.append(b['cohort_audit'])
                require(row['independence_scope'] in {'BOTH_TRAITS','DISEASE_ONLY','SLEEP_ONLY','NONE'},'Unknown replication independence scope')
                eligible=(b.get('independence_reviewed') is True and row['independence_scope']!='NONE'
                          and row['phenotype_match'] in {'exact','comparable'}
                          and b.get('qc_pass') is True)
                # Distinct study labels alone do not establish distinct individuals.
                overlap=set(b.get('discovery_cohorts',[])) & set(b.get('replication_cohorts',[]))
                if overlap:
                    eligible=False;row['reason']='Overlapping cohorts: '+','.join(sorted(overlap))
                if eligible:
                    if not lock['synthetic']:
                        require(b.get('result_file') and b.get('result_sha256'),'Empirical replication needs a hash-bound result table')
                        check_hash(b['result_file'],b['result_sha256'])
                        measured=list(read_tsv(b['result_file'],['rg','se','p']))
                        require(len(measured)==1,'Bind exactly one measured replication result')
                        require(all(math.isclose(float(measured[0][key]),float(b[key]),rel_tol=1e-10,abs_tol=1e-300)
                                    for key in ('rg','se','p')),'Manifest replication estimates differ from measured file')
                        inputs.append(b['result_file'])
                    rg,se,p,h2s,h2d=[float(b[k]) for k in ('rg','se','p','sleep_h2_z','disease_h2_z')]
                    require(all(map(math.isfinite,(rg,se,p,h2s,h2d))) and abs(rg)<=1 and se>0 and 0<=p<=1,
                            'Invalid replication estimates')
                    row.update(replication_rg=rg,se=se,p=p,
                               status='UNDERPOWERED' if min(h2s,h2d)<c.get('minimum_h2_z',4) else 'TESTED',
                               reason='')
        rows.append(row)
    # Only eligible well-powered tests can be promoted; all expected pairs
    # remain in the family denominator, including unavailable and underpowered.
    qs=bh([r['p'] if r['status']=='TESTED' else None for r in rows]) if rows else []
    for row,q in zip(rows,qs):
        row['replication_family_fdr']=q
        if row['status']=='TESTED':
            concordant=row['discovery_rg']*row['replication_rg']>0
            significant=q<c.get('fdr_max',.05)
            row['status']=('REPLICATED' if significant else 'DIRECTIONALLY_SUPPORTED') if concordant else (
                'DISCORDANT' if significant else 'INCONCLUSIVE_OPPOSITE_SIGN')
    with transaction(root,name,stage='replication_family',inputs=inputs,
                     parameters={'family_size':len(rows),'policy':c},synthetic=lock['synthetic']) as (work,meta):
        fields=['pair_id','status','discovery_rg','replication_rg','se','p','replication_family_fdr',
                'independence_scope','phenotype_match','source','reason']
        write_tsv(work/'replication.tsv',fields,({k:'NA' if v is None else v for k,v in r.items()} for r in rows))
        meta['scientific_status']='PASS' if rows and all(r['status'] in {'REPLICATED','DIRECTIONALLY_SUPPORTED','DISCORDANT','INCONCLUSIVE_OPPOSITE_SIGN'} for r in rows) else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'expected_pairs':len(rows),
                 'replicated':sum(r['status']=='REPLICATED' for r in rows),
                 'note':'Disease-only replication reuses the sleep GWAS; it is not an independent replication of both cohorts.'})
    return Path(root)/name
