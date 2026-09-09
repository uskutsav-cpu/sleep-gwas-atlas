"""Disk-backed joint PLACO family correction, and full-scope local-rg accounting.

No tail approximation and no method substitution. Consume hash-bound native
results; assign no discovered loci to a missing or failed analysis.
"""
from __future__ import annotations
import collections
import math
import sqlite3
from pathlib import Path
from .artifacts import transaction, verify_artifact
from .io import require, read_json, read_tsv, write_tsv, write_json, check_hash
from .stats import bh


def joint_placo(manifest, root, name):
    c=read_json(manifest)
    require(c.get('reviewed') is True,'Review joint family, legacy compatibility and numerical policy')
    require(c.get('analysis_scope')=='GENOMEWIDE_ALL_SELECTED_PAIRS','Candidate-only family is not allowed')
    expected=c['expected_pairs']
    require(len(expected)==len(set(expected))>0,'Freeze unique expected pairs')
    sources={}
    inputs=[manifest]
    for r in c['results']:
        require(r['pair_id'] in expected and r['pair_id'] not in sources,'Unexpected/duplicate pair source')
        require(r.get('method') in {'PLACO','PLACO_PLUS'},'Preserve actual method identity')
        check_hash(r['path'],r['sha256']);inputs.append(r['path']);sources[r['pair_id']]=r
        if not c.get('synthetic',False):
            require(r.get('receipt') and r.get('fingerprint'),'Native source receipt must be pinned')
            rec=verify_artifact(r['receipt'],r['fingerprint'])
            require(not rec['synthetic'] and rec.get('scientific_status')=='PASS','Rejected/synthetic native result')
            require(any((Path(r['receipt'])/o['path']).resolve()==Path(r['path']).resolve() and o['sha256']==r['sha256']
                        for o in rec['outputs']),'Result is not an output of its bound artifact')
            inputs.append(Path(r['receipt'])/'receipt.json')
    require(set(sources)==set(expected),'Missing selected pair; do not shrink the planned family')
    if len({r['method'] for r in sources.values()})>1:
        require(c.get('mixed_method_compatibility_reviewed') is True and c.get('mixed_method_justification'),
                'Original PLACO and PLACO+ cannot be silently pooled across changed model assumptions')
    cap=float(c.get('maximum_failure_rate',.001));require(0<=cap<1,'Invalid failure cap')
    floor=float(c.get('genomewide_base_threshold',5e-8));require(0<floor<1,'Invalid base threshold')
    with transaction(root,name,stage='joint_placo_family',inputs=inputs,parameters=c,
                     synthetic=c.get('synthetic',False)) as (work,meta):
        db=sqlite3.connect(work/'joint.sqlite')
        db.execute('PRAGMA cache_size=-65536');db.execute('PRAGMA temp_store=FILE')
        db.execute('CREATE TABLE r(pair TEXT,snp TEXT,chr INTEGER,bp INTEGER,p REAL,status TEXT,q REAL,PRIMARY KEY(pair,snp)) WITHOUT ROWID')
        counts=[]
        try:
            for pair in expected:
                source=sources[pair];n=failed=excluded=0;batch=[]
                for row in read_tsv(source['path'],['SNP','CHR','BP','P_PLACO','status']):
                    n+=1;status=row['status']
                    require(status in {'TESTED','NUMERICAL_FAILURE','EXCLUDED_EXTREME_Z'},'Unknown row status')
                    ch=int(row['CHR']);bp=int(row['BP'])
                    require(ch in range(1,23) and bp>=1 and row['SNP'],'Invalid autosomal variant')
                    p=None
                    if status=='TESTED':
                        p=float(row['P_PLACO'])
                        # PLACO's upstream extreme-Z policy bounds the input.
                        # A numerical zero must be investigated, not reported as exact evidence.
                        require(math.isfinite(p) and 0<p<=1,'Unresolved zero/nonfinite/out-of-range PLACO P')
                    failed+=status=='NUMERICAL_FAILURE';excluded+=status=='EXCLUDED_EXTREME_Z'
                    batch.append((pair,row['SNP'],ch,bp,p,status))
                    if len(batch)>=20000:
                        db.executemany('INSERT INTO r(pair,snp,chr,bp,p,status) VALUES(?,?,?,?,?,?)',batch)
                        db.commit();batch=[]
                db.executemany('INSERT INTO r(pair,snp,chr,bp,p,status) VALUES(?,?,?,?,?,?)',batch);db.commit()
                require(n==int(source['expected_rows']),f'Missing/extra variants for {pair}')
                m=n-excluded;require(m>0,'No eligible variants')
                counts.append({'pair_id':pair,'rows':n,'eligible':m,'failed':failed,'excluded':excluded,
                               'failure_rate':failed/m,'method':source['method']})
            denominator=sum(r['eligible'] for r in counts)
            db.execute('CREATE INDEX p_index ON r(p)')
            db.execute("CREATE TEMP TABLE ranking AS SELECT pair,snp,p,ROW_NUMBER() OVER(ORDER BY p,pair,snp) AS k FROM r WHERE status='TESTED'")
            q=1.;batch=[]
            for pair,snp,p,k in db.execute('SELECT pair,snp,p,k FROM ranking ORDER BY k DESC'):
                q=min(q,p*denominator/k);batch.append((q,pair,snp))
                if len(batch)>=20000:
                    db.executemany('UPDATE r SET q=? WHERE pair=? AND snp=?',batch);batch=[]
            db.executemany('UPDATE r SET q=? WHERE pair=? AND snp=?',batch);db.commit()
            alpha=floor/len(expected)
            fields=['pair_id','SNP','CHR','BP','P_PLACO','JOINT_FAMILY_BH','status','genomewide_threshold_pass']
            def export():
                for pair,snp,ch,bp,p,status,q in db.execute('SELECT pair,snp,chr,bp,p,status,q FROM r ORDER BY pair,chr,bp,snp'):
                    yield dict(zip(fields,[pair,snp,ch,bp,p if p is not None else 'NA',q if q is not None else 'NA',status,p is not None and p<alpha]))
            write_tsv(work/'joint_results.tsv.gz',fields,export())
            write_tsv(work/'denominators.tsv',list(counts[0]),counts)
            status='PASS' if all(r['failure_rate']<=cap for r in counts) else 'FAILED_QC_NOT_CONSUMED'
            meta['scientific_status']=status
            write_json(work/'status.json',{'status':status,'n_pairs':len(expected),'denominator':denominator,
                'headline_threshold':alpha,'joint_bh_interpretation':'BH requires appropriate dependence assumptions; genome-wide cross-pair threshold is separately reported.',
                'legacy_outputs_modified':False,'independent_loci':'NOT_COMPUTED_BY_THIS_STEP'})
        finally:db.close()
    return Path(root)/name


def local_family(manifest, root, name):
    c=read_json(manifest)
    require(c.get('reviewed') is True,'Local univariate and joint-family testing rules need review')
    plan=list(read_tsv(c['planned_units'],['pair_id','locus_id']))
    keys=[(r['pair_id'],r['locus_id']) for r in plan]
    require(keys and len(keys)==len(set(keys)),'Duplicate/empty planned local family')
    rows=list(read_tsv(c['results'],['pair_id','locus_id','status','p','local_rg']))
    got={}; keyset=set(keys)
    allowed={'TESTED','UNIVARIATE_UNDERPOWERED','FAILED','NOT_RUN','NO_OVERLAP'}
    for r in rows:
        key=r['pair_id'],r['locus_id'];require(key in keyset and key not in got,'Unexpected/duplicate local unit')
        require(r['status'] in allowed,'Unknown local result status')
        if r['status']=='TESTED':
            p=float(r['p']);rg=float(r['local_rg'])
            require(math.isfinite(p) and 0<=p<=1 and math.isfinite(rg) and abs(rg)<=1,'Invalid local estimate')
        got[key]=r
    output=[]
    for key in keys:
        r=got.get(key,{'pair_id':key[0],'locus_id':key[1],'status':'NOT_RUN','p':'NA','local_rg':'NA'})
        output.append({k:r[k] for k in ['pair_id','locus_id','status','p','local_rg']})
    # Entire prespecified pair x block family. No correction on successful loci only.
    q=bh([float(r['p']) if r['status']=='TESTED' else None for r in output])
    for r,value in zip(output,q):r['family_fdr']=value if value is not None else 'NA'
    failures=sum(r['status'] in {'FAILED','NOT_RUN','NO_OVERLAP'} for r in output)
    under=sum(r['status']=='UNIVARIATE_UNDERPOWERED' for r in output)
    cap=float(c.get('maximum_failure_rate',.01));require(0<=cap<1,'Invalid local failure policy')
    status='FAILED_QC_NOT_CONSUMED' if failures/len(output)>cap else 'INSUFFICIENT_EVIDENCE' if failures or not any(r['status']=='TESTED' for r in output) else 'PASS'
    with transaction(root,name,stage='local_rg_family',inputs=[manifest,c['planned_units'],c['results']],parameters=c,
                     synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'local_family.tsv',['pair_id','locus_id','status','p','local_rg','family_fdr'],output)
        meta['scientific_status']=status
        write_json(work/'status.json',{'status':status,'planned_units':len(keys),'execution_failures_or_missing':failures,
               'univariate_underpowered':under,'failure_rate':failures/len(keys),
               'note':'Low local heritability is not an execution error or evidence of zero local sharing; report separately.'})
    return Path(root)/name
