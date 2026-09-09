"""Indexed, full-cis QTL summary statistics and GWAS/QTL locus preparation.

Not a bulk-expression or differential-expression substitute for colocalization.
Requires a source card asserting complete cis summary data, not significant-only
QTL exports. Coordinates and alleles, not rsID spelling, determine matches.
"""
from __future__ import annotations
import collections
import math
import sqlite3
from pathlib import Path
from .artifacts import transaction, verify_artifact
from .gwas import align, PAIR_FIELDS
from .io import read_json, read_tsv, require, check_hash, ensure_free, write_json, write_tsv

REQUIRED=['FEATURE','CHR','BP','A1','A2','BETA','SE','P','N','EAF']


def index(source_path, root, name, *, synthetic=False):
    c=read_json(source_path)
    require(c.get('synthetic',False)==synthetic,'Synthetic QTL source mismatch')
    for key in ('study_id','genome_build','ancestry','molecular_trait','tissue','source_uri','citation','n_semantics'):
        require(c.get(key) not in (None,'','UNRESOLVED'),'Missing QTL metadata: '+key)
    require(c['genome_build'] in {'GRCh37','GRCh38'} and c['ancestry']=='EUR','Unsupported QTL build/ancestry')
    require(c['n_semantics']=='total','QTL N must be sourced total sample size')
    require(c.get('complete_cis_summary') is True,'Significant-only QTL exports cannot be fine-mapped or colocalized')
    require(synthetic or c.get('access_permitted') is True,'QTL access not authorized')
    mapping=c['column_map']
    require(set(mapping)==set(REQUIRED) and len(set(mapping.values()))==len(mapping),'Incomplete/ambiguous QTL column map')
    check_hash(c['path'],c['sha256'])
    ensure_free(root,0 if synthetic else c.get('minimum_free_bytes',2**30))
    with transaction(root,name,stage='qtl_index',inputs=[source_path,c['path']],
                     parameters=c,synthetic=synthetic) as (work,meta):
        db=sqlite3.connect(work/'qtl.sqlite')
        db.execute('PRAGMA cache_size=-65536');db.execute('PRAGMA temp_store=FILE')
        # Duplicate variant aliases at the same coordinate are safely quarantined.
        db.execute('''CREATE TABLE qtl (feature TEXT, chr INTEGER, bp INTEGER, a1 TEXT, a2 TEXT,
            beta REAL,se REAL,p REAL,n REAL,eaf REAL, seen INTEGER DEFAULT 1,
            PRIMARY KEY(feature,chr,bp)) WITHOUT ROWID''')
        counts=collections.Counter();batch=[]
        def flush():
            db.executemany('''INSERT INTO qtl(feature,chr,bp,a1,a2,beta,se,p,n,eaf) VALUES (?,?,?,?,?,?,?,?,?,?)
                 ON CONFLICT(feature,chr,bp) DO UPDATE SET seen=seen+1''',batch)
            db.commit();batch.clear()
        try:
            for raw in read_tsv(c['path'],mapping.values()):
                counts['input_rows']+=1
                r={k:raw[v] for k,v in mapping.items()}
                try:
                    ch=int(r['CHR'].lower().removeprefix('chr'));bp=int(r['BP'])
                    b,se,p,n,eaf=map(float,(r['BETA'],r['SE'],r['P'],r['N'],r['EAF']))
                    require(r['FEATURE'] and not any(x.isspace() for x in r['FEATURE']),'Invalid feature ID')
                    require(ch in range(1,23) and bp>0,'Invalid QTL coordinate')
                    require(all(map(math.isfinite,(b,se,p,n,eaf))) and se>0 and n>0 and 0<=p<=1 and 0<eaf<1,'Invalid QTL numbers')
                    require(align(r['A1'],r['A2'],r['A1'],r['A2']) is not None,'Ambiguous QTL allele')
                except (ValueError,TypeError):
                    counts['invalid_or_ambiguous_rows']+=1;continue
                batch.append((r['FEATURE'],ch,bp,r['A1'].upper(),r['A2'].upper(),b,se,p,n,eaf))
                if len(batch)>=20000:flush()
            flush()
            counts['duplicated_coordinates']=db.execute('SELECT count(*) FROM qtl WHERE seen>1').fetchone()[0]
            db.execute('DELETE FROM qtl WHERE seen>1')
            db.execute('CREATE INDEX qtl_region ON qtl(chr,bp,feature)');db.commit()
            counts['retained_rows']=db.execute('SELECT count(*) FROM qtl').fetchone()[0]
            require(counts['retained_rows']>0,'No usable QTL rows')
            feature_rows=({'feature_id':r[0],'n_variants':r[1]} for r in db.execute('SELECT feature,count(*) FROM qtl GROUP BY feature'))
            write_tsv(work/'features.tsv',['feature_id','n_variants'],feature_rows)
        finally:
            db.close()
        meta.update(scientific_status='PASS',qtl_source={k:c[k] for k in ('study_id','genome_build','ancestry','molecular_trait','tissue','source_uri','citation','n_semantics')})
        write_json(work/'qc.json',dict(counts))
    return Path(root)/name


def locus(gwas_dir,qtl_dir,feature,region,root,name,*,minimum_variants=20,minimum_overlap=.8,max_eaf_difference=.15):
    """Emit a GWAS/QTL pair consumable by index-pair -> prepare-locus -> SuSiE.

    This local pair must never enter a genome-wide PLACO parameter estimator.
    Both original N columns remain intact; no sample-size averaging is done.
    """
    gwas_dir,qtl_dir=Path(gwas_dir).resolve(),Path(qtl_dir).resolve()
    gr,qr=verify_artifact(gwas_dir),verify_artifact(qtl_dir)
    require(gr['stage']=='normalize' and qr['stage']=='qtl_index','Unexpected QTL locus inputs')
    require(gr['synthetic']==qr['synthetic'],'Synthetic QTL contamination')
    for key in ('genome_build','ancestry'):
        require(gr['source'][key]==qr['qtl_source'][key],'QTL/GWAS metadata mismatch: '+key)
    ch,lo,hi=[int(region[k]) for k in ('CHR','START','STOP')]
    require(ch in range(1,23) and 1<=lo<=hi and minimum_variants>=2,'Invalid QTL locus bounds')
    require(0<minimum_overlap<=1 and 0<=max_eaf_difference<=1,'Invalid QTL overlap policy')
    with transaction(root,name,stage='gwas_qtl_locus',
                     inputs=[gwas_dir/'variants.sqlite',qtl_dir/'qtl.sqlite',gwas_dir/'receipt.json',qtl_dir/'receipt.json'],
                     parameters={'feature':feature,'region':region,'minimum_variants':minimum_variants,
                                 'minimum_overlap':minimum_overlap,'max_eaf_difference':max_eaf_difference},
                     synthetic=gr['synthetic']) as (work,meta):
        db=sqlite3.connect(f'{(gwas_dir/"variants.sqlite").as_uri()}?mode=ro',uri=True)
        db.execute('ATTACH DATABASE ? AS mol',(f'{(qtl_dir/"qtl.sqlite").as_uri()}?mode=ro',))
        counts=collections.Counter()
        try:
            ng=db.execute('SELECT count(*) FROM variants WHERE chr=? AND bp BETWEEN ? AND ?',(ch,lo,hi)).fetchone()[0]
            nq=db.execute('SELECT count(*) FROM mol.qtl WHERE feature=? AND chr=? AND bp BETWEEN ? AND ?',(feature,ch,lo,hi)).fetchone()[0]
            def rows():
                sql='''SELECT g.snp,g.chr,g.bp,g.a1,g.a2,g.beta,g.se,g.p,g.n,g.eaf,
                     q.a1,q.a2,q.beta,q.se,q.p,q.n,q.eaf FROM variants g JOIN mol.qtl q
                     ON g.chr=q.chr AND g.bp=q.bp WHERE q.feature=? AND g.chr=? AND g.bp BETWEEN ? AND ? ORDER BY g.bp,g.snp'''
                for r in db.execute(sql,(feature,ch,lo,hi)):
                    counts['common_coordinates']+=1
                    a=align(r[3],r[4],r[10],r[11])
                    if a is None:
                        counts['allele_mismatch']+=1;continue
                    sign,label=a;eaf=r[16] if sign==1 else 1-r[16]
                    if r[9] is not None and abs(r[9]-eaf)>max_eaf_difference:
                        counts['frequency_mismatch']+=1;continue
                    b2=sign*r[12];counts['retained']+=1
                    vals=[*r[:10],b2,r[13],r[14],r[15],eaf,r[5]/r[6],b2/r[13],label]
                    yield dict(zip(PAIR_FIELDS,['NA' if v is None else v for v in vals]))
            write_tsv(work/'pair.tsv.gz',PAIR_FIELDS,rows())
            fraction=counts['retained']/min(ng,nq) if min(ng,nq)>0 else 0
            require(counts['retained']>=minimum_variants and fraction>=minimum_overlap,'BLOCKED_BY_QTL_COVERAGE')
        finally:db.close()
        meta.update(scientific_status='PASS',analysis_scope='GWAS_QTL_LOCUS',left=gr['source'],
                    right={**qr['qtl_source'],'trait_id':feature},qc=dict(counts))
        write_json(work/'qc.json',{**dict(counts),'gwas_region_variants':ng,'qtl_region_variants':nq,
                   'overlap_of_smaller':fraction,'note':'Prepared local summary data; not colocalization evidence yet.'})
    return Path(root)/name
