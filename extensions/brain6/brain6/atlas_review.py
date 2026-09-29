"""Review an actual archived/host atlas without substituting proposal examples.

This consumes summary results, not GWAS genotypes. It is a reanalysis of the
provided atlas; it does not rerun LDSC or imply independence of the studies.
"""
from __future__ import annotations
import collections
import math
from pathlib import Path
import numpy as np
from .artifacts import transaction
from .io import require, read_tsv, write_tsv, write_json, file_record, read_json
from .pairs import BRAIN, SLEEP, rank_pairs
from .stats import bh


def review(repo, root, name, *, archive_date=None, synthetic=False):
    repo=Path(repo).resolve()
    matrix=repo/'results/tables/rg_matrix.tsv'
    panel=repo/'config/analysis_panel.tsv'
    rows=list(read_tsv(matrix,['sleep_trait','disease_trait','rg','se','p','fdr',
                               'interpretation_status','analysis_tier']))
    traits=list(read_tsv(panel,['trait_id','domain','raw_file','source_id']))
    trait_ids=[r['trait_id'] for r in traits]
    require(len(trait_ids)==len(set(trait_ids))==45,'Expected locked 45-trait panel')
    sleeps={r['trait_id'] for r in traits if r['domain']=='sleep'}
    diseases=set(trait_ids)-sleeps
    require(sleeps==set(SLEEP) and len(diseases)==33,'Panel is not the known 12 x 33 family')
    keys=[(r['sleep_trait'],r['disease_trait']) for r in rows]
    require(len(keys)==len(set(keys))==396,'Missing/duplicate full atlas rows')
    require(set(keys)=={(s,d) for s in sleeps for d in diseases},'Wrong atlas Cartesian product')
    pv=[float(r['p']) for r in rows]; actual=bh(pv)
    for r,q in zip(rows,actual):
        f=float(r['fdr'])
        require(math.isfinite(f) and math.isclose(q,f,rel_tol=1e-8,abs_tol=1e-300),
                'Supplied full-family FDR does not reproduce from the 396 raw P values')
        require(r['interpretation_status'] in {'PRIMARY','EXCLUDED_FROM_PRIMARY_INFERENCE'},
                'Unknown QC label; add a reviewed mapping, not a guess')
        require((r['interpretation_status']=='PRIMARY')==(r['analysis_tier']=='PRIMARY_PHASE1'),
                'Inconsistent tier and primary interpretation labels')
    cfg={'expected_matrix_rows':396,'matrix_qc_column':'interpretation_status',
         'primary_qc_values':['PRIMARY'],'global_fdr_max':.05,'brain_traits':BRAIN,
         'sleep_traits':SLEEP,'synthetic':synthetic,'experiment':'brain6-measured-followup',
         'native_policy':{'pleiotropy':'PLACO_PLUS','local_rg':'REVIEW_REQUIRED',
                          'preserve_track_b':True}}
    ranking=rank_pairs(matrix,cfg)
    summaries=[]; drafts=[]
    for b in BRAIN:
        group=[r for r in ranking if r['disease_trait']==b]
        eligible=[r for r in group if r['eligible']]
        best=eligible[0] if eligible else group[0]
        # Preserve the previously chosen Pair B rather than retroactively
        # claiming the smallest-P partner was the original selection.
        anchor=next((r for r in eligible if r['pair_id']=='insomnia__adhd'),None)
        chosen=anchor or (eligible[0] if eligible else None)
        summaries.append({'disease_trait':b,'n_tested_sleep_traits':len(group),
            'n_fdr_significant':len(eligible),'statistical_top_pair':best['pair_id'],
            'top_rg':best['rg'],'top_se':best['se'],'top_p':best['p'],'top_family_fdr':best['fdr'],
            'proposed_pair':chosen['pair_id'] if chosen else '',
            'decision_status':'PROPOSED_NOT_FROZEN' if chosen else 'NO_ELIGIBLE_GLOBAL_PAIR',
            'archive_date':archive_date or 'HOST_SNAPSHOT'})
        reason=('Preserve previously selected Track-B Pair B; not the smallest-FDR ADHD partner'
                if anchor else 'Smallest original full-family FDR among eligible measured partners'
                if chosen else 'No global-FDR-significant sleep partner; retain disorder as global-null, not evidence of no local sharing')
        drafts.append({'disease_trait':b,'sleep_trait':chosen['sleep_trait'] if chosen else '',
                       'role':'PRIMARY' if chosen else 'NO_ELIGIBLE_PAIR','reason':reason})
    primary=[r for r in rows if r['interpretation_status']=='PRIMARY']
    sig=[r for r in primary if float(r['fdr'])<.05]
    report={'source':'archived_results' if archive_date else 'host_results',
            'archive_date':archive_date,'full_matrix_rows':len(rows),'primary_rows':len(primary),
            'primary_significant':len(sig),'primary_positive':sum(float(r['rg'])>0 for r in sig),
            'primary_negative':sum(float(r['rg'])<0 for r in sig),
            'brain_subset_rows':len(ranking),'brain_significant':sum(r['eligible'] for r in ranking),
            'new_gwas_or_ldsc_run':False,'pair_lock_created':False,'synthetic':synthetic,
            'fdr_source_column':'fdr','fdr_denominator':396,
            'warning':'Archival evidence is not automatically the current merged data state. Proposals require review and a freeze before native jobs.'}
    with transaction(root,name,stage='atlas_review',inputs=[matrix,panel],
                     parameters={'config':cfg,'archive_date':archive_date},synthetic=synthetic) as (work,meta):
        write_tsv(work/'brain72.tsv',list(ranking[0]),ranking)
        write_tsv(work/'brain_summary.tsv',list(summaries[0]),summaries)
        write_tsv(work/'pair_decisions.DRAFT.tsv',list(drafts[0]),drafts)
        write_json(work/'atlas_profile.json',cfg)
        write_json(work/'summary.json',report)
        meta['scientific_status']='PASS'
        meta['interpretation']='Verified supplied global-screen arithmetic; no new locus inference'
    return Path(root)/name


def inspect_covariance(repo, root, name, *, synthetic=False):
    """Recompute matrix diagnostics from actual numeric exports; preserve fit failures."""
    repo=Path(repo)
    s=repo/'results/tables/ldsc_covariance_45x45.tsv'
    v=repo/'results/tables/ldsc_sampling_covariance_1035x1035.tsv.gz'
    fit=repo/'results/tables/genomic_sem_model_fit.tsv'
    inputs=[s,v,fit]
    def load(path,expected):
        import csv
        from .io import open_text
        with open_text(path) as f:
            rr=csv.reader(f,delimiter='\t');header=next(rr); data=list(rr)
        labels=[r[0] for r in data]
        require(len(data)==expected and len(header)==expected+1,'Matrix dimension/label mismatch')
        require(len(set(labels))==expected,'Duplicate matrix row labels')
        require(labels==header[1:],'Row/column matrix order differs')
        a=np.asarray([r[1:] for r in data],dtype=float)
        require(a.shape==(expected,expected) and np.isfinite(a).all(),'Invalid numeric covariance')
        symmetry=float(np.max(np.abs(a-a.T)))
        require(symmetry<1e-8,'Nonsymmetric covariance; no silent symmetrization')
        eig=np.linalg.eigvalsh(a)
        return {'rows':expected,'columns':expected,'max_symmetry_error':symmetry,
                'min_eigenvalue':float(eig[0]),'max_eigenvalue':float(eig[-1]),
                'negative_eigenvalues_below_minus_1e_10':int(np.sum(eig < -1e-10)),
                'raw_psd':bool(eig[0]>=-1e-10),'matrix_altered':False}
    ds=load(s,45);dv=load(v,1035)
    fits=list(read_tsv(fit,['model_id','validation_status','validation_reason']))
    require(fits and len({r['model_id'] for r in fits})==len(fits),'Empty/duplicate fit rows')
    counts=dict(collections.Counter(r['validation_status'] for r in fits))
    result={'S':ds,'V':dv,'fit_status_counts':counts,'fits_run_in_this_session':False,
            'matrices_recomputed_from_gwas':False,'existing_models':len(fits),
            'model_status':'REVIEW_REQUIRED','note':'Negative eigenvalues in estimated S are reported, not automatically repaired. Recheck model fit and trait-scale specification before inference.'}
    with transaction(root,name,stage='covariance_review',inputs=inputs,parameters={},synthetic=synthetic) as (work,meta):
        write_json(work/'diagnostics.json',result)
        write_tsv(work/'model_validation.tsv',list(fits[0]),fits)
        meta['scientific_status']='INSUFFICIENT_EVIDENCE'
    return Path(root)/name


def review_replication(repo,root,name):
    """Audit an existing Track-B replication ledger; do not rerun or upgrade it."""
    from scipy.stats import norm
    source=Path(repo)/'results/track_b/02_independent_global_replication.tsv'
    rows=list(read_tsv(source,['pair_id','discovery_rg','discovery_SE','replication_rg',
        'replication_SE','replication_P','replication_status','sample_overlap','claim_limit']))
    require(rows and len({r['pair_id'] for r in rows})==len(rows),'Duplicate/empty replication ledger')
    checks=[]
    for r in rows:
        row={'pair_id':r['pair_id'],'original_status':r['replication_status'],
             'numeric_recheck':'NOT_APPLICABLE','rounded_z_p':None,'reported_p':r['replication_P'],
             'sample_overlap_caveat':r['sample_overlap'],'claim_limit':r['claim_limit']}
        if r['replication_rg'] not in {'NA',''}:
            rg,se,p=map(float,(r['replication_rg'],r['replication_SE'],r['replication_P']))
            require(all(math.isfinite(x) for x in (rg,se,p)) and abs(rg)<=1 and se>0 and 0<=p<=1,'Invalid archived replication estimate')
            ci=[rg-norm.ppf(.975)*se,rg+norm.ppf(.975)*se]
            require(all(math.isclose(a,float(r[k]),abs_tol=1e-10) for a,k in zip(ci,['replication_CI_lower','replication_CI_upper'])),
                    'Archived confidence interval does not match the displayed estimate/SE')
            p_rounded=2*norm.sf(abs(rg/se))
            # Published/exported P may use unrounded effects. Report, do not
            # replace, the native P with a value derived from rounded fields.
            row.update(numeric_recheck='CI_AND_DIRECTION_CHECKED',rounded_z_p=float(p_rounded),
                       reported_p=p,direction_checked='CONCORDANT' if rg*float(r['discovery_rg'])>0 else 'NONCONCORDANT')
        checks.append(row)
    fields=sorted({k for r in checks for k in r})
    with transaction(root,name,stage='archived_replication_review',inputs=[source],parameters={},synthetic=False) as (work,meta):
        write_tsv(work/'replication_ledger.tsv',list(rows[0]),rows)
        write_tsv(work/'arithmetic_checks.tsv',fields,({k:'NA' if r.get(k) is None else r[k] for k in fields} for r in checks))
        write_json(work/'summary.json',{'rows':len(rows),'new_replication_run':False,'statuses_preserved':True,
            'status_counts':dict(collections.Counter(r['replication_status'] for r in rows)),
            'note':'Archival evidence only. Directional support is not promoted to verified two-cohort independence or causal mechanism.'})
        meta['scientific_status']='PASS'
    return Path(root)/name
