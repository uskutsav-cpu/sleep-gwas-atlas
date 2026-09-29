"""Competitive locus-level pathway tests with an exact stratified randomization null.

Within each reviewed stratum the selected labels are exchangeable, with the
number of selected loci fixed. The number hitting a pathway is hypergeometric.
Convolution gives the exact distribution of total hits across strata, avoiding
Monte-Carlo noise. This is not MAGMA/S-LDSC or proof of a mechanism.
"""
from __future__ import annotations
import collections
import math
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom
from .io import require, read_tsv, write_tsv, read_json, write_json
from .stats import bh
from .artifacts import transaction


def conditional_overlap(strata):
    """strata: iterable of (N universe, M pathway-hit, K selected, X observed)."""
    distribution=np.array([1.]);offset=0;observed=0;expected=0.;variance=0.
    informative=0; nstrata=0
    for N,M,K,X in strata:
        nstrata+=1
        require(all(isinstance(z,(int,np.integer)) for z in (N,M,K,X)),'Integer counts required')
        require(N>0 and 0<=M<=N and 0<=K<=N,'Invalid stratum size')
        lo=max(0,K-(N-M));hi=min(M,K)
        require(lo<=X<=hi,'Observed overlap outside hypergeometric support')
        support=np.arange(lo,hi+1)
        prob=hypergeom.pmf(support,N,M,K)
        require(np.isfinite(prob).all() and prob.sum()>0,'Hypergeometric numeric failure')
        prob/=prob.sum()
        distribution=np.convolve(distribution,prob);distribution/=distribution.sum()
        offset+=lo;observed+=X;expected+=K*M/N
        variance+=K*(M/N)*(1-M/N)*((N-K)/(N-1)) if N>1 else 0.
        informative+=hi>lo
    require(nstrata>0 and observed>=offset,'Empty or invalid strata')
    p=float(np.clip(distribution[max(0,observed-offset):].sum(),0.,1.))
    return {'observed_loci':observed,'expected_loci':expected,'variance':variance,
            'informative_strata':informative,'p':p,
            'fold_enrichment':observed/expected if expected>0 else None}


def load_gmt(path):
    result={}
    with Path(path).open() as f:
        for line in f:
            values=line.rstrip('\n').split('\t')
            require(len(values)>=3 and values[0] and values[1],'GMT needs named, sourced gene sets')
            require(values[0] not in result,'Duplicate pathway identifier')
            genes={x for x in values[2:] if x}
            require(genes,'Empty pathway')
            result[values[0]]={'source':values[1],'genes':genes}
    require(result,'No pathways')
    return result


def run(manifest, root, name):
    c=read_json(manifest)
    require(c.get('reviewed') is True,'Gene mapping, LD-block independence and control matching require review')
    require(c.get('unit')=='independent_ld_block','Do not count correlated SNPs as independent observations')
    pathways=load_gmt(c['gmt'])
    expected=c['expected_pairs']
    require(len(expected)==len(set(expected))>0,'Freeze distinct pair family')
    min_size=int(c.get('min_genes',10));max_size=int(c.get('max_genes',1000))
    require(1<=min_size<=max_size,'Invalid pathway size bounds')
    rows=list(read_tsv(c['units'],['pair_id','locus_id','stratum','selected','eligible','genes']))
    seen=set();groups=collections.defaultdict(list)
    for r in rows:
        require(r['pair_id'] in expected and r['locus_id'] and r['stratum'],'Unknown pair/empty locus or stratum')
        key=r['pair_id'],r['locus_id'];require(key not in seen,'Duplicate LD block in pathway universe');seen.add(key)
        require(r['selected'] in {'0','1'} and r['eligible'] in {'0','1'},'Binary flags required')
        require(r['selected']=='0' or r['eligible']=='1','Selected locus was not eligible')
        if r['eligible']=='1':groups[r['pair_id']].append(r)
    output=[]
    for pair in expected:
        rs=groups[pair];background=set().union(*(set(r['genes'].split(';'))-{''} for r in rs)) if rs else set()
        for pid,pathway in sorted(pathways.items()):
            genes=pathway['genes'] & background
            row={'pair_id':pair,'pathway_id':pid,'source':pathway['source'],
                 'n_background_genes':len(background),'n_testable_pathway_genes':len(genes),
                 'n_eligible_loci':len(rs),'n_selected_loci':sum(r['selected']=='1' for r in rs),
                 'status':'TESTED','reason':'','p':None,'family_fdr':None}
            if not rs:row.update(status='NO_COVERAGE',reason='No eligible independent loci')
            elif not min_size<=len(genes)<=max_size:
                row.update(status='SIZE_FILTERED',reason='Pathway size outside frozen testable-universe range')
            elif row['n_selected_loci']==0:
                row.update(status='NO_SELECTED_LOCI',reason='No signal-selected loci; not a test of no biological involvement')
            else:
                strata=collections.defaultdict(list)
                for r in rs:strata[r['stratum']].append(r)
                counts=[]
                for group in strata.values():
                    hit=[bool((set(r['genes'].split(';'))-{''}) & genes) for r in group]
                    selected=[r['selected']=='1' for r in group]
                    counts.append((len(group),sum(hit),sum(selected),sum(a and b for a,b in zip(hit,selected))))
                row.update(conditional_overlap(counts))
                if row['informative_strata']==0:row.update(status='NO_EXCHANGEABLE_CONTROLS',p=None,reason='Null is degenerate in all strata')
            output.append(row)
    qs=bh([r['p'] for r in output])
    for r,q in zip(output,qs):r['family_fdr']=q
    fields=['pair_id','pathway_id','source','n_background_genes','n_testable_pathway_genes',
            'n_eligible_loci','n_selected_loci','status','reason','observed_loci','expected_loci',
            'variance','informative_strata','fold_enrichment','p','family_fdr']
    with transaction(root,name,stage='pathway_exact_stratified',inputs=[manifest,c['units'],c['gmt']],
                     parameters=c,synthetic=c.get('synthetic',False)) as (work,meta):
        write_tsv(work/'pathways.tsv',fields,({k:'NA' if r.get(k) is None else r[k] for k in fields} for r in output))
        meta['scientific_status']='PASS' if all(r['status']=='TESTED' for r in output) else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'family_size':len(output),
                   'tested':sum(r['status']=='TESTED' for r in output),'method':'exact stratified hypergeometric convolution',
                   'note':'Requires reviewed exchangeability/matching; FDR family retains untested planned slots.'})
    return Path(root)/name
