#!/usr/bin/env python3
"""Independent numerical implementation, with no import from audit generator.

Normal tails computed as regularized upper incomplete gamma Q(1/2,z²/2),
using a convergent series or continued fraction, instead of erfc. BH uses
forward ranks and explicit suffix minima. Checks historical source tables
against generated outputs. No native scientific replay is implied.
"""
import argparse
import csv
import json
import math
from pathlib import Path

def tail(z):
    if not math.isfinite(z):
        raise ValueError('nonfinite Z')
    x = z*z/2
    if x==0:
        return 1.0
    a=.5
    pref=math.exp(a*math.log(x)-x-math.lgamma(a))
    if x<a+1:
        delta=1/a; total=delta; ap=a
        for _ in range(10000):
            ap+=1;delta*=x/ap;total+=delta
            if abs(delta)<abs(total)*2e-16:
                return max(0,1-total*pref)
    else:
        b=x+1-a;c=1e300;d=1/b;h=d
        for i in range(1,10000):
            an=-i*(i-a);b+=2;d=an*d+b
            if abs(d)<1e-300:d=1e-300
            c=b+an/c
            if abs(c)<1e-300:c=1e-300
            d=1/d;delta=d*c;h*=delta
            if abs(delta-1)<3e-16:
                return pref*h
    raise ArithmeticError('incomplete gamma did not converge')

def adjustment(ps):
    if not ps or any(not math.isfinite(x) or not 0<=x<=1 for x in ps):
        raise ValueError('invalid P family')
    indexed=sorted(enumerate(ps),key=lambda x:x[1]);m=len(ps)
    candidates=[min(1,p*m/(rank+1)) for rank,(_,p) in enumerate(indexed)]
    suffix=[min(candidates[rank:]) for rank in range(m)]
    answer=[None]*m
    for (index,_),value in zip(indexed,suffix):answer[index]=value
    return answer

def rows(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def close(x,y):
    if not math.isclose(float(x),float(y),rel_tol=3e-12,abs_tol=1e-300):
        raise AssertionError(f'numerical mismatch: {x} versus {y}')

def audit_figures(root):
    """Cross-check all six figure data exports against the separately audited tables."""
    out=root/'discovery_extension/sleep_submission_evidence_v1'
    src=out/'figures/source_data'
    expected=['fig1_atlas','fig2_outcome_replication','fig3_flow','fig4_heterogeneity','fig5_domains','graphical_abstract']
    if any(not (src/(name+'.tsv')).is_file() for name in expected):
        return {'status':'BLOCKED_MISSING_FIGURE_SOURCE_DATA'}
    g=rows(out/'tables/global_1200.tsv');r=rows(out/'tables/replication_family_217.tsv')
    gd={(x['sleep_trait'],x['extension_trait_id']):x for x in g}
    rd={x['pair_id']:x for x in r}
    atlas=rows(src/'fig1_atlas.tsv')
    keys=[(x['sleep_trait'],x['extension_trait_id']) for x in atlas]
    if len(keys)!=len(set(keys)) or set(keys)!=set(gd):raise AssertionError('atlas pair universe mismatch')
    for x in atlas:
        y=gd[(x['sleep_trait'],x['extension_trait_id'])]
        for key in ('rg','se','p','extension_fdr'):close(x[key],y[key])
    for name,eligible in [('fig2_outcome_replication',{x['pair_id'] for x in r if x['replication_class']=='REPLICATED'}),
                          ('fig4_heterogeneity',{x['pair_id'] for x in r if x['replication_class']=='REPLICATED' and float(x['heterogeneity_p'])<.05})]:
        f=rows(src/(name+'.tsv'));ids=[x['pair_id'] for x in f]
        if len(ids)!=len(set(ids)) or set(ids)!=eligible:raise AssertionError('forest selection mismatch')
        for x in f:
            y=rd[x['pair_id']]
            for key in ('discovery_rg','discovery_se','replication_rg','replication_se','heterogeneity_p'):close(x[key],y[key])
            for key,verified in [('discovery_ci_low','discovery_ci95_low'),('discovery_ci_high','discovery_ci95_high'),('replication_ci_low','replication_ci95_low'),('replication_ci_high','replication_ci95_high')]:close(x[key],y[verified])
            if x['phenotype_match_status']!=y['phenotype_match_status']:raise AssertionError('forest definition qualification mismatch')
    domain=rows(src/'fig5_domains.tsv')
    if {x['phenotype_domain'] for x in domain}!={x['phenotype_domain'] for x in g}:raise AssertionError('domain set mismatch')
    for x in domain:
        selected=[y for y in g if y['phenotype_domain']==x['phenotype_domain']]
        if int(x['tests'])!=len(selected) or int(x['significant_pairs'])!=sum(y['historical_significant']=='True' for y in selected) or int(x['phenotypes'])!=len({y['extension_trait_id'] for y in selected}):raise AssertionError('domain count mismatch')
    candidate=len(rows(root/'discovery_extension/results/candidate_pool.tsv'))
    flow_expected={'Candidate phenotypes':candidate,'Selected phenotypes':len({x['extension_trait_id'] for x in g}),'Discovery tests':len(g),'FDR-significant pairs':sum(x['historical_significant']=='True' for x in g),'Replication candidates':len(r),'Outcome-side positives':sum(x['replication_class']=='REPLICATED' for x in r)}
    for name in ('fig3_flow','graphical_abstract'):
        f=rows(src/(name+'.tsv'))
        if {x['stage']:int(x['count']) for x in f}!=flow_expected:raise AssertionError('flow counts mismatch')
    result={'status':'PASS','source_exports_checked':len(expected),'atlas_rows_checked':len(atlas),'forest_rows_checked':len(rows(src/'fig2_outcome_replication.tsv')),'heterogeneity_rows_checked':len(rows(src/'fig4_heterogeneity.tsv')),'domain_rows_checked':len(domain),'scope':'NUMERICAL_SOURCE_MAPPING_ONLY_VISUAL_QA_SEPARATE'}
    (out/'logs/numerical_figure_source_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

def run(root):
    out=root/'discovery_extension/sleep_submission_evidence_v1'
    original=rows(root/'discovery_extension/results/ldsc/extension_rg_matrix.tsv')
    audited=rows(out/'tables/global_1200.tsv')
    if len(original)!=len(audited):raise AssertionError('family size changed')
    ps=[float(r['p']) for r in original]; q=adjustment(ps)
    independent_q=adjustment([tail(float(r['z'])) for r in original])
    for src,dst,fdr,qz in zip(original,audited,q,independent_q):
        if src['sleep_trait']!=dst['sleep_trait'] or src['extension_trait_id']!=dst['extension_trait_id']:raise AssertionError('identity changed')
        close(fdr,dst['fdr_recomputed_original_p']);close(qz,dst['fdr_from_printed_z'])
        close(tail(float(src['rg'])/float(src['se'])),dst['p_from_printed_rg_se'])
    rep=rows(out/'tables/replication_family_217.tsv');n=sum(r['replication_rg']!='NA' for r in rep)
    for r in rep:
        close(.05/len(rep),r['replication_alpha_exact'])
        if r['replication_rg']=='NA':continue
        a,sa,b,sb=[float(r[k]) for k in ('discovery_rg','discovery_se','replication_rg','replication_se')]
        delta=b-a; se=math.hypot(sa,sb)
        close(tail(float(r['replication_z'])),r['p_from_printed_z'])
        close(tail(b/sb),r['p_from_printed_rg_se'])
        close(tail(delta/se),r['heterogeneity_p_recomputed_cov0'])
        close(b-1.959963984540054*sb,r['replication_ci95_low'])
        close(b+1.959963984540054*sb,r['replication_ci95_high'])
        if r['replication_class']=='REPLICATED' and (a*b<=0 or tail(float(r['replication_z']))>=.05/len(rep)):raise AssertionError('invalid replicated label')
    receipt={'status':'PASS','implementation':'incomplete_gamma_series_continued_fraction_and_explicit_suffix_BH','global_rows_checked':len(original),'replication_family_rows_checked':len(rep),'replication_estimates_checked':n,'native_ldsc_reruns':0,'tolerance_relative':3e-12,'tolerance_absolute':1e-300,'figure_source_audit':audit_figures(root)}
    (out/'logs/numerical_independent.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    return receipt

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[3])
    print(json.dumps(run(parser.parse_args().root.resolve()),sort_keys=True))
