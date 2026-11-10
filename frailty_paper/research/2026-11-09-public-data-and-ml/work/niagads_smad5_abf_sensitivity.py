#!/usr/bin/env python3
"""Postselected single-signal ABF sensitivity, not a clinical prediction."""
import csv, gzip, hashlib, json, math, sys
from pathlib import Path

root=Path(__file__).resolve().parent.parent
out=root/'outputs'
assert len(sys.argv) <= 2 and (len(sys.argv) == 1 or sys.argv[1] == '--tad')
suffix='_tad' if len(sys.argv)==2 else ''
def logsum(v):
    a=max(v); return a+math.log(sum(math.exp(x-a) for x in v))
def bf(beta,se,sd):
    v=se*se; w=sd*sd; r=w/(w+v); return .5*(math.log1p(-r)+r*(beta/se)**2)
def posterior(a,b,p12):
    a1,b1=logsum(a),logsum(b); both=logsum([x+y for x,y in zip(a,b)])
    delta=a1+b1-both
    assert delta>0
    distinct=a1+b1+math.log1p(-math.exp(-delta))
    h=[0,math.log(1e-4)+a1,math.log(1e-4)+b1,math.log(1e-8)+distinct,math.log(p12)+both]
    norm=logsum(h); return [math.exp(x-norm) for x in h]
# Check distinct-versus-shared term against explicit enumeration.
a=[-2.,1.3,.7];b=[.1,2.1,-.5]; p=1e-5
terms=[1,1e-4*sum(map(math.exp,a)),1e-4*sum(map(math.exp,b)),
       1e-8*sum(math.exp(x+y) for i,x in enumerate(a) for j,y in enumerate(b) if i!=j),
       p*sum(math.exp(x+y) for x,y in zip(a,b))]
assert max(abs(x-y/sum(terms)) for x,y in zip(posterior(a,b,p),terms))<1e-13
result=[];summary={}
qfile=out/f'stage6_niagads_smad5_exc_region{suffix}.tsv.gz'
with gzip.open(qfile,'rt') as f:
    all_qtl=list(csv.DictReader(f,delimiter='\t'))
qtl_eligible=[]
for q in all_qtl:
    c,pos,ref,alt=q['variant_id'].split(':'); af=float(q['non_ref_af'])
    if len(ref)==len(alt)==1 and {ref,alt} not in ({'A','T'},{'C','G'}) and .01<=af<=.99:
        qtl_eligible.append(q)
all_qbf=[bf(float(q['beta_non_ref']),float(q['beta_se_non_ref']),.15) for q in qtl_eligible]
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',.1232461359076205),('FinnGen_insomnia',.2)]:
    path=out/f'stage6_niagads_smad5_exc_{trait}_aligned{suffix}.csv'
    with path.open() as f: full=list(csv.DictReader(f))
    summary[trait]={}
    for filter_name,maxdelta in [('unfiltered_frequency',1.0),('abs_af_difference_le_0.1',.1)]:
        rows=[r for r in full if abs(float(r['qtl_alt_af'])-float(r['gwas_af_alt']))<=maxdelta]
        assert len(rows)>100 and len({r['variant_GRCh38'] for r in rows})==len(rows)
        matched_qbf=[bf(float(r['qtl_beta_alt']),float(r['qtl_se']),.15) for r in rows]
        support=math.exp(logsum(matched_qbf)-logsum(all_qbf))
        summary[trait][filter_name]={'n':len(rows),'qtl_abf_support_fraction_within_eligible_window':support,
                                     'qtl_lead_present':any(r['variant_GRCh38']=='chr5:136206558:C:T' for r in rows)}
        for qsd in [.075,.15,.30]:
            qbf=[bf(float(r['qtl_beta_alt']),float(r['qtl_se']),qsd) for r in rows]
            gbf=[bf(float(r['gwas_beta_alt']),float(r['gwas_se']),sd) for r in rows]
            for p12 in [1e-6,1e-5,5e-5]:
                h=posterior(qbf,gbf,p12)
                result.append({'trait':trait,'filter':filter_name,'n':len(rows),'qtl_prior_sd':qsd,
                               'gwas_prior_sd':sd,'p12':p12,'H0':h[0],'H1':h[1],'H2':h[2],'H3':h[3],'H4':h[4]})
path=out/f'stage6_niagads_smad5_exc_abf_sensitivity{suffix}.csv'
with path.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(result[0]));writer.writeheader();writer.writerows(result)
receipt={'inputs_sha256':{qfile.name:hashlib.sha256(qfile.read_bytes()).hexdigest(),
         **{f'stage6_niagads_smad5_exc_{t}_aligned{suffix}.csv':hashlib.sha256((out/f'stage6_niagads_smad5_exc_{t}_aligned{suffix}.csv').read_bytes()).hexdigest() for t in ['Mak_FI','Mak_FP','FinnGen_insomnia']}},
         'qtl_eligible_in_window':len(qtl_eligible),'coverage':summary,'result_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'limitations':'Exploratory postselected single-signal model; disease-window support unknown; no causal or predictive validation.'}
(out/f'stage6_niagads_smad5_exc_abf_receipt{suffix}.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'coverage':summary,'base':[r for r in result if r['filter']=='abs_af_difference_le_0.1' and r['qtl_prior_sd']==.15 and r['p12']==1e-5]},indent=2))
