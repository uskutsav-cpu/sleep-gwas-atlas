import csv,hashlib,json,math
from pathlib import Path
from scipy.stats import norm
p=Path(__file__).resolve().parent;root=p.parent; atlas=p.parents[1]
inputs={'matched':root/'matched.tsv','ffs':root/'ffs_readiness/matched_source_stats.tsv','ld':atlas/'placo_diagnostic/completed_audit/ffs_reference_ld/chr5_cross_project/pairs.tsv','ld_independent':root/'published_ffs_ld/pairwise_ld.tsv','ld_fi':root/'pairwise_ld.tsv'}
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def rows(k):
 with inputs[k].open() as f:return list(csv.DictReader(f,delimiter='\t'))
ids={'rs7703746':('G','A'),'rs1363103':('T','C'),'rs2431108':('T','C')};m={r['rsid']:r for r in rows('matched') if r['rsid'] in ids};f={r['rsid']:r for r in rows('ffs') if r['rsid'] in ids};assert set(m)==set(f)==set(ids)
source={}
for snp,(a,b) in ids.items():
 x=m[snp];y=f[snp];assert {x['effect_allele'],x['other_allele']}=={a,b}=={y['FFS_source_A1'],y['FFS_source_A2']} and x['bp37']==y['bp37']
 for t,beta,se,pv,src_a in [('FI',x['fi_beta'],x['fi_se'],x['fi_p'],x['effect_allele']),('R13_insomnia',x['r13_beta'],x['r13_se'],x['r13_p'],x['effect_allele']),('FFS',y['FFS_source_beta'],y['FFS_source_se'],y['FFS_source_p'],y['FFS_source_A1'])]:
  beta=float(beta)*(1 if src_a==a else -1);se=float(se);pv=float(pv);assert se>0 and 0<pv<=1 and beta!=0
  source[(snp,t)]={'z_p':math.copysign(norm.isf(pv/2),beta),'z_wald':beta/se,'p':pv,'oriented_beta':beta,'se':se}
ld={}
for x in rows('ld'):
 if x['FFS_marker']!='rs1363103' or x['anchor'] not in ('rs7703746','rs2431108') or x['population'] not in ('EUR','FIN'):continue
 assert x['FFS_counted_allele']=='C' and x['anchor_counted_allele']==ids[x['anchor']][0] and x['status']=='OK'
 ld[(x['anchor'],x['population'])]=-float(x['r'])
assert len(ld)==4
check1=next(x for x in rows('ld_independent') if x['variant_2']=='rs1363103');check2=next(x for x in rows('ld_fi') if x['variant_2']=='rs1363103')
assert abs(ld[('rs7703746','EUR')]-float(check1['r']))<1e-12 and abs(ld[('rs2431108','EUR')]-float(check2['r']))<1e-12
out=[]
for anchor in ('rs7703746','rs2431108'):
 for pop in ('EUR','FIN'):
  r=ld[(anchor,pop)];assert abs(r)<1
  for t in ('FI','FFS','R13_insomnia'):
   for target,condition in [(anchor,'rs1363103'),('rs1363103',anchor)]:
    for method in ('z_p','z_wald'):
     zt=source[(target,t)][method];zc=source[(condition,t)][method];res=(zt-r*zc)/math.sqrt(1-r*r)
     out.append({'trait':t,'panel':pop,'target':target,'target_allele':ids[target][0],'condition':condition,'condition_allele':ids[condition][0],'signed_r':r,'z_method':method,'target_z':zt,'condition_z':zc,'residual_z':res,'nominal_normal_p':2*norm.sf(abs(res))})
assert len(out)==48
with (p/'results.tsv').open('w',newline='') as f0:w=csv.DictWriter(f0,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'input_sha256':{k:sha(v) for k,v in inputs.items()},'counts':{'markers':3,'pairs':2,'panels':2,'traits':3,'directions':2,'z_methods':2,'rows':len(out)},'eur_reference_correlations_independently_checked':True,'warning':'Simplified post-selected diagnostic; no full conditional model or inference claim.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
