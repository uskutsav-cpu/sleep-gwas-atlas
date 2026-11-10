import csv,json,math,pathlib,hashlib

def logsum(x):
 m=max(x);return m+math.log(sum(math.exp(y-m) for y in x))
def bf(beta,se,prior_sd):
 v=se*se;w=prior_sd**2;r=w/(w+v);z=beta/se
 return .5*(math.log1p(-r)+r*z*z)
def posterior(l1,l2,p1,p2,p12):
 s1=logsum(l1);s2=logsum(l2);s4=logsum([x+y for x,y in zip(l1,l2)]);d=s1+s2-s4
 assert d>0
 s3=s1+s2+math.log1p(-math.exp(-d))
 h=[0,math.log(p1)+s1,math.log(p2)+s2,math.log(p1*p2)+s3,math.log(p12)+s4];z=logsum(h);return [math.exp(x-z) for x in h]
# Test H3 algebra against explicit enumeration of distinct-SNP configurations.
a=[-2.0,1.3,.7];b=[.1,2.1,-.5];h=posterior(a,b,1e-4,1e-4,1e-5);terms=[1,1e-4*sum(map(math.exp,a)),1e-4*sum(map(math.exp,b)),1e-8*sum(math.exp(x+y) for i,x in enumerate(a) for j,y in enumerate(b) if i!=j),1e-5*sum(math.exp(x+y) for x,y in zip(a,b))];assert max(abs(x-y/sum(terms)) for x,y in zip(h,terms))<1e-13
p=pathlib.Path('work/stage3_smim32_FinnGen_insomnia_aligned.csv');r=list(csv.DictReader(p.open()));assert len(r)>100 and len({x['variant'] for x in r})==len(r)
res=[]
for qtl_sd in (.075,.15,.30):
 l1=[bf(float(x['qtl_beta_alt']),float(x['qtl_se']),qtl_sd) for x in r];l2=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),.20) for x in r]
 for p12 in (1e-6,1e-5,5e-5):
  h=posterior(l1,l2,1e-4,1e-4,p12);res.append({'qtl_prior_effect_sd':qtl_sd,'insomnia_prior_log_odds_sd':.2,'p1':.0001,'p2':.0001,'p12':p12,'H0':h[0],'H1':h[1],'H2':h[2],'H3':h[3],'H4':h[4],'H4_cond_H3_H4':h[4]/(h[3]+h[4]),'H3_plus_H4':h[3]+h[4]})
out=pathlib.Path('outputs/stage3_smim32_finngen_coloc_abf_sensitivity.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=res[0].keys());w.writeheader();w.writerows(res)
receipt={'method':'Single-causal-variant approximate Bayes factor / coloc-like H0–H4','n_variants':len(r),'inputs_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'outputs_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'implementation_check':'H3 differs from H4 configurations; exact tiny brute-force enumeration matches within 1e-13','qtl_prior_sd_values':[.075,.15,.3],'gwas_prior_sd':.2,'p1':.0001,'p2':.0001,'p12':[1e-6,1e-5,5e-5],'warning':'Exploratory custom ABF implementation; sample overlap and source model assumptions, linkage, window, variant density, phenotype scales, prior sensitivity and multiple causal signals not fully established; output not validated causal inference.'}
pathlib.Path('outputs/stage3_smim32_finngen_coloc_abf_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('variants',len(r));print(*res,sep='\n')
