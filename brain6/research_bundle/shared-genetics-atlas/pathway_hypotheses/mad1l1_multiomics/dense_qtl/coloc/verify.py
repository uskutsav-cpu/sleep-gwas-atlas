import csv,math,pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent
with (p/'posteriors.tsv').open() as f:post=list(csv.DictReader(f,delimiter='\t'))
def lsum(x):
 m=max(x);return m+math.log(math.fsum(math.exp(y-m) for y in x))
def bf(beta,se,sd):
 v=se*se;r=sd*sd/(sd*sd+v);return .5*(math.log1p(-r)+r*(beta/se)**2)
maxerr=0;count=0
for row in post:
 if row['prior']=='NOT_ESTIMABLE':continue
 with (p/(row['gene']+'_'+row['trait']+'_matched.tsv')).open() as f:d=list(csv.DictReader(f,delimiter='\t'))
 prior=row['prior'];sd1={'default':.15,'narrow':.1,'wide':.3}[prior];sd2=({'default':.2,'narrow':.1,'wide':.3} if row['trait']=='insomnia' else {'default':.15,'narrow':.1,'wide':.3})[prior]
 a=[bf(float(x['qtl_beta']),float(x['qtl_se']),sd1) for x in d];b=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),sd2) for x in d];la=lsum(a);lb=lsum(b);lab=lsum([x+y for x,y in zip(a,b)]);ldiff=la+lb+math.log1p(-math.exp(lab-la-lb));h=[0,math.log(1e-4)+la,math.log(1e-4)+lb,math.log(1e-8)+ldiff,math.log(float(row['p12']))+lab];den=lsum(h);prob=[math.exp(v-den) for v in h]
 err=max(abs(v-float(row['H'+str(i)])) for i,v in enumerate(prob));maxerr=max(maxerr,err);count+=1;assert err<1e-10,(row['gene'],row['trait'],prior,err)
assert count==24
# Independently check displayed exact-marker values against prior Portal significant-pair records.
with (p.parent.parent/'gtex/associations.tsv').open() as f:gtex=list(csv.DictReader(f,delimiter='\t'))
checks=[]
for gene,id_,tissue in [('AC110781_3','ENSG00000176349','Brain_Nucleus_accumbens_basal_ganglia'),('MRM2','ENSG00000122687','Brain_Caudate_basal_ganglia')]:
 with (p/(gene+'_insomnia_matched.tsv')).open() as f:matched=list(csv.DictReader(f,delimiter='\t'))
 for rs in ('rs28528997','rs56305291'):
  a=[x for x in matched if x['rsid']==rs and x['alt']=='G'];b=[x for x in gtex if x['rsid']==rs and x['type']=='Eqtl' and x['gencode_id'].split('.')[0]==id_ and x['tissue']==tissue];assert len(a)==len(b)==1
  assert math.isclose(float(a[0]['qtl_p']),float(b[0]['p']),rel_tol=1e-8) and math.isclose(float(a[0]['qtl_beta']),float(b[0]['nes']),rel_tol=1e-8)
  checks.append([gene,rs,'p_and_beta_match'])
result={'posteriors_verified':count,'max_absolute_posterior_difference':maxerr,'portal_checks':checks,'posteriors_sha256':hashlib.sha256((p/'posteriors.tsv').read_bytes()).hexdigest(),'interpretation':'arithmetic and source-value verification only; not scientific calibration'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
