import pathlib,csv,math,json,hashlib
p=pathlib.Path(__file__).resolve().parent
def lsum(x):
 m=max(x);return m+math.log(math.fsum(math.exp(y-m) for y in x))
def bf(b,se,sd):
 r=sd*sd/(sd*sd+se*se);return .5*(math.log1p(-r)+r*(b/se)**2)
with (p/'posteriors.tsv').open() as f:results=list(csv.DictReader(f,delimiter='\t'))
maxerr=0
for row in results:
 assert row['prior']!='NOT_ESTIMABLE'
 with (p/(row['tissue']+'_'+row['gene']+'_'+row['trait']+'_matched.tsv')).open() as f:d=list(csv.DictReader(f,delimiter='\t'))
 prior=row['prior'];sd1={'default':.15,'narrow':.1,'wide':.3}[prior];sd2={'default':.2,'narrow':.1,'wide':.3}[prior]
 a=[bf(float(x['qtl_beta']),float(x['qtl_se']),sd1) for x in d];b=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),sd2) for x in d];la=lsum(a);lb=lsum(b);lab=lsum([x+y for x,y in zip(a,b)]);ld=la+lb+math.log1p(-math.exp(lab-la-lb));h=[0,math.log(1e-4)+la,math.log(1e-4)+lb,math.log(1e-8)+ld,math.log(float(row['p12']))+lab];den=lsum(h);q=[math.exp(x-den) for x in h]
 err=max(abs(q[i]-float(row['H'+str(i)])) for i in range(5));assert err<1e-10;maxerr=max(maxerr,err)
assert len(results)==48
with (p.parent/'gtex_v8/rs62264764_dynamic_summary.tsv').open() as f:dynamic=list(csv.DictReader(f,delimiter='\t'))
checks=[]
for tissue in ('Adrenal_Gland','Brain_Frontal_Cortex_BA9'):
 for gene in ('ENSG00000185565','ENSG00000239268'):
  with (p/(tissue+'_'+gene+'_insomnia_matched.tsv')).open() as f:a=[x for x in csv.DictReader(f,delimiter='\t') if x['variant']=='chr3_117920728_G_A']
  b=[x for x in dynamic if x['gencode_id'].split('.')[0]==gene and x['tissue']==tissue and x['status']=='OK'];assert len(a)==len(b)==1
  for s,t in [('qtl_beta','nes'),('qtl_se','se'),('qtl_p','nominal_p')]:assert math.isclose(float(a[0][s]),float(b[0][t]),rel_tol=1e-5,abs_tol=1e-9),(tissue,gene,s,a[0][s],b[0][t])
  checks.append({'tissue':tissue,'gene':gene,'p_beta_se':'match_within_source_rounding'})
output={'posterior_rows_verified':len(results),'maximum_absolute_difference':maxerr,'dynamic_gtex_comparisons':checks,'posterior_sha256':hashlib.sha256((p/'posteriors.tsv').read_bytes()).hexdigest(),'note':'arithmetic and source value checks; not model calibration'}
(p/'verification.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
