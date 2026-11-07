import pathlib,csv,math,json,hashlib
p=pathlib.Path(__file__).resolve().parent
def lsum(x):
 m=max(x);return m+math.log(math.fsum(math.exp(y-m) for y in x))
def bf(b,se,sd):
 r=sd*sd/(sd*sd+se*se);return .5*(math.log1p(-r)+r*(b/se)**2)
with (p/'posteriors.tsv').open() as f:results=list(csv.DictReader(f,delimiter='\t'))
maxerr=0
for row in results:
 if row['prior']=='NOT_ESTIMABLE':
  assert row['n']=='0' and row['H4']=='NA'; continue
 with (p/(row['tissue']+'_'+row['gene']+'_'+row['trait']+'_matched.tsv')).open() as f:d=list(csv.DictReader(f,delimiter='\t'))
 prior=row['prior'];sd1={'default':.15,'narrow':.1,'wide':.3}[prior];sd2={'default':.2,'narrow':.1,'wide':.3}[prior]
 a=[bf(float(x['qtl_beta']),float(x['qtl_se']),sd1) for x in d];b=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),sd2) for x in d];la=lsum(a);lb=lsum(b);lab=lsum([x+y for x,y in zip(a,b)]);ld=la+lb+math.log1p(-math.exp(lab-la-lb));h=[0,math.log(1e-4)+la,math.log(1e-4)+lb,math.log(1e-8)+ld,math.log(float(row['p12']))+lab];den=lsum(h);q=[math.exp(x-den) for x in h]
 err=max(abs(q[i]-float(row['H'+str(i)])) for i in range(5));assert err<1e-10;maxerr=max(maxerr,err)
assert len(results)==14

with (p.parent/'gtex_v8/significant_associations.tsv').open() as f: static=list(csv.DictReader(f,delimiter='\t'))
with (p/'Testis_ENSG00000251574_insomnia_matched.tsv').open() as f: matched=list(csv.DictReader(f,delimiter='\t'))
checks=[]
for rsid in ('rs2431108','rs77960'):
 a=[x for x in matched if x['gwas_rsid']==rsid]; b=[x for x in static if x['snp']==rsid and x['type']=='Eqtl' and x['tissue']=='Testis']
 assert len(a)==len(b)==1,(rsid,len(a),len(b))
 for q,t in [('qtl_beta','normalized_effect_size'),('qtl_p','p_value')]: assert math.isclose(float(a[0][q]),float(b[0][t]),rel_tol=1e-5,abs_tol=1e-10)
 checks.append(rsid)
output={'posterior_rows_verified':12,'not_estimable_rows':2,'maximum_absolute_difference':maxerr,'static_significant_pair_checks':checks,'posterior_sha256':hashlib.sha256((p/'posteriors.tsv').read_bytes()).hexdigest(),'note':'arithmetic and source value checks; not model validation'}
(p/'verification.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
