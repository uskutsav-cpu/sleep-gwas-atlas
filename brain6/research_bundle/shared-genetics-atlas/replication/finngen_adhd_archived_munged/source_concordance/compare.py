#!/usr/bin/env python3
import collections,csv,hashlib,json,math,re
from pathlib import Path
from statistics import NormalDist
p=Path(__file__).resolve().parent
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
sample=read(p/'sample.tsv');source=read(p/'raw_matches.tsv');lookup=collections.defaultdict(list)
for r in source:
 for rid in re.split(r'[,;]',r['rsids']):lookup[rid].append(r)
out=[]
for r in sample:
 z={'SNP':r['SNP'],'sample_Z':r['Z'],'status':'','raw_p':'','raw_beta_over_se_aligned':'','raw_p_to_z_aligned':'','difference_pz':'','difference_wald':''}
 alleles={r['A1'],r['A2']}
 if len(alleles)!=2 or not alleles.issubset(set('ACGT')) or alleles in ({'A','T'},{'C','G'}):z['status']='ambiguous_or_invalid_alleles'
 else:
  matches=[x for x in lookup[r['SNP']] if {x['ref'],x['alt']}==alleles]
  if len(matches)!=1:z['status']='missing_match' if not matches else 'multiple_matches'
  else:
   x=matches[0]
   try:
    beta=float(x['beta']);se=float(x['sebeta']);pv=float(x['pval']);az=float(r['Z']);assert 0<pv<=1 and se>0 and all(map(math.isfinite,(beta,se,pv,az))) and beta!=0
    sign=(1 if beta>0 else -1)*(1 if x['alt']==r['A1'] else -1)
    zp=-NormalDist().inv_cdf(pv/2)*sign;zw=beta/se*(1 if x['alt']==r['A1'] else -1)
    z.update(status='matched',raw_p=pv,raw_beta_over_se_aligned=zw,raw_p_to_z_aligned=zp,difference_pz=az-zp,difference_wald=az-zw)
   except (ValueError,AssertionError):z['status']='invalid_statistics'
 out.append(z)
with (p/'comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
matched=[x for x in out if x['status']=='matched'];d1=sorted(abs(float(x['difference_pz'])) for x in matched);d2=sorted(abs(float(x['difference_wald'])) for x in matched)
assert matched
result={'counts':dict(collections.Counter(x['status'] for x in out)),'raw_matched_records':len(source),'abs_difference_p_to_z':{'median':d1[len(d1)//2],'p95':d1[int(.95*(len(d1)-1))],'max':max(d1),'within_1e-4':sum(x<=1e-4 for x in d1)},'abs_difference_beta_over_se':{'median':d2[len(d2)//2],'p95':d2[int(.95*(len(d2)-1))],'max':max(d2),'within_1e-4':sum(x<=1e-4 for x in d2)},'input_hashes':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('sample.tsv','raw_matches.tsv')},'limitations':'Cross-file numerical consistency, not an ingest record or cohort independence audit.'}
(p/'comparison_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
