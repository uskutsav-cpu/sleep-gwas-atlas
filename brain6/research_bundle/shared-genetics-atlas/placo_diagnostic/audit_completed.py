#!/usr/bin/env python3
"""Audit completed PLACO+ outputs only; never consume a live partial file."""
import argparse,collections,csv,gzip,hashlib,json,math,pathlib
ap=argparse.ArgumentParser();ap.add_argument('--input',type=pathlib.Path,required=True);ap.add_argument('--results',type=pathlib.Path,required=True);ap.add_argument('--ffs-results',type=pathlib.Path,default=None);ap.add_argument('--out',type=pathlib.Path,required=True);ap.add_argument('--expected-n',type=int,default=1138262);ap.add_argument('--controls',type=pathlib.Path,required=True);ap.add_argument('--pairs',nargs='+',choices=('FI','FFS'),default=['FI','FFS']);a=ap.parse_args();assert len(set(a.pairs))==len(a.pairs);assert a.expected_n>0
family=2*a.expected_n;threshold=.05/family;all_hits=[];meta={'expected_rows_per_pair':a.expected_n,'family_size':family,'bonferroni_threshold':threshold,'pair_counts':{},'files_sha256':{}}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
with a.controls.open() as f:control_rows=list(csv.DictReader(f,delimiter='\t'))
controls={r['SNP']:r for r in control_rows};assert len(controls)==len(control_rows) and len(controls)>0
meta['control_file_sha256']=sha(a.controls);meta['n_control_markers']=len(controls)
for pair in a.pairs:
 sf=a.input/'S.sumstats.gz';tf=a.input/(pair+'.sumstats.gz');rf=(a.ffs_results if pair=='FFS' and a.ffs_results else a.results)/(pair+'_screen.tsv.gz')
 if not rf.is_file():raise FileNotFoundError(f'Completed output absent: {rf}; no partial results admitted')
 for f in (sf,tf,rf):meta['files_sha256'][str(f)]=sha(f)
 counts=collections.Counter();seen=set()
 with gzip.open(sf,'rt') as s,gzip.open(tf,'rt') as t,gzip.open(rf,'rt') as r:
  for sr,tr,rr in zip(*(csv.DictReader(x,delimiter='\t') for x in (s,t,r)),strict=True):
   k=tuple(sr[x] for x in ('SNP','A1','A2'));assert all(tuple(v[x] for x in ('SNP','A1','A2'))==k for v in (tr,rr));assert k[0] not in seen;seen.add(k[0])
   z1,z2=float(sr['Z']),float(tr['Z']);assert z1==float(rr['Z_insomnia']) and z2==float(rr['Z_frailty']); assert all(math.isfinite(z) for z in (z1,z2))
   p1,p2=(math.erfc(abs(z)/math.sqrt(2)) for z in (z1,z2));assert math.isclose(p1,float(rr['P_insomnia']),rel_tol=2e-12,abs_tol=1e-300);assert math.isclose(p2,float(rr['P_frailty']),rel_tol=2e-12,abs_tol=1e-300)
   status=rr['status'];assert status in ('TESTED','EXCLUDED_EXTREME_Z','NUMERICAL_FAILURE');counts[status]+=1;counts['rows']+=1
   if status=='TESTED':
    pp=float(rr['P_PLACO_plus']);assert math.isfinite(pp) and 0<pp<=1;assert z1*z1<=80 and z2*z2<=80
    if k[0] in controls:
     ref=controls[k[0]];assert tuple(ref[x] for x in ('SNP','A1','A2'))==k;assert float(ref['Z_insomnia'])==z1 and float(ref['Z_'+pair])==z2
     assert math.isclose(pp,float(ref['P_PLACO_'+pair]),rel_tol=2e-12,abs_tol=1e-300);counts['control_matches']+=1
    if pp<=threshold:
     counts['threshold_crossing']+=1;conj=max(p1,p2);all_hits.append({'pair':pair,'SNP':k[0],'A1':k[1],'A2':k[2],'Z_insomnia':z1,'Z_frailty':z2,'p_insomnia':p1,'p_frailty':p2,'p_placo_plus':pp,'p_placo_bonferroni':min(1,pp*family),'p_max_marginal':conj,'p_max_marginal_bonferroni':min(1,conj*family)})
   else:
    assert rr['P_PLACO_plus']=='NA'
    if status=='EXCLUDED_EXTREME_Z':assert z1*z1>80 or z2*z2>80
 assert counts['rows']==a.expected_n and counts['control_matches']==len(controls);meta['pair_counts'][pair]=dict(counts)
a.out.mkdir(exist_ok=True,parents=True);fields=['pair','SNP','A1','A2','Z_insomnia','Z_frailty','p_insomnia','p_frailty','p_placo_plus','p_placo_bonferroni','p_max_marginal','p_max_marginal_bonferroni'];path=a.out/'audited_threshold_crossings.tsv'
with path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(sorted(all_hits,key=lambda x:x['p_placo_plus']))
meta['hits_sha256']=sha(path);(a.out/'completed_audit.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
