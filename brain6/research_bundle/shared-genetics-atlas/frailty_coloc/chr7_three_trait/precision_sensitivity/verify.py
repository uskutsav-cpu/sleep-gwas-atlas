import csv,hashlib,json,math,pathlib,collections
p=pathlib.Path(__file__).resolve().parent
with (p/'replicates.tsv').open(newline='') as h:r=list(csv.DictReader(h,delimiter='\t'))
with (p.parent/'pairwise.tsv').open(newline='') as h:orig=list(csv.DictReader(h,delimiter='\t'))
assert len(r)==201 and [int(z['replicate']) for z in r]==list(range(201))
for name,pair,prior in [('FI_insomnia_1e5','FI_insomnia',1e-5),('FI_insomnia_1e6','FI_insomnia',1e-6),('FI_ADHD_1e5','FI_ADHD',1e-5),('FI_ADHD_1e6','FI_ADHD',1e-6)]:
 t=[q for q in orig if q['pair']==pair and float(q['p12'])==prior];assert len(t)==1 and math.isclose(float(r[0][name]),float(t[0]['H4']),abs_tol=1e-12)
for x in r:
 assert x['cluster']=='FI, insomnia' and 0<=float(x['posterior'])<=1
summary={'realizations':200,'baseline_matches_prior_pairwise':True,'cluster_counts':dict(collections.Counter(z['cluster'] for z in r[1:])),'candidate_counts':dict(collections.Counter(z['candidate'] for z in r[1:])), 'replicate_sha256':hashlib.sha256((p/'replicates.tsv').read_bytes()).hexdigest(),'source_input_sha256':hashlib.sha256(pathlib.Path('/Volumes/Brain6Work/chr7-three-trait-input.tsv').read_bytes()).hexdigest(),'note':'Consistency checks only; independent uniform jitter is an illustrative stress test, not an uncertainty distribution.'}
(p/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
