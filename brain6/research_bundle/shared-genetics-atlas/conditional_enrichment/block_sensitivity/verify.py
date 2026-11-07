import csv,collections,math,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent
with (p/'summary.tsv').open() as f:s=list(csv.DictReader(f,delimiter='\t'))
with (p/'block_counts.tsv').open() as f:b=list(csv.DictReader(f,delimiter='\t'))
assert len(s)==8
for x in s:
 rows=[r for r in b if (r['pair'],r['direction'],r['block_bp'])==(x['pair'],x['direction'],x['block_bp'])];assert len(rows)==int(x['n_blocks'])
 n,t,c,j=[sum(int(r[k]) for r in rows) for k in ('n','n_target','n_condition','n_joint')]; assert n==int(x['n_matched'])
 delta=j/c-t/n;leave=[(j-int(r['n_joint']))/(c-int(r['n_condition']))-(t-int(r['n_target']))/(n-int(r['n'])) for r in rows];mean=sum(leave)/len(leave);se=math.sqrt((len(leave)-1)/len(leave)*sum((v-mean)**2 for v in leave))
 assert math.isclose(delta,float(x['difference']),rel_tol=1e-12);assert math.isclose(se,float(x['jackknife_se']),rel_tol=1e-12)
ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  while q:=z.read(4*1024*1024):h.update(q)
 return h.hexdigest()
out={'reference_path':str(ref),'reference_sha256':sha(ref),'source_receipt_sha256':sha(p.parent/'receipt.json'),'summary_sha256':sha(p/'summary.tsv'),'block_counts_sha256':sha(p/'block_counts.tsv'),'script_sha256':sha(p/'run.R'),'checked_combinations':len(s),'result':'PASS block counts reproduce differences and jackknife SEs','caveat':'Does not establish physical block independence or calibration.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(out['result'])
