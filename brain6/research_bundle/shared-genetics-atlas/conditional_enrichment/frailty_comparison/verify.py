import csv,math,json,hashlib,pathlib
p=pathlib.Path(__file__).resolve().parent
with (p/'summary.tsv').open() as f:summary=list(csv.DictReader(f,delimiter='\t'))
with (p/'block_counts.tsv').open() as f:blocks=list(csv.DictReader(f,delimiter='\t'))
assert len(summary)==6
for s in summary:
 a,b=s['contrast'].split(' minus ');rows=[r for r in blocks if r['block_bp']==s['block_bp']];assert len(rows)==int(s['blocks'])
 n=sum(int(r['n']) for r in rows);na=sum(int(r['n_insomnia']) for r in rows);ca=sum(int(r[a+'_cond']) for r in rows);cb=sum(int(r[b+'_cond']) for r in rows);ja=sum(int(r[a+'_joint']) for r in rows);jb=sum(int(r[b+'_joint']) for r in rows)
 assert n==int(s['matched_n']) and ca==int(s['selected_a']) and cb==int(s['selected_b']) and ja==int(s['joint_a']) and jb==int(s['joint_b'])
 value=ja/ca-jb/cb;leave=[(ja-int(r[a+'_joint']))/(ca-int(r[a+'_cond']))-(jb-int(r[b+'_joint']))/(cb-int(r[b+'_cond'])) for r in rows];mean=sum(leave)/len(leave);se=math.sqrt((len(leave)-1)/len(leave)*sum((v-mean)**2 for v in leave))
 assert math.isclose(value,float(s['difference']),rel_tol=1e-12);assert math.isclose(se,float(s['jackknife_se']),rel_tol=1e-12);assert math.isclose(na/n,float(s['baseline_insomnia']),rel_tol=1e-12)
root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation');ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  while q:=z.read(4*1024*1024):h.update(q)
 return h.hexdigest()
files={k:root/'common_support'/f'{k}.sumstats.gz' for k in ('S','FI','FFS')};files['HFRS']=root/'filtered/HFRS_HM3.tsv';files['reference']=ref
j={'files':{k:{'path':str(v),'sha256':sha(v)} for k,v in files.items()},'output_sha256':{'summary':sha(p/'summary.tsv'),'blocks':sha(p/'block_counts.tsv'),'script':sha(p/'run.R')},'comparisons_checked':len(summary),'crosscheck':'PASS paired block count, fraction, difference and jackknife SE reproducibility','caveat':'Does not establish block independence or genetic inference validity.'}
(p/'verification.json').write_text(json.dumps(j,indent=2)+'\n');print(j['crosscheck'])
