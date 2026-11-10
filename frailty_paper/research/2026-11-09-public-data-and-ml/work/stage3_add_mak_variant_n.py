import csv,pathlib,collections,json,hashlib
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort'); summary={}
for src,name in [('fi.fastGWA','Mak_FI'),('fp.fastGWA','Mak_FP')]:
 p=pathlib.Path('work/stage3_smim32_'+name+'_aligned.csv');rows=list(csv.DictReader(p.open()));ids={r['variant'] for r in rows};found=collections.defaultdict(list)
 with (root/src).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['SNP'] in ids:found[r['SNP']].append(r)
 for r in rows:
  hits=found[r['variant']];assert len(hits)==1;u=hits[0];assert {u['A1'],u['A2']}=={r['ref'],r['alt']};r['gwas_n']=int(u['N']);assert r['gwas_n']>100000
 out=pathlib.Path('work/stage3_smim32_'+name+'_with_n.csv')
 with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 summary[name]={'n_variants':len(rows),'n_min':min(x['gwas_n'] for x in rows),'n_max':max(x['gwas_n'] for x in rows),'sha256':hashlib.sha256(out.read_bytes()).hexdigest()};print(name,summary[name],flush=True)
pathlib.Path('outputs/stage3_mak_variant_sample_sizes.json').write_text(json.dumps(summary,indent=2)+'\n')
