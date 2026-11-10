#!/usr/bin/env python3
"""Stream public PGS Catalog files; retain only candidate-level matches."""
import csv, gzip, io, subprocess
from pathlib import Path
p=Path('outputs'); models={}
for label in ('15','22'):
 with (p / f'unvalidated_{label}_variant_score_prototype.csv').open() as f:
  for r in csv.DictReader(f): models.setdefault(r['variant'],[]).append((label,r))
positions={}
for fn in ('ncds_15_grch38_reference_audit.csv','ncds_22_grch38_ensembl_audit.csv'):
 with (p/fn).open() as f:
  for row in csv.DictReader(f): positions[row['variant']]=(row['chr38'],row['pos38'])
by_pos = {}
for rs,pos in positions.items(): by_pos.setdefault(pos, set()).add(rs)
output=[]
for pgs in ('PGS019928','PGS019929'):
 url=f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs}/ScoringFiles/{pgs}.txt.gz'
 proc = subprocess.Popen(['curl','-fsSL','--max-time','90',url], stdout=subprocess.PIPE)
 with gzip.GzipFile(fileobj=proc.stdout) as z, io.TextIOWrapper(z) as stream:
  d=csv.DictReader((line for line in stream if not line.startswith('#')),delimiter='\t')
  n=0
  for r in d:
   n+=1
   matches=set(by_pos.get((r['chr_name'],r['chr_position']),set()))
   if r.get('rsID') in models: matches.add(r['rsID'])
   for rs in sorted(matches):
    for label,m in models[rs]:
     pos=positions.get(rs); pair={r['effect_allele'],r['other_allele']}; modelpair={m['counted_allele'],m['other_allele']}
     output.append({'pgs':pgs,'model':label,'variant':rs,'pgs_rsid':r.get('rsID',''),'pgs_chr38':r['chr_name'],'pgs_pos38':r['chr_position'],'position_matches_project_derived':str(pos==(r['chr_name'],r['chr_position'])),'allele_pair_match':str(pair==modelpair),'pgs_effect_allele':r['effect_allele'],'pgs_other_allele':r['other_allele'],'pgs_beta':r['effect_weight'],'project_counted':m['counted_allele'],'project_other':m['other_allele']})
  print(pgs,'records',n,'candidate_rows',sum(x['pgs']==pgs for x in output))
 if proc.wait() != 0: raise RuntimeError('Download failed: '+pgs)
with (p/'stage6_wyss_pgs_candidate_overlap.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
