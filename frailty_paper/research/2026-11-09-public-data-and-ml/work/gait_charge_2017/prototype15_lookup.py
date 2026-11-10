import csv,hashlib,json,requests
from pathlib import Path
src='https://www.nhlbi.nih.gov/grasp/downloads/FullResults/2017/2017_Ben_Avrahim/Gait-Aging-All.txt'
score=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
lookup={x['variant']:[] for x in score}; digest=hashlib.sha256(); count=0
with requests.get(src,stream=True,timeout=(30,180)) as response:
 response.raise_for_status(); stream=response.iter_lines(); head=next(stream);digest.update(head+b'\n');names=head.decode().split()
 assert names[0]=='all' and names[2:4]==['Allele1','Allele2']
 for line in stream:
  if not line: continue
  digest.update(line+b'\n'); cells=line.decode().split(); assert len(cells)==len(names);count+=1
  if cells[0] in lookup:lookup[cells[0]].append(dict(zip(names,cells)))
rows=[]
for item in score:
 a,b=item['counted_allele'],item['other_allele']; raw=lookup[item['variant']]; matches=[r for r in raw if {r['Allele1'].upper(),r['Allele2'].upper()}=={a.upper(),b.upper()}]
 result={'variant':item['variant'],'candidate_allele':a,'other_allele':b,'rsid_matches':len(raw),'pair_matches':len(matches),'palindromic':{a,b} in [{'A','T'},{'C','G'}]}
 if len(matches)==1:
  r=matches[0];result.update({'source_allele1':r['Allele1'],'source_allele2':r['Allele2'],'source_beta':r['Effect'],'source_allele1_frequency':r['Freq1'],'se':r['StdErr'],'p':r['P-value'],'heterogeneity_p':r['HetPVal']})
 rows.append(result)
fields=['variant','candidate_allele','other_allele','rsid_matches','pair_matches','palindromic','source_allele1','source_allele2','source_beta','source_allele1_frequency','se','p','heterogeneity_p']
with open('outputs/stage6_gait_15_candidate_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
receipt={'source':src,'source_rows':count,'sha256_lf_normalized':digest.hexdigest(),'score_sha256':hashlib.sha256(Path('outputs/unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'ids':len(score),'unique_pair_matches':sum(x['pair_matches']==1 for x in rows),'nominal_p_below_005':sum(float(x['p'])<.05 for x in rows if 'p' in x),'caution':'No source coordinates; no genetic score association; no allele frequency or strand verification.'}
Path('outputs/stage6_gait_15_candidate_lookup.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
