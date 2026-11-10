"""Retrieve the 15/22 union for BMI from the published 2022 KoGES PheWeb."""
import csv,json,urllib.request,re,hashlib,concurrent.futures,time
from pathlib import Path
root=Path('.')
files=['outputs/unvalidated_15_variant_score_prototype.csv','outputs/unvalidated_22_shared_direction_score.csv']
candidates={}
for file in files:
 with open(file,newline='') as f:
  for row in csv.DictReader(f):
   variant=row['variant'];chrom=row.get('chr') or row['chr_GRCh37_assumed'];pos=row.get('position_GRCh37_inferred') or row['position_GRCh37_assumed']; pair=frozenset([row['counted_allele'],row['other_allele']])
   if variant in candidates and (candidates[variant]['chrom'],candidates[variant]['pos'],candidates[variant]['pair'])!=(chrom,pos,pair):raise RuntimeError('inconsistent candidate '+variant)
   candidates[variant]={'chrom':chrom,'pos':pos,'pair':pair}
assert len(candidates)==33
requests=[]
for rs,x in candidates.items():
 a,b=sorted(x['pair']); requests += [(rs,x['chrom'],x['pos'],a,b),(rs,x['chrom'],x['pos'],b,a)]
def fetch(t):
 rs,ch,pos,ref,alt=t;url=f'https://koges.leelabsg.org/variant/{ch}:{pos}-{ref}-{alt}'
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'academic-genetics-research/1.0'}),timeout=30) as res: data=res.read().decode()
 except urllib.error.HTTPError as ex:
  if ex.code==404:return {'query':t,'state':'404'}
  raise
 m=re.search(r'window\.variant = (.*);',data)
 if not m:raise ValueError('missing variant JSON '+url)
 obj=json.loads(m.group(1));
 if (str(obj['chrom']),str(obj['pos']),obj['ref'],obj['alt'])!=(ch,pos,ref,alt):raise ValueError('identity mismatch')
 return {'query':t,'state':'ok','bmi':[p for p in obj['phenos'] if p['phenocode']=='KoGES_BMI']}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(fetch,requests))
Path('work/koges/bmi_candidate_raw.json').write_text(json.dumps(results,indent=2)+'\n')
rows=[]
for rs,x in candidates.items():
 matched=[z for z in results if z['query'][0]==rs and z.get('bmi')]
 if len(matched)>1:raise ValueError('ambiguous BMI allele representation '+rs)
 row={'variant':rs,'chr_GRCh37':x['chrom'],'pos_GRCh37':x['pos'],'allele_pair':'/'.join(sorted(x['pair'])),'status':'no exact BMI record'}
 if matched:
  z=matched[0];p=z['bmi'][0];row.update(status='exact pair found',ref=z['query'][3],alt=z['query'][4],p=p.get('pval'),beta_alt=p.get('beta'),se=p.get('sebeta'),af=p.get('af'),n=p.get('num_samples'))
 rows.append(row)
# Conservative multiplicity correction over all 33 candidates; unobserved variants treated as p=1.
rank=sorted(range(len(rows)), key=lambda i:float(rows[i].get('p',1)))
prev=1
for j in range(len(rank)-1,-1,-1):
 i=rank[j];q=min(prev,float(rows[i].get('p',1))*len(rows)/(j+1));rows[i]['BH_q_33']=q;prev=q
fields=['variant','chr_GRCh37','pos_GRCh37','allele_pair','status','ref','alt','beta_alt','se','af','n','p','BH_q_33']
with open('outputs/stage6_koges_bmi_candidate_followup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
print(json.dumps({'candidate_count':len(rows),'exact_matches':sum(x['status']=='exact pair found' for x in rows),'nominal_p_lt_005':sum(float(x.get('p',1))<0.05 for x in rows),'BH_q_lt_005':sum(x['BH_q_33']<0.05 for x in rows),'source':'https://koges.leelabsg.org/about','phenotype':'KoGES_BMI','raw_sha256':hashlib.sha256(Path('work/koges/bmi_candidate_raw.json').read_bytes()).hexdigest()},indent=2))
