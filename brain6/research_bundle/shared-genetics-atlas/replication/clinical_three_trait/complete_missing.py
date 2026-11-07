import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent; root=p.parent
source=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90707291.tsv'); adhd=root/'finngen_adhd_raw_lookup/indexed_results.tsv'; mvp=root/'mvp_six_leads/six_leads.tsv'; prior=p/'results.tsv'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
a=next(x for x in read(adhd) if x['variant']=='rs2582894'); b=next(x for x in read(mvp) if x['variant']=='rs2582894');assert a['status']=='MATCHED' and b['status']=='MATCHED'
h=hashlib.md5();found=[]; count=0;header=None
with source.open('rb') as f:
 for line in f:
  h.update(line)
  if header is None:header=line.decode().rstrip('\n').split('\t');continue
  count+=1
  if line.rstrip(b'\r\n').endswith((b'\trs2582894',b'\trs6452785')):found.append(dict(zip(header,line.decode().strip().split('\t'))))
assert h.hexdigest()=='38d9434f50a4493def9d2c9d1cc1b4b3' and count==21294560
rows=[x for x in found if x['rs_id']=='rs2582894'];assert len(rows)==1
x=rows[0];assert x['chromosome']==a['chr']==b['chr'] and x['base_pair_location']==a['position38']==b['pos38'] and {x['effect_allele'],x['other_allele']}=={a['ref'],a['alt']}=={a['effect_allele'],a['other_allele']} and b['effect_allele']==a['effect_allele']
ps=[float(b['mvp_p']),float(a['source_p']),float(x['p_value'])];assert all(math.isfinite(z) and 0<=z<=1 for z in ps)
alt=a['alt'];orient=lambda val,allele:float(val)*(1 if allele==alt else -1)
out={'variant':x['rs_id'],'grch38':f"chr{x['chromosome']}:{x['base_pair_location']}:{a['ref']}>{alt}",'mvp_insomnia_p':ps[0],'finngen_adhd_p':ps[1],'finngen_hfrs_p':ps[2],'mvp_beta_alt':orient(b['mvp_logor_aligned'],b['effect_allele']),'adhd_beta_alt':orient(a['beta_aligned'],a['effect_allele']),'hfrs_beta_alt':orient(x['beta'],x['effect_allele']),'max_p':max(ps),'max_p_times_7':min(1,7*max(ps))}
with (p/'rs2582894_completion.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out.keys(),delimiter='\t');w.writeheader();w.writerow(out)
(p/'completion_receipt.json').write_text(json.dumps({'source_md5':h.hexdigest(),'source_rows':count,'candidate_matches':found,'input_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (adhd,mvp,prior)},'output':out},indent=2)+'\n');print(json.dumps(out,indent=2))
