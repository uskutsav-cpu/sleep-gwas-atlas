import csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
p=Path(__file__).resolve().parent;follow=p.parent/'results.tsv';source=p.parents[3]/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
r=read(follow);s=read(source);by={x['snp']:x for x in s};assert len(r)==len(by)==425
out=[]
for x in r:
 y=by[x['rsid']]; assert x['grch38']==y['atlas_variant'] and x['region']==y['candidate_locus_id']; ref,alt=x['grch38'].split(':')[-1].split('>'); assert {y['source_a1'],y['source_a2']}=={ref,alt}
 row={'rsid':x['rsid'],'region':x['region'],'grch38':x['grch38'],'status':x['adhd_status'],'discovery_beta_alt':'','followup_beta_alt':x['adhd_beta_alt'],'followup_se':x['adhd_se'],'followup_p':x['adhd_p'],'p_from_beta_se':'','abs_log10_p_difference':'','same_direction':''}
 if x['adhd_status']=='matched':
  b=float(x['adhd_beta_alt']);se=float(x['adhd_se']);pv=float(x['adhd_p']);assert math.isfinite(b) and se>0 and 0<pv<=1
  a=float(y['adhd_beta_a1'])*(1 if y['source_a1']==alt else -1) if y['adhd_beta_a1'] else None
  pp=math.erfc(abs(b/se)/math.sqrt(2));assert pp>0
  row.update(p_from_beta_se=pp,abs_log10_p_difference=abs(math.log10(pp)-math.log10(pv)))
  if a is not None and math.isfinite(a) and a!=0 and b!=0: row.update(discovery_beta_alt=a,same_direction=(a*b>0))
 out.append(row)
matched=[x for x in out if x['same_direction']!='']; arithmetic=[x for x in out if x['abs_log10_p_difference']!=''];assert len(arithmetic)==420
regions=[]
for region in sorted(set(x['region'] for x in out)):
 z=[x for x in matched if x['region']==region];regions.append({'region':region,'matched':len(z),'same_direction':sum(x['same_direction'] for x in z),'nominal_p_lt_005':sum(float(x['followup_p'])<.05 for x in z),'same_direction_and_nominal':sum(x['same_direction'] and float(x['followup_p'])<.05 for x in z)})
def write(name,z):
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=z[0],delimiter='\t');w.writeheader();w.writerows(z)
write('candidate_directions.tsv',out);write('region_summary.tsv',regions)
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (follow,source)},'matched':len(matched),'same_direction':sum(x['same_direction'] for x in matched),'max_abs_log10_p_difference':max(x['abs_log10_p_difference'] for x in arithmetic),'warnings':'correlated, previously selected SNPs; p reconstructed via normal approximation'},indent=2)+'\n')
print(json.dumps(regions,indent=2));print('max log10 difference',max(x['abs_log10_p_difference'] for x in arithmetic))
