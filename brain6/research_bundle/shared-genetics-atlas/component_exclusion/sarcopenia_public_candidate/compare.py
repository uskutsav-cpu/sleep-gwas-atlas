import csv,json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[1]
source=json.loads((p/'selected_source_rows.json').read_text());assert len(source)==7 and len({x['rsid'] for x in source})==7;source={r['rsid']:r for r in source}
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
leads=read(atlas/'replication/six_lead_joint/six_leads.tsv');sleep={x['variant']:x for x in read(atlas/'replication/mvp_six_leads/six_leads.tsv')};assert len(leads)==6 and len(sleep)==6
rows=[]
for x in leads:
 r=source[x['variant']];m=sleep[x['variant']];assert x['chr']==r['chromosome'] and x['position37']==r['base_pair_location'] and {x['effect_allele'],x['other_allele']}=={r['effect_allele'],r['other_allele']} and x['position38']==m['pos38'] and x['effect_allele']==m['effect_allele']
 a=float(m['mvp_p']);b=float(r['p_value']);assert math.isfinite(a) and math.isfinite(b) and 0<=a<=1 and 0<=b<=1
 rows.append({'variant':x['variant'],'chr':x['chr'],'position37':x['position37'],'mvp_insomnia_p':a,'sarcopenia_factor_p':b,'max_p':max(a,b),'max_p_times_six':min(1,6*max(a,b))})
with (p/'selected_comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print((p/'selected_comparison.tsv').read_text())
