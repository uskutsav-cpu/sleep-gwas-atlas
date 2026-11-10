import csv,json
from pathlib import Path
studies=[('GSE117525','Muscle, frail/pre-frail vs healthy older',134,'outputs/stage6_gse117525_candidate_expression.json'),('GSE144304','Muscle, frail/pre-frail vs fit older',54,'outputs/stage6_gse144304_candidate_expression.json'),('GSE167186','Muscle, sarcopenia vs healthy older',53,'outputs/stage6_gse167186_candidate_expression.json'),('GSE40366','PBMC, frail vs non-frail nonagenarians',67,'outputs/stage6_gse40366_frailty_expression.json')]
rows=[]
for name,phen,n,path in studies:
 j=json.loads(Path(path).read_text())
 for r in j['results']:
  diff=r.get('log2cpm_adjusted_difference',r.get('adjusted_difference'));rows.append({'study':name,'outcome_contrast':phen,'n':n,'gene':r['gene'],'adjusted_difference':diff,'hc3_se':r['hc3_se'],'p_unadjusted':r['p_unadjusted'],'q_within_four_genes':r['bh_q_four_genes']})
with open('outputs/stage6_expression_four_study_summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('studies',len(studies),'rows',len(rows))
