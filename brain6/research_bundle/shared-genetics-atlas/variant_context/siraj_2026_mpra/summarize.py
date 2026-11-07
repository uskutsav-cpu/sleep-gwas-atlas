import json,csv,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent
records=json.loads((p/'sknsh_matches.json').read_text());rows=[]
for rec in records:
 r=rec['row']; f=lambda key:'' if r.get(key) in (None,'NA','') else 10**(-float(r[key]))
 rows.append({'rsid':rec['rsid'],'variant37':rec['query_variant'],'library':rec['source'].split('/')[-1].split('_')[0],'log2_skew':r['Log2Skew'],'skew_se':r['Skew_SE'],'p':f('Skew_logP'),'source_fdr':f('Skew_logFDR'),'source_fdr_activity_subset':f('Skew_logFDR_act'),'a_log2_activity':r['A_log2FC'],'b_log2_activity':r['B_log2FC'],'a_activity_bonferroni_p':f('A_logPadj_BF'),'b_activity_bonferroni_p':f('B_logPadj_BF')})
assert len(rows)==38 and len({r['rsid'] for r in rows})==32
with (p/'sknsh_results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('rs1456193',next(r for r in rows if r['rsid']=='rs1456193'));print('rs12658032',[r for r in rows if r['rsid']=='rs12658032'])
