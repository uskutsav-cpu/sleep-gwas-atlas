import json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent; root=p.parents[3]
original=json.loads((p/'original_rows.json').read_text());assert len(original)==7 and len({x['SNP_ID'] for x in original})==7
m=json.loads((root/'mvp_six_leads/receipt.json').read_text());lookup={x['variant']:x['mvp_exact_rows'][0] for x in m['audit']};q=json.loads((root/'rs62264764_mvp/query_receipt.json').read_text());lookup['rs62264764']=q['exact_records'][0]
rows=[]
for r in original:
 x=lookup[r['SNP_ID']];checks={'chr':r['chrom']==x['chromosome'],'pos':r['pos']==x['base_pair_location'],'effect_allele':r['ea']==x['effect_allele'],'alleles':{r['ref'],r['alt']}=={x['effect_allele'],x['other_allele']},'or':r['or']==x['odds_ratio'],'ci':r['ci']==x['ci_lower']+','+x['ci_upper'],'p':r['pval']==x['p_value'],'cases':r['num_cases']==x['num_cases'],'controls':r['num_controls']==x['num_controls']};assert all(checks.values()),(r['SNP_ID'],checks);rows.append({'variant':r['SNP_ID'],'all_checked_fields_match':True})
(p/'comparison.json').write_text(json.dumps({'comparison':rows,'fields':['chr','pos','effect_allele','alleles','or','ci','p','cases','controls'],'source_rows_sha256':hashlib.sha256((p/'original_rows.json').read_bytes()).hexdigest(),'scope':'seven selected variants only; full-file identity not established'},indent=2)+'\n');print('seven variants, all selected fields match')
