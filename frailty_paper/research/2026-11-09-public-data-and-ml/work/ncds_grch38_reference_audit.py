import csv, json, hashlib
from pathlib import Path
p=Path('outputs/ncds_score_variant_request.csv'); q=Path('work/variant_reference_grch38.csv')
a=[r for r in csv.DictReader(p.open()) if r['score']=='15_weighted']; e={r['variant']:r for r in csv.DictReader(q.open())}
rows=[]
for r in a:
 x=e[r['variant']]; ms=json.loads(x['Ensembl_GRCh38_mappings']); same=[m for m in ms if str(m['chr'])==r['chr38'] and str(m['pos'])==r['pos38'] and m['assembly']=='GRCh38']; pair={r['counted_allele_38_candidate'],r['other_allele_38_candidate']}; alleles=set(same[0]['alleles'].split('/')) if len(same)==1 else set()
 rows.append({'variant':r['variant'],'chr38':r['chr38'],'pos38':r['pos38'],'same_position_unique':len(same)==1,'reference_alleles':same[0]['alleles'] if len(same)==1 else '', 'both_score_alleles_listed':pair<=alleles,'reference_multiallelic':len(alleles)>2,'source_url':x['url'],'saved_response_sha256':x['sha256']})
o=Path('outputs/ncds_15_grch38_reference_audit.csv')
with o.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps({'rows':len(rows),'unique_location':sum(x['same_position_unique'] for x in rows),'both_alleles_listed':sum(x['both_score_alleles_listed'] for x in rows),'reference_multiallelic':sum(x['reference_multiallelic'] for x in rows),'reference_inventory_sha256':hashlib.sha256(q.read_bytes()).hexdigest()},indent=2))
