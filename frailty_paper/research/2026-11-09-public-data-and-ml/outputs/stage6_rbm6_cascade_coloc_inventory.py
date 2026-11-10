#!/usr/bin/env python3
"""Query CASCADE public selected colocalizations near two postselected markers."""
from pathlib import Path
from urllib.parse import quote
import json,hashlib,requests,csv,collections
root=Path(__file__).resolve().parents[1];out=root/'outputs';regions=['chr3:49758613-50258613','chr3:49963404-50463404'];allrows={};receipt=[]
for region in regions:
 u='https://cascade.finngen.fi/api/cascade/coloc/by_region/'+quote(region,safe='');r=requests.get(u,timeout=45);r.raise_for_status();j=r.json();assert isinstance(j.get('pairs'),list);p=out/('stage6_cascade_'+region.split(':')[1]+'.json');p.write_text(json.dumps(j,separators=(',',':')));receipt.append({'region':region,'url':u,'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(j['pairs'])})
 for x in j['pairs']:
  key=json.dumps(x,sort_keys=True,separators=(',',':'));allrows[key]=x
pairs=list(allrows.values()); selected=[x for x in pairs if x.get('trait2_ensg')=='ENSG00000004534' and x.get('data_type2')=='eQTL'];columns=['trait1','cell_type2','PP.H4.abf','hit1','hit1_pos','hit2','hit2_pos','cs1_id','cs2_id','cs1_size','cs2_size','cs_overlap','clpp','qtl_beta','gwas_beta','dataset1','dataset2']
with open(out/'stage6_rbm6_cascade_selected_coloc.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows({k:x.get(k) for k in columns} for x in selected)
sleep=[x for x in pairs if any(y in str(x.get('trait1','')).lower() for y in ['insomn','sleep','frail'])];summary={'requests':receipt,'distinct_exact_pair_rows':len(pairs),'rbm6_eqtl_pairs':len(selected),'rbm6_trait_codes':dict(collections.Counter(x['trait1'] for x in selected)),'explicit_sleep_or_frailty_code_rows':len(sleep),'note':'browser-visible selected results only; trait code string check cannot prove phenotypes tested or absent'};(out/'stage6_rbm6_cascade_coloc_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
