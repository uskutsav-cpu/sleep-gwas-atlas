#!/usr/bin/env python3
"""Check the manually transcribed TPMI positions against Ensembl mappings."""
import hashlib,json,time
from pathlib import Path
import pandas as pd, requests
O=Path('outputs');p=O/'stage6_tpmi_sleep_candidate_inventory.csv';d=pd.read_csv(p);out=[]
for r in d.itertuples(index=False):
 if r.status!='browser_lookup':out.append({'variant':r.variant,'status':'not_checked'});continue
 u=f'https://grch37.rest.ensembl.org/map/human/GRCh37/{r.prototype_chr37}:{r.prototype_pos37}..{r.prototype_pos37}/GRCh38?content-type=application/json'
 try:
  resp=requests.get(u,timeout=20);resp.raise_for_status();m=resp.json().get('mappings',[])
  ok=len(m)==1 and str(m[0]['mapped']['seq_region_name'])==str(int(r.chr38)) and m[0]['mapped']['start']==int(r.pos38) and m[0]['mapped']['strand']==1
  out.append({'variant':r.variant,'status':'unique_positive_match' if ok else 'unresolved','url':u,'response_sha256':hashlib.sha256(resp.content).hexdigest(),'n_mappings':len(m)})
 except (requests.RequestException,ValueError,KeyError) as e:out.append({'variant':r.variant,'status':'lookup_error','error':type(e).__name__,'url':u})
 time.sleep(.12)
a=pd.DataFrame(out);a.to_csv(O/'stage6_tpmi_liftover_audit.csv',index=False);d=d.drop(columns=['ensembl_liftover_check'],errors='ignore').merge(a[['variant','status']].rename(columns={'status':'ensembl_liftover_check'}),on='variant',validate='one_to_one');d.to_csv(p,index=False)
receipt={'source':'Ensembl GRCh37 assembly mapping API','results':d.ensembl_liftover_check.value_counts().to_dict(),'scope':'Position mapping only; effect allele, source data quality and score validation are not established.'}
(O/'stage6_tpmi_liftover_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
