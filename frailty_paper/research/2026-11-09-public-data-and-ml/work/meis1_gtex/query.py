import json,hashlib,subprocess,urllib.parse
from pathlib import Path
root=Path('work/meis1_gtex');base='https://gtexportal.org/api/v2/';variant='chr2_66523432_G_T_b38'
queries={'variant':('dataset/variant',{'datasetId':'gtex_v10','variantId':variant}),'eqtl':('association/singleTissueEqtl',{'datasetId':'gtex_v10','variantId':variant,'itemsPerPage':250,'page':0}),'sqtl':('association/singleTissueSqtl',{'datasetId':'gtex_v10','variantId':variant,'itemsPerPage':250,'page':0})}
receipt={}
for name,(endpoint,args) in queries.items():
 url=base+endpoint+'?'+urllib.parse.urlencode(args);raw=subprocess.check_output(['curl','-fsSL','--max-time','30',url]);data=json.loads(raw);p=data['paging_info'];rows=data['data'];assert p['page']==0 and p['totalNumberOfItems']==len(rows) and p['numberOfPages']<=1
 if name=='variant':
  assert len(rows)==1 and rows[0]['snpId']=='rs113851554' and rows[0]['variantId']==variant and rows[0]['b37VariantId']=='2_66750564_G_T_b37'
 path=root/(name+'.json');path.write_bytes(raw);receipt[name]={'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'rows':len(rows),'paging':p}
Path('outputs/stage6_meis1_gtex_v10_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
