import json,hashlib,subprocess,urllib.parse
from pathlib import Path
root=Path('work/meis1_gtex');base='https://gtexportal.org/api/v2/'
variants={'rs13033745':'chr2_66433271_G_T_b38','rs1519102':'chr2_66450684_G_C_b38'}
receipt={}
for rs,variant in variants.items():
 receipt[rs]={}
 for kind,endpoint in [('variant','dataset/variant'),('eqtl','association/singleTissueEqtl'),('sqtl','association/singleTissueSqtl')]:
  args={'datasetId':'gtex_v10','variantId':variant,'itemsPerPage':250,'page':0};url=base+endpoint+'?'+urllib.parse.urlencode(args)
  try:
   raw=subprocess.check_output(['curl','-fsSL','--max-time','25',url]);obj=json.loads(raw);info=obj['paging_info'];rows=obj['data'];assert info['page']==0 and info['totalNumberOfItems']==len(rows) and info['numberOfPages']<=1
   if kind=='variant':assert len(rows)==1 and rows[0]['snpId']==rs and rows[0]['variantId']==variant
   path=root/(rs+'_'+kind+'.json');path.write_bytes(raw);receipt[rs][kind]={'url':url,'rows':len(rows),'sha256':hashlib.sha256(raw).hexdigest(),'status':'ok'}
  except Exception as e:receipt[rs][kind]={'url':url,'status':'error','error':str(e)}
Path('outputs/stage6_meis1_proxy_gtex_v10_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
