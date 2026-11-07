import csv,json,pathlib,subprocess,urllib.parse,hashlib
out=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/variant_context/gtex_v8'); rec=[]
for gene in ['ENSG00000185565.11','ENSG00000239268.2']:
 for tissue in ['Brain_Cortex','Brain_Frontal_Cortex_BA9','Adrenal_Gland']:
  params={'datasetId':'gtex_v8','tissueSiteDetailId':tissue,'gencodeId':gene,'variantId':'chr3_117920728_G_A_b38'};url='https://gtexportal.org/api/v2/association/dyneqtl?'+urllib.parse.urlencode(params)
  try:
   result=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','45',url],check=True,capture_output=True); data=result.stdout;j=json.loads(data)
   assert j['gencodeId']==gene and j['tissueSiteDetailId']==tissue and j['variantId']==params['variantId']
   row={'status':'OK','gene':j.get('geneSymbol'),'gencode_id':gene,'tissue':tissue,'variant_id':j['variantId'],'nominal_p':j.get('pValue'),'nes':j.get('nes'),'se':j.get('error'),'t':j.get('tStatistic'),'n':sum(j[k] for k in ['hetCount','homoAltCount','homoRefCount']),'ref_homozygotes':j.get('homoRefCount'),'heterozygotes':j.get('hetCount'),'alt_homozygotes':j.get('homoAltCount'),'maf':j.get('maf'),'p_threshold':j.get('pValueThreshold'),'url':url,'response_sha256':hashlib.sha256(data).hexdigest()}
  except Exception as e:row={'status':'ERROR '+str(e),'gencode_id':gene,'tissue':tissue,'variant_id':params['variantId'],'url':url}
  rec.append(row);print(gene,tissue,row['status'],row.get('nominal_p'),row.get('nes'),flush=True)
fields=[]
for row in rec:
 for k in row:
  if k not in fields:fields.append(k)
with (out/'rs62264764_dynamic_summary.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rec)
