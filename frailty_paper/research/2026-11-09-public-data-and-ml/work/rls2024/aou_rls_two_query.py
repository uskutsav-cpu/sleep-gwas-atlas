import csv,json,hashlib,subprocess,urllib.parse
from pathlib import Path
O=Path('outputs');W=Path('work/rls2024');meta=Path('work/aou_complete/analysis_metadata_check.json');m=json.loads(meta.read_text());entry=[x for x in m if x['analysis_id']=='NS_333.64' and x['ancestry_group']=='EUR'];assert len(entry)==1
queries=[('rs113851554',66523432,'G','T'),('rs78206187',66795102,'A','G')];rows=[];receipts=[]
for rs,pos,ref,alt in queries:
 q=urllib.parse.urlencode({'analysis_id':'NS_333.64','ancestry_group':'eur','sequencing_type':'genome','query_mode':'slow'});url=f'https://allbyall.researchallofus.org/api/variants/associations/interval/chr2:{pos-1}-{pos+1}?{q}'
 raw=subprocess.check_output(['curl','-fsSL','--max-time','30',url]);j=json.loads(raw);assert j['storage_source']=='hail_gcs';vid=f'chr2-{pos}-{ref}-{alt}';hits=[r for r in j['data'] if r['variant_id']==vid];assert len(hits)==1
 h=hits[0];assert (h['phenotype'],h['ancestry'],h['sequencing_type'])==('NS_333.64','eur','genome') and h['ref']==ref and h['alt']==alt
 rows.append({'candidate':rs,'variant_id':vid,'phenotype':h['phenotype'],'effect_assumed':'ALT','beta':h['beta'],'se':h['se'],'pvalue':h['pvalue'],'af':h['af'],'n_returned_at_position':j['count']});(W/(rs+'_aou_rls.json')).write_bytes(raw);receipts.append({'url':url,'raw_sha256':hashlib.sha256(raw).hexdigest()})
p=O/'stage6_aou_rls_two_candidate_followup.csv'
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
r={'metadata':entry[0],'metadata_sha256':hashlib.sha256(meta.read_bytes()).hexdigest(),'queries':receipts,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'caveat':'Current public CDRv8 EUR group; ALT direction supported by pipeline audit; not same definition, period or sample as 2024 Akcimen study; related overlap uncertain; no individual score or validation.'}
(O/'stage6_aou_rls_two_candidate_followup_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(rows,indent=2))
