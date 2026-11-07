#!/usr/bin/env python3
import csv,hashlib,json,math,subprocess
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;base=Path('/Volumes/Brain6Work/xqtl-chr5')
expected={'ast.bed.gz':'0773eafce97111a3c5f995f6c12cc28a','ast.bed.gz.tbi':'768aae3a32c8c90abc04afb8c4622997','opc.bed.gz':'5efae548fea817d4ffa783f9442bcbfc','opc.bed.gz.tbi':'f777b9c92d047a0e0a5595f58e7500b6'}
def digest(path,algo):
 h=hashlib.new(algo)
 with path.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
hashes={}
for name,md5 in expected.items():
 actual=digest(base/name,'md5');assert actual==md5;hashes[name]={'md5':actual,'sha256':digest(base/name,'sha256')}
variants={'rs62369151':(88587685,'T','C'),'rs4916723':(88558577,'A','C')};targets={'ast':'ENSG00000250156','opc':'ENSG00000247828'}
fields=['context','rsid','variant_id','target_gene_symbol','target_ensembl_id','pval','FDR','beta_non_ref','beta_se_non_ref','non_ref_af','qtl_dist_to_target','target_info','user_input']
out=[];slots={}
for ctx,gene in targets.items():
 tab=pysam.TabixFile(str(base/(ctx+'.bed.gz')));header=tab.header[0].lstrip('#').split('\t')
 for rs,(pos,ref,alt) in variants.items():
  matched=[]
  for line in tab.fetch('chr5',pos-1,pos):
   r=dict(zip(header,line.split('\t')))
   if r['variant_id']==f'chr5:{pos}:{ref}:{alt}' and r['target_ensembl_id']==gene:matched.append(r)
  slots[ctx+'_'+rs]={'row_count':len(matched),'status':'matched' if len(matched)==1 else 'absent' if not matched else 'duplicate'}
  for r in matched:
   assert 0<=float(r['pval'])<=1 and float(r['beta_se_non_ref'])>0
   out.append({'context':ctx,'rsid':rs,**{k:r[k] for k in fields if k not in ('context','rsid')}})
 tab.close()
with (p/'exact_qtl_rows.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
url='https://xqtl.niagads.org/download/rsid/rs62369151';content=subprocess.run(['curl','-fsSL','--max-time','30',url],capture_output=True,check=True).stdout
(p/'portal_response').write_bytes(content)
receipt={'dataset':'NIAGADS NG00184.v1 ROSMAP_CUIMC1_2_MIT nominal snuc-eQTL','source_hashes':hashes,'slots':slots,'padded_slot_family':4,'bonferroni_threshold':.0125,'portal':{'url':url,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'no_records_message':b'No records found for' in content},'rows_sha256':digest(p/'exact_qtl_rows.tsv','sha256'),'limitations':'Selected variants, genes and cells; four-slot correction does not cover prior searches. Different cell contexts share aging-brain donors. No causal direction, sleep–ADHD colocalization, or assay-to-gene mediation established.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
