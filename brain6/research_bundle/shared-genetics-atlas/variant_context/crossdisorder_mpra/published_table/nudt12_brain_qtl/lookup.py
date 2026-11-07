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
variants={'rs6874138':(104563895,'T','C'),'rs325485':(104659667,'A','G')};target='ENSG00000112874'
fields=['context','rsid','variant_id','target_gene_symbol','target_ensembl_id','pval','FDR','beta_non_ref','beta_se_non_ref','non_ref_af','qtl_dist_to_target','target_info','user_input']
results=[];slots={}
for ctx in ('ast','opc'):
 tab=pysam.TabixFile(str(base/(ctx+'.bed.gz')));header=tab.header[0].lstrip('#').split('\t')
 for rs,(pos,ref,alt) in variants.items():
  matches=[]
  for line in tab.fetch('chr5',pos-1,pos):
   r=dict(zip(header,line.split('\t')))
   if r['variant_id']==f'chr5:{pos}:{ref}:{alt}' and r['target_ensembl_id']==target:matches.append(r)
  slots[ctx+'_'+rs]={'row_count':len(matches),'status':'matched' if len(matches)==1 else 'absent' if not matches else 'duplicate'}
  for r in matches:
   pv=float(r['pval']);assert 0<=pv<=1 and math.isfinite(pv) and float(r['beta_se_non_ref'])>0
   results.append({'context':ctx,'rsid':rs,**{key:r[key] for key in fields if key not in ('context','rsid')}})
 tab.close()
with (p/'exact_qtl_rows.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(results)
portal={}
for rs in variants:
 url='https://xqtl.niagads.org/download/rsid/'+rs
 content=subprocess.run(['curl','-fsSL','--max-time','30',url],check=True,capture_output=True).stdout;ctype='not recorded'
 (p/(rs+'_portal_response')).write_bytes(content)
 portal[rs]={'url':url,'content_type':ctype,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'no_records_message':('No records found for' in content.decode('utf8','replace'))}
receipt={'dataset':'NIAGADS NG00184.v1, ROSMAP_CUIMC1_2_MIT nominal snuc-eQTL astrocytes and oligodendrocyte progenitors','source_hashes':hashes,'slots':slots,'slot_family_count':4,'bonferroni_threshold':.0125,'portal':portal,'rows_sha256':digest(p/'exact_qtl_rows.tsv','sha256'),'limitations':'Post-selected two-variant two-cell descriptive family, not correction for all genomic searches; absent and filtered records are not negative tests. The same aging-brain donors may contribute across cell types; QTL is not disease colocalization or gene mediation.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
