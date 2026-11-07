import csv,gzip,json,hashlib,math,collections
from pathlib import Path
p=Path(__file__).resolve().parent;src=p.parents[1]/'chr5_developmental_caqtl/retrieval_receipt.json';inputs=json.loads(src.read_text());out=[];receipt=[]
for item in inputs['members']:
 path=Path(item['file']);sha=hashlib.sha256(path.read_bytes()).hexdigest();assert sha==item['inner_gzip_sha256'];n=0;matched=0
 with gzip.open(path,'rt') as f:
  assert f.readline().split()==['snp','beta_coef','pvalue','peak','EffectAllele']
  for line in f:
   n+=1
   if not line.startswith('chr5:104568525:'):continue
   snp,beta,pv,peak,allele=line.split();parts=snp.split(':');matched+=1
   if frozenset(parts[2:])!=frozenset(['A','G']):continue
   b=float(beta);q=float(pv);assert allele in ['A','G'] and 0<=q<=1 and math.isfinite(b)
   _,a,z=peak.split('_');a=int(a);z=int(z)
   out.append({'cell':item['cell'],'variant':snp,'beta_alt_G':b if allele=='G' else -b,'p':q,'peak':peak,'source_effect_allele':allele,'variant_inside_peak_assumed_bed':a<=104568524<z})
 receipt.append({'cell':item['cell'],'rows_read':n,'position_matches':matched,'compatible':sum(x['cell']==item['cell'] for x in out),'gzip_sha256':sha,'gzip_full_read':True})
for x in out:x['bonferroni_across_returned_rows']=min(1,x['p']*len(out))
fields=['cell','variant','beta_alt_G','p','bonferroni_across_returned_rows','peak','source_effect_allele','variant_inside_peak_assumed_bed']
with (p/'matches.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'datasets':receipt,'matched_rows':len(out),'correction':'post-selected returned rows only; not broader project family','source_receipt_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256((p/'matches.tsv').read_bytes()).hexdigest()},indent=2)+'\n');print(json.dumps(receipt,indent=2));print('rows',len(out),'best',sorted(out,key=lambda x:x['p'])[:3])
