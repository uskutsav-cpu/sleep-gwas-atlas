#!/usr/bin/env python3
import collections,csv,hashlib,json,re,zipfile
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent;parent=p.parent;file=p/'table_s1.xlsx'
with zipfile.ZipFile(file) as z:assert z.testzip() is None
with (parent/'provisional_allelic_42.tsv').open() as f: candidates=list(csv.DictReader(f,delimiter='\t'))
assert len(candidates)==42 and len({x['rsid'] for x in candidates})==42
w=openpyxl.load_workbook(file,read_only=True,data_only=True)
assert w.sheetnames==['Description','Table S1A','Table S1B','Table S1C']
descriptions=list(w['Description'].iter_rows(values_only=True));assert 'MPRA summary statistics for enhancer' in descriptions[0][0]
tables={};counts={}
for sheet in ('Table S1A','Table S1B'):
 it=w[sheet].iter_rows(values_only=True); header=next(it);d=collections.defaultdict(list);count=0
 for row in it:
  if not any(x is not None for x in row):continue
  record=dict(zip(header,row)); count+=1;d[str(record['RSID'])].append(record)
 tables[sheet]=d;counts[sheet]=count
fields=['rsid','region','source_key','grch38','element_status','element_class','element_fdr','allelic_status','allelic_class','allelic_fdr','published_mpra_logfc_unoriented','v124_fdr','v134_fdr','listed_excluded']
out=[];summary=collections.Counter()
for c in candidates:
 m=re.fullmatch(r'(chr\d+):(\d+)_([ACGT]+):([ACGT]+)',c['source_key']);m38=re.fullmatch(r'(chr\d+):(\d+):([ACGT]+)>([ACGT]+)',c['grch38']);assert m and m38
 chrom,pos,ref,alt=m.groups();chr38,pos38,ref38,alt38=m38.groups();assert chrom==chr38 and {ref,alt}=={ref38,alt38}
 result={k:c[k] for k in ('rsid','region','source_key','grch38','listed_excluded')}
 for sheet,tag in [('Table S1A','element'),('Table S1B','allelic')]:
  records=tables[sheet].get(c['rsid'],[])
  hits=[r for r in records if r['CHR']==chrom and str(r['BP_HG19'])==pos and {str(r['A1']).upper(),str(r['A2']).upper()}=={ref,alt} and (sheet!='Table S1B' or str(r['BP_HG38'])==pos38)]
  status='unique_match' if len(hits)==1 else 'absent_rsid' if not records else 'ambiguous_or_mismatch'
  result[tag+'_status']=status;summary[tag+'_'+status]+=1
  if len(hits)==1:
   row=hits[0];result[tag+'_class']=row['Element_Class' if tag=='element' else 'Variant_Class'];result[tag+'_fdr']=row['MPRA_FDR']
   if tag=='allelic':result['published_mpra_logfc_unoriented']=row['MPRA_logFC']
 result['v124_fdr']=c['v124_adj.P.Val'];result['v134_fdr']=c['v134_adj.P.Val'];out.append(result)
with (p/'candidate_comparison.tsv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(out)
for r in out:
 if r.get('allelic_status')=='unique_match':
  pub=float(r['allelic_fdr'])<.05;summary['published_fdr_lt_005']+=pub
  for version in ('v124','v134'):
   val=r[version+'_fdr'];assert val not in ('','NA');pred=float(val)<.05
   summary[version+'_published_fdr_disagreements']+=pub!=pred
  summary['variant_class:'+str(r['allelic_class'])]+=1
 if r.get('element_status')=='unique_match':summary['element_class:'+str(r['element_class'])]+=1
 if r.get('element_status')=='unique_match' and r.get('allelic_status')=='unique_match':summary['joint_classes:'+str(r['element_class'])+'/'+str(r['allelic_class'])]+=1
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'publisher_url':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867424014351-mmc2.xlsx','publisher_article':'https://doi.org/10.1016/j.cell.2024.12.022','xlsx_zip_integrity':'passed','xlsx_bytes':file.stat().st_size,'xlsx_sha256':sha(file),'candidate_source_sha256':sha(parent/'provisional_allelic_42.tsv'),'data_rows_by_sheet':counts,'candidate_count':len(out),'summary':dict(summary),'candidate_comparison_sha256':sha(p/'candidate_comparison.tsv'),'limitations':'Post-selected candidate comparison; matching and FDR arithmetic do not validate genomic target genes, genetic replication, MPRA model interpretation or native chromatin effects. MPRA sign orientation is unverified.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
