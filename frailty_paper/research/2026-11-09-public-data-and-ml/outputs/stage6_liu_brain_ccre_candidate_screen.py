#!/usr/bin/env python3
"""Overlap 22 candidate positions with public Liu et al. processed peak BEDs.
Assembly and exact BED start semantics remain conditional; boundary hits excluded.
"""
import csv,hashlib,json,re,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote
p=Path(__file__).resolve().parent
base='https://compbio.mit.edu/AD_Multiomic_MultiRegion/03_cCREs/'
def get(url):return subprocess.run(['curl','-fsSL','--retry','2','--max-time','60',url],stdout=subprocess.PIPE,check=True).stdout
def read(n):
 with (p/n).open(newline='') as f:return list(csv.DictReader(f))
class Parser(HTMLParser):
 def __init__(self):super().__init__();self.paths=[]
 def handle_starttag(self,tag,attrs):
  if tag=='a':
   for k,v in attrs:
    if k=='href' and v and re.fullmatch(r'[A-Za-z0-9_\-.]+\.bed',v):self.paths.append(v)
index=get(base);parser=Parser();parser.feed(index.decode());paths=sorted(set(parser.paths));assert len(paths)==68 and 'Union_PeakSet.bed' in paths
score=read('unvalidated_22_shared_direction_score.csv');inventory={r['variant']:r for r in read('stage6_tpmi_shared22_insomnia_inventory.csv')};atlas=read('stage5_functional_atlas_with_priority_splice.csv');variants=[]
for row in score:
 name=row['variant'];v=inventory[name]
 if v.get('pos38') and v['queue_identity_match']=='True':chrom='chr'+v['chr38'];pos=int(v['pos38']);ref=v['ref'];alt=v['alt']
 else:
  hits={x['variant_GRCh38_ref_alt'] for x in atlas if x['input_classes']=='exploratory_marker' and x['input_labels']==name};assert len(hits)==1
  m=re.fullmatch(r'(chr\d+):(\d+):([ACGT]+)>([ACGT]+)',next(iter(hits)));assert m;chrom,pos,ref,alt=m.groups();pos=int(pos)
 if {ref,alt}!={row['counted_allele'],row['other_allele']}:raise ValueError('allele mismatch '+name)
 variants.append({'variant':name,'chrom':chrom,'pos':pos,'ref':ref,'alt':alt})
assert len(variants)==22
bychr={}
for v in variants:bychr.setdefault(v['chrom'],[]).append(v)
def scan(file):
 data=get(base+quote(file));hits=[];lines=0;boundary=[]
 for raw in data.splitlines():
  lines+=1;parts=raw.split(b'\t');assert len(parts)>=3
  chrom=parts[0].decode();start=int(parts[1]);end=int(parts[2]);assert 0<=start<end
  for v in bychr.get(chrom,[]):
   pos=v['pos']
   if start<pos<end:hits.append((v['variant'],chrom,start,end,pos))
   elif pos==start or pos==end:boundary.append((v['variant'],chrom,start,end,pos))
 return {'file':file,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data),'rows':lines,'hits':hits,'boundary':boundary}
results=[]
with ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(scan,f):f for f in paths}
 for job in as_completed(jobs):
  results.append(job.result());
  if len(results)%10==0:print(f'scanned {len(results)}/{len(paths)} files',flush=True)
results.sort(key=lambda r:r['file']);matches=[]
for r in results:
 for name,chrom,start,end,pos in r['hits']:
  matches.append({'variant':name,'chromosome':chrom,'position_grch38':pos,'source_file':r['file'],'cell_subtype':r['file'][:-4],'peak_start':start,'peak_end':end,'status':'interior_overlap_conditional_on_source_build'})
fields=['variant','chromosome','position_grch38','source_file','cell_subtype','peak_start','peak_end','status']
with (p/'stage6_liu_brain_ccre_candidate_overlaps.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(matches)
receipt={'index_url':base,'index_sha256':hashlib.sha256(index).hexdigest(),'files':[{'file':r['file'],'sha256':r['sha256'],'bytes':r['bytes'],'rows':r['rows'],'interior_hits':len(r['hits']),'boundary_hits':len(r['boundary'])} for r in results],'candidate_count':len(variants),'cell_subtype_files':len(paths)-1,'matched_candidate_subtype_pairs':len([m for m in matches if m['source_file']!='Union_PeakSet.bed']),'matched_variants':sorted({m['variant'] for m in matches if m['source_file']!='Union_PeakSet.bed'}),'boundary_hit_count':sum(len(r['boundary']) for r in results),'source_version':'https://doi.org/10.1016/j.cell.2025.06.031','interpretation':'Exploratory position overlap only; assume source BED GRCh38 pending source-file-specific confirmation. Require start<one-based position<end to exclude both boundary conventions. Does not imply allele direction, cell-type enrichment, target gene, disease causality, or score prediction.'}
(p/'stage6_liu_brain_ccre_candidate_overlaps.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['candidate_count','cell_subtype_files','matched_candidate_subtype_pairs','matched_variants','boundary_hit_count']},indent=2))
