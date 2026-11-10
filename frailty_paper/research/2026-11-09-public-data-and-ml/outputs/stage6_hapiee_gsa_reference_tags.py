"""Exploratory 1000G-EUR tagging of candidate loci by standard GSA probes."""
import csv,io,json,zipfile,hashlib
from pathlib import Path
import numpy as np,pysam
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';targets={}
for path,chrcol,poscol in [(o/'unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),(o/'unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]:
 for r in csv.DictReader(path.open()):
  key=r['variant'];x=(str(r[chrcol]),int(r[poscol]),frozenset([r['counted_allele'],r['other_allele']]))
  if key in targets and targets[key]!=x:raise ValueError('candidate contrast inconsistency')
  targets[key]=x
array={str(i):set() for i in range(1,23)}
with zipfile.ZipFile(w/'gsa_v3_a1_hg19.zip') as z:
 f=io.TextIOWrapper(z.open(z.namelist()[0]),encoding='utf-8-sig',newline='')
 for line in f:
  if line.strip()=='[Assay]':break
 for r in csv.DictReader(f):
  if r['IlmnID']=='[Controls]':break
  try:p=int(r['MapInfo'])
  except (TypeError,ValueError):continue
  if r['Chr'] in array and len(r['SNP'])==5 and r['SNP'][0]=='[' and r['SNP'][2]=='/' and r['SNP'][4]==']' and set(r['SNP'][1:4:2])<=set('ACGT') and r['Intensity_Only']=='0':array[r['Chr']].add(p)
with (w/'1kg_phase3_panel.tsv').open() as f:ids={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
def dosage(v,samples):
 a=[]
 for s in samples:
  g=v.samples[s].get('GT')
  if g is None or len(g)!=2 or None in g:return None
  a.append(sum(x==1 for x in g))
 a=np.array(a,dtype=float);return a if a.std()>0 else None
out=[]
for chrom in sorted({v[0] for v in targets.values()},key=int):
 f=pysam.VariantFile(url.format(chrom=chrom));samples=sorted(ids & set(f.header.samples))
 for name,(ch,pos,alleles) in targets.items():
  if ch!=chrom:continue
  a=[v for v in f.fetch(ch,pos-1,pos) if v.pos==pos and v.id and name in v.id.split(';') and len(v.alleles)==2 and set(v.alleles)==set(alleles)]
  row={'variant':name,'chr':ch,'pos37':pos,'reference_n':len(samples),'status':'candidate_missing_or_nonunique','tags_evaluated':0,'best_tag_id':'','best_tag_pos37':'','best_tag_r2':''}
  if len(a)==1:
   g=dosage(a[0],samples)
   if g is None:row['status']='candidate_incomplete_or_invariant'
   else:
    row['status']='no_eligible_reference_tag';best=(-1,None)
    for v in f.fetch(ch,max(0,pos-100001),pos+100000):
     if v.pos not in array[ch] or len(v.alleles)!=2 or v.pos==pos:continue
     x=dosage(v,samples)
     if x is None:continue
     r=float(np.corrcoef(g,x)[0,1]);r2=r*r;row['tags_evaluated']+=1
     if r2>best[0]:best=(r2,v)
    if best[1] is not None:
     r2,v=best;row.update(status='found',best_tag_id=v.id or '.',best_tag_pos37=v.pos,best_tag_r2=r2)
  out.append(row);print(name,row['status'],row['tags_evaluated'],row['best_tag_r2'],flush=True)
with (o/'stage6_hapiee_gsa_reference_tags.csv').open('w',newline='') as f:wr=csv.DictWriter(f,fieldnames=out[0]);wr.writeheader();wr.writerows(out)
(o/'stage6_hapiee_gsa_reference_tags_receipt.json').write_text(json.dumps({'source_vcf':url,'reference':'1000 Genomes phase 3 EUR (503 people)','manifest_sha256':hashlib.sha256((w/'gsa_v3_a1_hg19.zip').read_bytes()).hexdigest(),'window_bp':100000,'tag_eligibility':'manifest non-intensity-only SNV position, not allele/strand validated','candidate_eligibility':'unique rsid, position, exact biallelic pair, all complete noninvariant genotype calls','limits':'in-sample tag selection, same reference scoring, no cohort genotypes or imputation accuracy, no independent validation'},indent=2)+'\n')
