#!/usr/bin/env python3
import csv,hashlib,json,math,re
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
with panel.open() as f:eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
assert len(eur)==503
with (p.parent/'ranked_coverage.tsv').open() as f:variants={r['rsid']:r for r in csv.DictReader(f,delimiter='\t')}
ids=('rs62369151','rs4916723');url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);calls={};records={}
for rs in ids:
 m=re.fullmatch(r'chr5:(\d+)_([ACGT]+):([ACGT]+)',variants[rs]['source_key']);assert m
 pos,ref,alt=m.groups();pos=int(pos);hits=[r for r in vcf.fetch('5',pos-1,pos) if r.pos==pos and r.id==rs and ref in r.alleles and alt in r.alleles];assert len(hits)==1
 rec=hits[0];ir=rec.alleles.index(ref);ia=rec.alleles.index(alt);data={}
 for sample in sorted(eur & set(vcf.header.samples)):
  s=rec.samples[sample];gt=s.get('GT')
  if gt and s.phased and len(gt)==2 and all(g is not None and g in (ir,ia) for g in gt):data[sample]=tuple(int(g==ia) for g in gt)
 calls[rs]=data;records[rs]={'position37':pos,'ref':ref,'alt':alt,'n_phased':len(data)}
common=sorted(set(calls[ids[0]])&set(calls[ids[1]]));a=[];b=[]
for sample in common:a.extend(calls[ids[0]][sample]);b.extend(calls[ids[1]][sample])
fa=sum(a)/len(a);fb=sum(b)/len(b);f11=sum(x*y for x,y in zip(a,b))/len(a);den=math.sqrt(fa*(1-fa)*fb*(1-fb));assert den>0
r=(f11-fa*fb)/den
out={'reference_url':url,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'variants':records,'n_individuals':len(common),'haplotype_r':r,'r_squared':r*r,'limitations':'Post-hoc selected pair in 503 EUR labeled individuals; indexed remote reference not fully authenticated, and LD is not from the ADHD GWAS or assay donor.'}
(p/'rs62369151_ld_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
