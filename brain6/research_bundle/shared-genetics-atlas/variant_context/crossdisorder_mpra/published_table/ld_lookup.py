#!/usr/bin/env python3
import csv,hashlib,json,math,re
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent
panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
with panel.open() as f:eur={x['sample'] for x in csv.DictReader(f,delimiter='\t') if x['super_pop']=='EUR'}
assert len(eur)==503
with (p.parent/'ranked_coverage.tsv').open() as f: rows={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
ids=('rs325485','rs6874138','rs77960','rs2431108')
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);calls={};metadata={}
for rs in ids:
 match=re.fullmatch(r'chr5:(\d+)_([ACGT]+):([ACGT]+)',rows[rs]['source_key']);assert match
 pos,ref,alt=match.groups();pos=int(pos)
 found=[x for x in vcf.fetch('5',pos-1,pos) if x.pos==pos and x.id==rs and ref in x.alleles and alt in x.alleles]
 assert len(found)==1
 rec=found[0];ir=rec.alleles.index(ref);ia=rec.alleles.index(alt);d={}
 for sample in sorted(eur & set(vcf.header.samples)):
  x=rec.samples[sample];gt=x.get('GT')
  if gt and len(gt)==2 and x.phased and all(y is not None and y in (ir,ia) for y in gt):d[sample]=tuple(int(y==ia) for y in gt)
 calls[rs]=d;metadata[rs]={'pos37':pos,'reference_alleles':list(rec.alleles),'n_phased':len(d)}
out=[]
for a in ids[:2]:
 for b in ids[2:]:
  common=sorted(set(calls[a])&set(calls[b])); x=[];y=[]
  for sample in common:x.extend(calls[a][sample]);y.extend(calls[b][sample])
  fa=sum(x)/len(x);fb=sum(y)/len(y);fab=sum(u*v for u,v in zip(x,y))/len(x);den=math.sqrt(fa*(1-fa)*fb*(1-fb));assert den>0
  corr=(fab-fa*fb)/den;out.append({'reporter_variant':a,'trait_lead':b,'n_individuals':len(common),'r':corr,'r2':corr*corr})
with (p/'pairwise_ld.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['reporter_variant','trait_lead','n_individuals','r','r2'],delimiter='\t');w.writeheader();w.writerows(out)
receipt={'reference_url':url,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'records':metadata,'pairwise':out,'limitations':'Post-selected comparisons, 503 European-labeled samples, phased haplotype correlation; remote indexed VCF not whole-file authenticated, reference not matched to GWAS or MPRA donor.'}
(p/'ld_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
