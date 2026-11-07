#!/usr/bin/env python3
import csv,hashlib,json,math,urllib.request
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent
panel=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
with panel.open() as f:eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
assert len(eur)==503
variants=[('depression_hippocampus',117897197,'C','T'),('depression_cingulate',117907302,'C','T'),
 ('rna_hippocampus_lead',117905108,'T','G'),('rna_cingulate_lead',117908698,'A','C'),
 ('insomnia_lead_rs1456193',117918850,'T','C'),('alphagenome_rs62264764',117920728,'G','A')]
source='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(source);calls={};meta={}
for label,pos,ref,alt in variants:
 u=f'https://rest.ensembl.org/map/human/GRCh38/3:{pos}..{pos}/GRCh37?content-type=application/json'
 with urllib.request.urlopen(u,timeout=45) as f:mapping=json.load(f)
 mapped=mapping['mappings'];assert len(mapped)==1
 m=mapped[0]['mapped'];assert m['seq_region_name']=='3' and m['strand']==1 and m['start']==m['end'];pos37=m['start']
 records=[r for r in vcf.fetch('3',pos37-1,pos37) if r.pos==pos37 and r.ref==ref and alt in r.alts]
 item={'label':label,'pos38':pos,'ref':ref,'alt':alt,'mapping_url':u,'mapping':mapping,'pos37':pos37,'matching_records':len(records)}
 if len(records)!=1:meta[label]=item;continue
 r=records[0];ai=r.alleles.index(alt);valid={};excluded=0
 for s in sorted(eur & set(vcf.header.samples)):
  call=r.samples[s];g=call.get('GT')
  if not g or len(g)!=2 or not call.phased or any(x is None or x not in (0,ai) for x in g):excluded+=1;continue
  valid[s]=tuple(int(x==ai) for x in g)
 calls[label]=valid;item.update(rsid=r.id,vcf_alleles=list(r.alleles),eligible=len(valid),excluded=excluded);meta[label]=item
 print(label,pos37,'eligible',len(valid),flush=True)
rows=[]
for i,(a,*_) in enumerate(variants):
 for b,*_ in variants[i+1:]:
  if a not in calls or b not in calls:continue
  shared=sorted(set(calls[a])&set(calls[b]));x=[];y=[]
  for s in shared:x.extend(calls[a][s]);y.extend(calls[b][s])
  n=len(x);fa=sum(x)/n;fb=sum(y)/n;fab=sum(t*u for t,u in zip(x,y))/n;den=math.sqrt(fa*(1-fa)*fb*(1-fb));r=(fab-fa*fb)/den if den else None
  if r is not None:
   assert abs(r)<=1+1e-12
   r=max(-1.0,min(1.0,r))
  rows.append({'variant1':a,'variant2':b,'n_individuals':len(shared),'n_haplotypes':n,'alt1_frequency':fa,'alt2_frequency':fb,'r':r,'r2':r*r if r is not None else None})
with (p/'ld.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['variant1','variant2','n_individuals','n_haplotypes','alt1_frequency','alt2_frequency','r','r2'],delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'variants':meta,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'eur_panel_n':len(eur),'vcf_url':source,
 'ld_table_sha256':hashlib.sha256((p/'ld.tsv').read_bytes()).hexdigest(),'pairs':len(rows),
 'note':'Post-selected, phased EUR reference; remote VCF and index not whole-file hashed; not cohort LD or causal evidence.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
