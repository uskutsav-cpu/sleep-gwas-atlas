import csv,pysam,numpy as np,hashlib,json,pathlib
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
wanted=[('rs2838788',46542374,'T','C'),('rs8129315',46498863,'G','A'),('rs35424628',46701438,'G','A')]
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'}
f=pysam.VariantFile(url); data={};results=[]
for rs,pos,a,b in wanted:
 hits=[v for v in f.fetch('21',pos-1,pos) if v.id==rs and v.pos==pos and len(v.alleles)==2 and set(v.alleles)=={a,b}]
 if len(hits)!=1:
  results.append({'variant':rs,'status':'exact biallelic match count '+str(len(hits))});continue
 v=hits[0];ix=v.alleles.index(a);calls={}
 for s in sorted(eur & set(f.header.samples)):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(z is not None for z in gt) and v.samples[s].phased:calls[s]=[int(z==ix) for z in gt]
 data[rs]=calls;results.append({'variant':rs,'alleles':list(v.alleles),'counted':a,'n_called':len(calls),'status':'matched'})
pairs=[]
for rs in ['rs8129315','rs35424628']:
 if rs not in data or 'rs2838788' not in data: pairs.append({'variant':rs,'status':'not calculated'});continue
 samples=sorted(data[rs].keys() & data['rs2838788'].keys());x=np.array([z for s in samples for z in data['rs2838788'][s]]);y=np.array([z for s in samples for z in data[rs][s]])
 assert x.std()>0 and y.std()>0
 r=float(np.corrcoef(x,y)[0,1]);pairs.append({'variant':rs,'n_individuals':len(samples),'n_haplotypes':len(x),'r':r,'r2':r*r,'status':'matched phased EUR'})
out={'source':url,'build':'GRCh37 1000 Genomes reference; Mak build assumed','reference':'503 EUR-labeled phase 3 samples, phased haplotypes','variants':results,'pairs':pairs,'limitations':'No proof of same causal variants; reference-panel correlation differs from conditional GWAS; source alleles are Mak A1 and candidate T; no individual data saved.'}
p=pathlib.Path('outputs/adarb1_mak_regional_ld.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
