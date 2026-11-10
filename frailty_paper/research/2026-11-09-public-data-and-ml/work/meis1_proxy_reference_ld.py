import csv,json,math
from pathlib import Path
import pysam
import numpy as np
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
variants={'rs13033745':(66660403,'T','G'),'rs1519102':(66677816,'C','G'),'rs113851554':(66750564,'T','G'),'rs78206187':(67022234,'G','A')}
with open('work/1kg_phase3_panel.tsv') as f: ids={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
v=pysam.VariantFile(url); samples=sorted(ids & set(v.header.samples)); data={}; info={}
for rs,(pos,a,b) in variants.items():
 matches=[x for x in v.fetch('2',pos-1,pos) if x.pos==pos and x.id and rs in x.id.split(';')]
 info[rs]={'GRCh37_pos':pos,'records':len(matches),'required_contrast':[a,b]}
 if len(matches)!=1: continue
 x=matches[0];info[rs]['vcf_alleles']=x.alleles
 if set(x.alleles)!={a,b}:continue
 k=x.alleles.index(a);d={}
 for sample in samples:
  call=x.samples[sample];g=call.get('GT')
  if g and len(g)==2 and all(y is not None for y in g) and call.phased:d[sample]=[int(y==k) for y in g]
 data[rs]=d
pairs=[]
for a,b in [('rs13033745','rs1519102'),('rs13033745','rs113851554'),('rs1519102','rs113851554'),('rs13033745','rs78206187'),('rs1519102','rs78206187'),('rs113851554','rs78206187')]:
 row={'variants':[a,b],'status':'unavailable'}
 if a in data and b in data:
  keys=sorted(data[a].keys() & data[b].keys()); x=np.array([z for k in keys for z in data[a][k]]); y=np.array([z for k in keys for z in data[b][k]])
  if len(x)>0 and np.std(x)>0 and np.std(y)>0:
   r=float(np.corrcoef(x,y)[0,1]);row.update(status='available',n_individuals=len(keys),n_haplotypes=len(x),r=r,r2=r*r)
 pairs.append(row)
out={'source':url,'population':'1000 Genomes phase 3 EUR; phased haplotypes','variants':info,'pairs':pairs,'limitations':'LD in a reference is not target-study LD, causality, score validation, or proof of source multiallelic contrast.'}
Path('outputs/stage6_meis1_npc_proxy_reference_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
