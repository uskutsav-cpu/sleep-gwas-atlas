import csv,json,hashlib
from pathlib import Path
import pysam,numpy as np
from pyliftover import LiftOver
variants={'eqtl_lead':('rs6738144',66425823,'C','A'),'insomnia_lead':('rs11693221',66572854,'C','T'),'score_marker':('rs113851554',66523432,'G','T'),'hfrs_lead':('rs111736174',66413206,'C','A'),'fi_lead':('rs2194667',66847522,'A','G')}
lo=LiftOver('work/hg38ToHg19.over.chain.gz');url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';f=pysam.VariantFile(url)
with open('work/1kg_phase3_panel.tsv') as h:ids={r['sample'] for r in csv.DictReader(h,delimiter='\t') if r['super_pop']=='EUR'}
samples=sorted(ids & set(f.header.samples));data={};status={}
for tag,(rs,pos,ref,alt) in variants.items():
 m=lo.convert_coordinate('chr2',pos-1);status[tag]={'rsid':rs,'GRCh38_pos':pos,'alleles':[ref,alt],'mapping':m}
 if len(m)!=1 or m[0][0]!='chr2' or m[0][2]!='+':continue
 pos37=m[0][1]+1;xs=[x for x in f.fetch('2',pos37-1,pos37) if x.pos==pos37 and x.id and rs in x.id.split(';')];status[tag]['vcf_count']=len(xs)
 if len(xs)!=1:continue
 x=xs[0];status[tag]['vcf_alleles']=x.alleles
 if set(x.alleles)!={ref,alt}:continue
 ix=x.alleles.index(alt);vals={}
 for s in samples:
  c=x.samples[s];g=c.get('GT')
  if g and len(g)==2 and None not in g and c.phased:vals[s]=[int(z==ix) for z in g]
 data[tag]=vals
out=[]
for t in ['insomnia_lead','score_marker','hfrs_lead','fi_lead']:
 row={'variants':['eqtl_lead',t],'status':'unavailable'}
 if 'eqtl_lead' in data and t in data:
  ids=sorted(data['eqtl_lead'].keys() & data[t].keys());x=np.array([z for s in ids for z in data['eqtl_lead'][s]]);y=np.array([z for s in ids for z in data[t][s]])
  if x.std()>0 and y.std()>0:r=float(np.corrcoef(x,y)[0,1]);row.update(status='available',n_individuals=len(ids),r=r,r2=r*r)
 out.append(row)
res={'reference':url,'population':'1000 Genomes phase 3 EUR','variants':status,'pairs':out,'caveat':'postselected marginal leads; reference LD only; no conditional association, colocalization or clinical prediction'}
Path('outputs/stage6_meis1_lead_reference_ld.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
