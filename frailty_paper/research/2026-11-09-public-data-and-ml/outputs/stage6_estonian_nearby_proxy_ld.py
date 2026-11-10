"""EUR 1000G phase 3 LD for three nearby candidate/Estonian lead pairs."""
import csv,json
from pathlib import Path
import numpy as np,pysam
root=Path(__file__).resolve().parents[1]
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}
pairs=[('3',('rs817509',117750861),('rs1093463',117744762)),('7',('rs17312861',114125510),('rs1989903',114137940)),('11',('rs2155292',112910783),('rs745631',112945892))]
output=[]
for chrom,*sites in pairs:
 url=f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
 t=pysam.TabixFile(url,index=str(root/f'work/stage6_1000g_chr{chrom}_phase3.tbi'))
 samples=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:]; idx=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR'];assert len(idx)==503
 vectors=[]; info=[]
 for rs,pos in sites:
  hits=[x.split('\t') for x in t.fetch(chrom,pos-1,pos) if x.split('\t')[2]==rs]
  if len(hits)!=1 or ',' in hits[0][4]:
   info.append({'rsid':rs,'position_GRCh37':pos,'problem':'not exactly one biallelic record'});continue
  row=hits[0]; vec=[]
  for i in idx:
   gt=row[9+i].split(':')[0].replace('|','/').split('/')
   vec.append(sum(int(a)==1 for a in gt) if len(gt)==2 and all(a!='.' for a in gt) else np.nan)
  vectors.append(np.array(vec));info.append({'rsid':rs,'position_GRCh37':pos,'ref':row[3],'alt':row[4]})
 result={'sites':info,'reference':url}
 if len(vectors)==2:
  mask=np.isfinite(vectors[0])&np.isfinite(vectors[1]); r=float(np.corrcoef(vectors[0][mask],vectors[1][mask])[0,1]);result.update(n=int(mask.sum()),alt_dosage_r=r,r_squared=r*r)
 output.append(result)
(root/'outputs/stage6_estonian_nearby_proxy_ld.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
