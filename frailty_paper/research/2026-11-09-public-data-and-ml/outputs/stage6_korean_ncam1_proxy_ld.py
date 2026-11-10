from pathlib import Path
import csv,json,numpy as np,pysam
root=Path(__file__).resolve().parents[1]
u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
t=pysam.VariantFile(u)
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}
results={}
for group in ['EAS','EUR']:
 samples=[s for s in t.header.samples if s in panel and panel[s]['super_pop']==group];vec=[];records=[]
 for pos,rs,allele in [(112910783,'rs2155292','G'),(112885527,'rs10891490','C')]:
  rows=[x for x in t.fetch('11',pos-1,pos) if x.pos==pos and x.id==rs];assert len(rows)==1
  row=rows[0];assert len(row.alleles)==2 and allele in row.alleles
  ai=row.alleles.index(allele);v=[]
  for name in samples:
   gt=row.samples[name].get('GT');v.append(sum(x==ai for x in gt) if gt and len(gt)==2 and None not in gt else np.nan)
  vec.append(np.array(v));records.append({'rsid':rs,'position_grch37':pos,'reference_alleles':list(row.alleles),'counted':allele,'frequency':float(np.nanmean(v)/2)})
 mask=np.isfinite(vec[0])&np.isfinite(vec[1]);r=float(np.corrcoef(vec[0][mask],vec[1][mask])[0,1]);results[group]={'n_called':int(mask.sum()),'r':r,'r2':r*r,'records':records}
output={'reference':u,'results':results,'caveat':'Panel LD, not cohort LD; rs2155292 has other alleles in other resources. The Korean reported effect-allele/other-allele semantics, build and direction need source confirmation.'}
(root/'outputs/stage6_korean_ncam1_proxy_ld.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
