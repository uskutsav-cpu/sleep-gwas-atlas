import csv,json,hashlib
import numpy as np,pysam
from pathlib import Path
root=Path(__file__).resolve().parents[1]
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}
u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
t=pysam.TabixFile(u,index=str(root/'work/stage5_1000g_chr5_phase3.tbi'));samples=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:]; rows=[]
for pos,pair in [(135538341,{'A','C'}),(135467309,{'C','G'})]:
 hits=[r.split('\t') for r in t.fetch('5',pos-1,pos) if set((r.split('\t')[3],r.split('\t')[4]))==pair];assert len(hits)==1;rows.append(hits[0])
out={'source':u,'records':[{'grch37_position':int(r[1]),'id':r[2],'ref':r[3],'alt':r[4]} for r in rows],'groups':{},'caveat':'External reference LD is not target-study LD. The top SNP is palindromic; do not infer target allele effect direction.'}
for pop in ['EUR','EAS']:
 idx=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']==pop];vec=[]
 for row in rows:
  x=[]
  for i in idx:
   gt=row[9+i].split(':')[0].replace('|','/').split('/');x.append(sum(int(a)==1 for a in gt) if len(gt)==2 and all(a!='.' for a in gt) else np.nan)
  vec.append(np.array(x))
 mask=np.isfinite(vec[0])&np.isfinite(vec[1]);r=float(np.corrcoef(vec[0][mask],vec[1][mask])[0,1]);out['groups'][pop]={'n':int(mask.sum()),'r_squared':r*r}
(root/'outputs/stage6_tenk10k_smad5_reference_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
