import csv,json
import numpy as np
import pysam
from pathlib import Path
root=Path(__file__).resolve().parents[1]
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}
u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
t=pysam.TabixFile(u,index=str(root/'work/stage6_1000g_chr7_phase3.tbi'))
samples=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:]
idx=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR']; assert len(idx)==503
records=[];vec=[]
for pos,rs,count in [(2074122,'rs28528997','A'),(2052318,'rs62444881','T')]:
 rows=[x.split('\t') for x in t.fetch('7',pos-1,pos) if x.split('\t')[2]==rs]; assert len(rows)==1
 row=rows[0]; alleles=[row[3]]+row[4].split(','); assert count in alleles and len(alleles)==2
 ai=alleles.index(count);dos=[]
 for i in idx:
  gt=row[9+i].split(':')[0].replace('|','/').split('/')
  dos.append(sum(int(a)==ai for a in gt) if len(gt)==2 and all(a!='.' for a in gt) else np.nan)
 a=np.array(dos,dtype=float);vec.append(a);records.append({'rsid':rs,'position_grch37':pos,'reference_alleles':alleles,'counted_reference_allele':count,'counted_frequency':float(np.nanmean(a)/2),'missing':int(np.isnan(a).sum())})
mask=np.isfinite(vec[0])&np.isfinite(vec[1]); r=float(np.corrcoef(vec[0][mask],vec[1][mask])[0,1])
out={'source':u,'population':'1000 Genomes phase 3 EUR','n_complete':int(mask.sum()),'records':records,'signed_r_candidate_A_vs_lead_T':r,'r_squared':r*r,'caveat':'T is the complement of published HRS effect allele A, but the original published table lacks the other allele; strand alignment and direction remain unconfirmed. Association is from selected published hits, and reference LD does not validate the full score.'}
(root/'outputs/stage6_mad1l1_hrs_proxy_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
