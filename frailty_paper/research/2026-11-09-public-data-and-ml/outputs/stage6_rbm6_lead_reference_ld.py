#!/usr/bin/env python3
import json,csv,hashlib,itertools
from pathlib import Path
import numpy as np,pandas as pd,pysam
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs';summary=json.loads((out/'stage6_rbm6_five_gwas_region_summary.json').read_text());leads={}
for c in summary['comparisons'].values():
 for k in ['sleep_best_shared_variant','frailty_best_shared_variant']:
  r=c[k];leads[r['rsid']]=r
leads['rs2353586']={'rsid':'rs2353586','pos':50008613,'ref':'G','alt':'T'}
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';index=work/'stage6_1000g_chr3_phase3.tbi';t=pysam.TabixFile(url,index=str(index));samples=[x for x in t.header if x.startswith('#CHROM')][0].split('\t')[9:];panel={x['sample']:x for x in csv.DictReader((work/'1kg_phase3_panel.tsv').open(),delimiter='\t')};groups={name:[i for i,s in enumerate(samples) if s in panel and pred(panel[s])] for name,pred in [('EUR',lambda x:x['super_pop']=='EUR'),('FIN',lambda x:x['pop']=='FIN')]};assert len(groups['EUR'])==503 and len(groups['FIN'])==99;lift=LiftOver(str(work/'hg38ToHg19.over.chain.gz'));vec={};inventory=[]
for rs,r in sorted(leads.items()):
 m=lift.convert_coordinate('chr3',r['pos']-1);status='';match=[]
 if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':status='nonunique_or_nonpositive_liftover'
 else:
  pos=m[0][1]+1
  for l in t.fetch('3',pos-1,pos):
   z=l.split('\t')
   if int(z[1])==pos and z[3]==r['ref'] and r['alt'] in z[4].split(','):match.append(z)
  if len(match)!=1:status='no_unique_exact_reference_allele'
  else:
   z=match[0];idx=z[4].split(',').index(r['alt'])+1;gi=z[8].split(':').index('GT');d=[]
   for s in z[9:]:
    gt=s.split(':')[gi].replace('|','/').split('/');d.append(sum(int(a)==idx for a in gt) if len(gt)==2 and all(a.isdigit() for a in gt) else float('nan'))
   vec[rs]=np.array(d,dtype=float);status='matched'
 inventory.append({'rsid':rs,'pos_grch38':r['pos'],'ref':r['ref'],'alt':r['alt'],'status':status,'pos_grch37':pos if len(m)==1 and m[0][0]=='chr3' else ''})
rows=[]
for a,b in itertools.combinations(sorted(vec),2):
 for group,idx in groups.items():
  x=vec[a][idx];y=vec[b][idx];valid=np.isfinite(x)&np.isfinite(y);r=float(np.corrcoef(x[valid],y[valid])[0,1]);rows.append({'snp1':a,'snp2':b,'reference_group':group,'n':int(sum(valid)),'r':r,'r2':r*r})
pd.DataFrame(inventory).to_csv(out/'stage6_rbm6_lead_reference_coverage.csv',index=False);pd.DataFrame(rows).to_csv(out/'stage6_rbm6_lead_reference_ld.csv',index=False);(out/'stage6_rbm6_lead_reference_receipt.json').write_text(json.dumps({'source':url,'index_sha256':hashlib.sha256(index.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256((work/'1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'reference_sample_counts':{k:len(v) for k,v in groups.items()},'matched':len(vec),'candidates':len(leads)},indent=2)+'\n');print('matched',len(vec),'of',len(leads));print(pd.DataFrame(rows)[(pd.DataFrame(rows).snp1=='rs10865954')&(pd.DataFrame(rows).snp2=='rs2282749')].to_string(index=False))
