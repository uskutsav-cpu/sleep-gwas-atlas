#!/usr/bin/env python3
"""Build and check public 1000G EUR reference LD on the five-GWAS common support."""
from pathlib import Path
import json,csv,hashlib
import numpy as np,pandas as pd,pysam
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs';data=pd.read_csv(out/'stage6_rbm6_five_gwas_common_region.csv.gz');assert data.rsid.is_unique and len(data)==2256
lift=LiftOver(str(work/'hg38ToHg19.over.chain.gz'));lookup={};status={}
for r in data.itertuples():
 m=lift.convert_coordinate('chr3',int(r.pos)-1)
 if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':status[r.rsid]='nonunique_liftover';continue
 key=(int(m[0][1])+1,r.ref,r.alt)
 if key in lookup:raise ValueError('duplicate variant key')
 lookup[key]=r.rsid;status[r.rsid]='absent_reference'
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';index=work/'stage6_1000g_chr3_phase3.tbi';t=pysam.TabixFile(url,index=str(index)); samples=[s for s in t.header if s.startswith('#CHROM')][0].split('\t')[9:];panel={x['sample']:x for x in csv.DictReader((work/'1kg_phase3_panel.tsv').open(),delimiter='\t')};eur=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR'];assert len(eur)==503
vec={}; mi=min(k[0] for k in lookup);ma=max(k[0] for k in lookup)
for line in t.fetch('3',mi-1,ma):
 z=line.split('\t'); key=(int(z[1]),z[3],z[4]); rs=lookup.get(key)
 if not rs:continue
 if rs in vec:status[rs]='duplicate_vcf_records';vec.pop(rs,None);continue
 gt_idx=z[8].split(':').index('GT');vals=[]
 for i in eur:
  gt=z[9+i].split(':')[gt_idx].replace('|','/').split('/')
  vals.append(sum(int(a)==1 for a in gt) if len(gt)==2 and all(a in ('0','1') for a in gt) else np.nan)
 v=np.asarray(vals,dtype=np.float64);missing=float(np.mean(~np.isfinite(v)));af=float(np.nanmean(v)/2)
 if missing>.01:status[rs]='missing_genotypes_gt_1pct';continue
 if min(af,1-af)<.01:status[rs]='reference_maf_lt_1pct';continue
 v[~np.isfinite(v)]=np.nanmean(v);vec[rs]=v;status[rs]='matched'
ids=[r for r in data.rsid if r in vec and status[r]=='matched'];G=np.stack([vec[r] for r in ids]);G-=G.mean(axis=1,keepdims=True);G/=np.sqrt((G*G).sum(axis=1,keepdims=True));R=G@G.T;assert np.isfinite(R).all() and np.max(abs(np.diag(R)-1))<1e-5
np.savez_compressed(out/'stage6_rbm6_1000g_eur_reference_ld.npz',ld=R,rsid=np.array(ids)); pd.DataFrame({'rsid':data.rsid,'status':data.rsid.map(status)}).to_csv(out/'stage6_rbm6_ld_variant_coverage.csv',index=False)
summary={'n_common_gwas_snps':len(data),'n_ld':len(ids),'n_samples':len(eur),'reference_maf_min':.01,'max_missing':.01,'ref_vcf':url,'index_sha256':hashlib.sha256(index.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256((work/'1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'status_counts':pd.Series(status).value_counts().to_dict(),'ld_file': 'stage6_rbm6_1000g_eur_reference_ld.npz','reference_rank_max':502,'note':'Source VCF whole file not checksum verified; external ancestry reference, sample overlap and LD mismatch possible.'};(out/'stage6_rbm6_ld_reference_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
