#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib
import numpy as np,pandas as pd,pysam
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs';d=pd.read_csv(out/'stage6_foxp2_mvp_hfrs_exact_aligned.csv.gz'); assert d.rs.is_unique
lift=LiftOver(str(work/'hg38ToHg19.over.chain.gz')); lookup={};status={}
for r in d.itertuples():
 m=lift.convert_coordinate('chr7',int(r.pos)-1)
 if len(m)!=1 or m[0][0]!='chr7' or m[0][2]!='+':status[r.rs]='bad_liftover';continue
 key=(int(m[0][1])+1,''.join(sorted([r.anchor_effect_allele,r.anchor_other_allele])))
 if key in lookup:raise ValueError('duplicate key')
 lookup[key]=r.rs;status[r.rs]='absent_reference'
u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';index=work/'stage6_1000g_chr7_phase3.tbi';t=pysam.TabixFile(u,index=str(index));samples=[x for x in t.header if x.startswith('#CHROM')][0].split('\t')[9:];panel={x['sample']:x for x in csv.DictReader((work/'1kg_phase3_panel.tsv').open(),delimiter='\t')};idx=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR'];assert len(idx)==503
vec={};duplicate=set();anchors=d.set_index('rs').anchor_effect_allele.to_dict()
for line in t.fetch('7',min(k[0] for k in lookup)-1,max(k[0] for k in lookup)):
 z=line.split('\t');key=(int(z[1]),''.join(sorted([z[3],z[4]])));rs=lookup.get(key)
 if rs is None:continue
 if rs in duplicate:continue
 if rs in vec:duplicate.add(rs);vec.pop(rs);status[rs]='duplicate_reference_record';continue
 gi=z[8].split(':').index('GT');a=[]
 for i in idx:
  gt=z[9+i].split(':')[gi].replace('|','/').split('/');a.append(sum(int(x)==1 for x in gt) if len(gt)==2 and all(x in ('0','1') for x in gt) else np.nan)
 v=np.array(a,dtype=np.float32);v=2-v if anchors[rs]==z[3] else v; missing=np.mean(~np.isfinite(v));af=float(np.nanmean(v)/2)
 if missing>.01:status[rs]='missing_gt_1pct';continue
 if min(af,1-af)<.01:status[rs]='reference_maf_lt_1pct';continue
 v[~np.isfinite(v)]=np.nanmean(v);vec[rs]=v;status[rs]='matched'
ids=[x for x in d.rs if x in vec];g=np.stack([vec[x] for x in ids]);np.savez_compressed(work/'stage6_foxp2_1000g_eur_genotypes.npz',rsid=np.array(ids),dosage=g)
a=pd.DataFrame({'rs':d.rs,'status':d.rs.map(status)});a.to_csv(out/'stage6_foxp2_reference_coverage.csv',index=False);receipt={'source':u,'n_gwas_eligible':len(d),'n_reference':len(ids),'sample_count':503,'reference_rank_bound':502,'counts':a.status.value_counts().to_dict(),'index_sha256':hashlib.sha256(index.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256((work/'1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'whole_vcf_checksum_verified':False};(out/'stage6_foxp2_reference_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
