import csv,pysam,numpy as np,json,pathlib,hashlib
from pyliftover import LiftOver
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
variants=[('rs2838788',45122459,'T','C'),('rs188632112',45103160,'A','G'),('rs1320256433',45022799,'A','G')]
chain=pathlib.Path('work/hg38ToHg19.over.chain.gz');lo=LiftOver(str(chain));panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'};f=pysam.VariantFile(url);data={};records=[]
for rs,pos38,a,b in variants:
 maps=lo.convert_coordinate('chr21',pos38-1) or []; matches=[v for v in maps if v[0]=='chr21' and v[2]=='+']; assert len(maps)==len(matches)==1
 pos37=matches[0][1]+1; hits=[v for v in f.fetch('21',pos37-1,pos37) if v.id==rs and v.pos==pos37 and len(v.alleles)==2 and set(v.alleles)=={a,b}]
 row={'variant':rs,'GRCh38_pos':pos38,'GRCh37_pos_lifted':pos37,'counted':a,'other':b,'reference_exact_matches':len(hits)}
 if len(hits)==1:
  v=hits[0];ix=v.alleles.index(a);calls={}
  for s in sorted(eur & set(f.header.samples)):
   gt=v.samples[s].get('GT')
   if gt and len(gt)==2 and all(z is not None for z in gt) and v.samples[s].phased:calls[s]=[int(z==ix) for z in gt]
  data[rs]=calls;vals=[z for v in calls.values() for z in v];row.update(n_called=len(calls),counted_allele_frequency=sum(vals)/len(vals),other_allele_haplotypes=len(vals)-sum(vals))
 records.append(row)
pairs=[]
for rs in ('rs188632112','rs1320256433'):
 row={'variant':rs}
 if rs not in data or 'rs2838788' not in data:row['status']='reference exact match missing or ambiguous'
 else:
  ids=sorted(data[rs].keys() & data['rs2838788'].keys()); x=np.array([v for s in ids for v in data['rs2838788'][s]]);y=np.array([v for s in ids for v in data[rs][s]]);row['n_individuals']=len(ids)
  if not len(x) or x.std()==0 or y.std()==0:row['status']='no variation or calls in EUR reference; r undefined'
  else:
   r=float(np.corrcoef(x,y)[0,1]);row.update(status='matched phased EUR',r=r,r2=r*r)
 pairs.append(row)
out={'reference':url,'chain_sha256':hashlib.sha256(chain.read_bytes()).hexdigest(),'records':records,'pairs':pairs,'source_rows':{'MVP':'21 45103160 A G OR 1.507 EAF 0.9988 p 0.0001799 rs188632112 imputation_r2 0.652','FinnGen':'21 45022799 ref G alt A rs1320256433 p 1.15854e-6 beta_alt 0.320041 se 0.0658169 af_alt 0.00210788'},'caveat':'Reference rare-allele sampling and source imputation limit LD estimation; no identity or causal inference. MVP OR rounded, EAF near 1.'}
pathlib.Path('outputs/adarb1_insomnia_regional_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
