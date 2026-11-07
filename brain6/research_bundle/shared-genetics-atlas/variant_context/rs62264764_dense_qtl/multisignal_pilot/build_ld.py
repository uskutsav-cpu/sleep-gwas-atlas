import pathlib,csv,hashlib,collections,json,math,gzip
import numpy as np,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;par=p.parent
paths={t:par/('Adrenal_Gland_ENSG00000185565_'+t+'_matched.tsv') for t in ('insomnia','adhd')}
def read(f):
 with f.open() as x:return list(csv.DictReader(x,delimiter='\t'))
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
a={t:read(f) for t,f in paths.items()};indices={t:{x['variant']:x for x in r} for t,r in a.items()};assert all(len(indices[t])==len(a[t]) for t in a)
common=sorted(set(indices['insomnia'])&set(indices['adhd']));assert len(common)==1556
chain=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz');panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');lo=LiftOver(str(chain))
with panel.open() as f:eur=sorted(x['sample'] for x in csv.DictReader(f,delimiter='\t') if x['super_pop']=='EUR')
assert len(eur)==503
positions=collections.defaultdict(list);counts=collections.Counter()
for key in common:
 x=indices['insomnia'][key];other=indices['adhd'][key]
 assert (x['ref'],x['alt'],x['pos38'],x['qtl_beta'],x['qtl_se'])==(other['ref'],other['alt'],other['pos38'],other['qtl_beta'],other['qtl_se'])
 m=lo.convert_coordinate('chr3',int(x['pos38'])-1)
 if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':counts['lift_failed']+=1;continue
 positions[m[0][1]+1].append(key)
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);assert set(eur).issubset(set(vcf.header.samples));records=collections.defaultdict(list)
for q in vcf.fetch('3',min(positions)-1,max(positions)):
 if q.pos not in positions:continue
 for key in positions[q.pos]:
  x=indices['insomnia'][key];rs=x['gwas_rsid']
  if q.id and rs in q.id.split(';') and x['ref'] in q.alleles and x['alt'] in q.alleles:records[key].append(q)
retained=[];genotypes=[]
for key in common:
 if not any(key in v for v in positions.values()):continue
 records_for_key=records.get(key,[])
 if len(records_for_key)!=1:counts['no_unique_vcf_record']+=1;continue
 q=records_for_key[0];x=indices['insomnia'][key];ref,alt=x['ref'],x['alt'];ri=q.alleles.index(ref);ai=q.alleles.index(alt);arr=[];bad=False
 for person in eur:
  gt=q.samples[person].get('GT')
  if not gt or len(gt)!=2 or any(g not in (ri,ai) for g in gt):bad=True;break
  arr.append(sum(g==ai for g in gt))
 if bad:counts['missing_or_other_allele']+=1;continue
 freq=sum(arr)/(2*len(arr));maf=min(freq,1-freq)
 if maf<.01:counts['reference_maf_below_001']+=1;continue
 geno=np.asarray(arr,dtype=np.float64)
 if geno.std()==0:counts['monomorphic']+=1;continue
 y=indices['adhd'][key];retained.append({'variant':key,'rsid':x['gwas_rsid'],'pos38':x['pos38'],'pos37':q.pos,'ref':ref,'alt':alt,'reference_alt_frequency':freq,'qtl_beta':x['qtl_beta'],'qtl_se':x['qtl_se'],'insomnia_beta':x['gwas_beta_alt'],'insomnia_se':x['gwas_se'],'adhd_beta':y['gwas_beta_alt'],'adhd_se':y['gwas_se']});genotypes.append(geno)
counts['retained']=len(retained);assert sum(counts.values())==len(common)
G=np.array(genotypes);assert G.shape==(len(retained),503);G-=G.mean(axis=1,keepdims=True);G/=np.linalg.norm(G,axis=1,keepdims=True);ld=G@G.T
assert np.isfinite(ld).all() and np.max(np.abs(ld-ld.T))<1e-12 and np.max(np.abs(np.diag(ld)-1))<1e-12
with (p/'variants.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(retained[0]),delimiter='\t');w.writeheader();w.writerows(retained)
with gzip.open(p/'ld_matrix.tsv.gz','wt') as f:np.savetxt(f,ld,delimiter='\t',fmt='%.12g')
receipt={'source_input_hashes':{k:sha(f) for k,f in paths.items()},'panel_sha256':sha(panel),'chain_sha256':sha(chain),'reference_url':url,'input_intersection':len(common),'counts':dict(counts),'retention':len(retained)/len(common),'ld_shape':list(ld.shape),'n_individuals':len(eur),'max_abs_diagonal_error':float(np.max(np.abs(np.diag(ld)-1))),'max_abs_asymmetry':float(np.max(np.abs(ld-ld.T))),'sample_rank_upper_bound':len(eur)-1,'variants_sha256':sha(p/'variants.tsv'),'ld_sha256':sha(p/'ld_matrix.tsv.gz'),'caveat':'non-study-matched out-of-sample LD; no individual genotype output; no source EAF; no independent replication'}
(p/'ld_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
