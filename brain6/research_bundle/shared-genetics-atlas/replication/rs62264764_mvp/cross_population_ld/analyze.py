import collections,csv,hashlib,json,math,pathlib,random
import pysam
p=pathlib.Path(__file__).resolve().parent
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
meta=list(csv.DictReader(panel.open(),delimiter='\t'));groups={g:{r['sample'] for r in meta if r['super_pop']==g} for g in sorted({r['super_pop'] for r in meta})}
assert set(groups)=={'AFR','AMR','EAS','EUR','SAS'} and sum(map(len,groups.values()))==2504
spec=[('rs62264764',117639575,'A','G'),('rs62264767',117642005,'C','A')]
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);calls=[];info=[]
for rs,pos,allele,other in spec:
 matches=[x for x in vcf.fetch('3',pos-1,pos) if x.pos==pos and x.id==rs and allele in x.alleles and other in x.alleles];assert len(matches)==1
 rec=matches[0];ia=rec.alleles.index(allele);ib=rec.alleles.index(other);calls.append({});excluded=collections.Counter()
 for group,ids in groups.items():
  for sample in sorted(ids):
   if sample not in vcf.header.samples:excluded[group+':absent']+=1;continue
   c=rec.samples[sample];gt=c.get('GT')
   if not gt or len(gt)!=2 or any(t is None for t in gt) or not c.phased:excluded[group+':missing_or_unphased']+=1;continue
   if any(t not in (ia,ib) for t in gt):excluded[group+':other_allele']+=1;continue
   calls[-1][sample]=tuple(int(t==ia) for t in gt)
 info.append(dict(rsid=rs,position_grch37=pos,counted=allele,other=other,source_alleles=list(rec.alleles),exclusions=dict(excluded)))
rows=[]
for group,ids in groups.items():
 shared=sorted(ids & calls[0].keys() & calls[1].keys());counts=collections.Counter()
 for sample in shared:
  for a,b in zip(calls[0][sample],calls[1][sample]):counts[(a,b)]+=1
 n=2*len(shared);p1=(counts[1,0]+counts[1,1])/n;p2=(counts[0,1]+counts[1,1])/n;p11=counts[1,1]/n
 def r2(a,b,c):
  den=a*(1-a)*b*(1-b)
  return (c-a*b)**2/den if den>0 else None
 v=r2(p1,p2,p11);rng=random.Random(20261028);per=[(sum(calls[0][s]),sum(calls[1][s]),sum(a*b for a,b in zip(calls[0][s],calls[1][s]))) for s in shared];boot=[]
 for i in range(2000):
  draw=[per[rng.randrange(len(per))] for j in per];a=sum(x[0] for x in draw)/n;b=sum(x[1] for x in draw)/n;c=sum(x[2] for x in draw)/n;z=r2(a,b,c)
  if z is not None:boot.append(z)
 boot.sort();row=dict(group=group,panel_n=len(ids),complete_n=len(shared),haplotypes=n,n00=counts[0,0],n01=counts[0,1],n10=counts[1,0],n11=counts[1,1],freq1=p1,freq2=p2,r2=v,bootstrap_valid=len(boot),bootstrap_025=boot[int(.025*(len(boot)-1))] if boot else None,bootstrap_975=boot[int(.975*(len(boot)-1))] if boot else None)
 assert sum(row[k] for k in ('n00','n01','n10','n11'))==n;rows.append(row)
 if group=='EUR':assert len(shared)==503 and abs(v-1)<1e-12
with (p/'group_ld.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps(dict(url=url,source='1000 Genomes Phase 3 phased GRCh37, indexed VCF',whole_file_hashed=False,panel_sha256=hashlib.sha256(panel.read_bytes()).hexdigest(),pysam=pysam.__version__,variants=info,bootstrap_seed=20261028,bootstrap_replicates=2000,groups=rows),indent=2)+'\n')
print('PASS: exact allele matches, EUR reproduction and aggregated haplotype arithmetic; no genotypes saved')
