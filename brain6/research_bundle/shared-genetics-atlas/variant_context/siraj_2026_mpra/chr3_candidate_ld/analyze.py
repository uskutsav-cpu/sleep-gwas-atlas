import collections,csv,hashlib,json,math,pathlib,random
import pysam
p=pathlib.Path(__file__).resolve().parent;src=p.parent/'design_matches.tsv';panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
source=list(csv.DictReader(src.open(),delimiter='\t'));cand=[r for r in source if 'chr3_' in r['region']];assert {r['rsid'] for r in cand}=={'rs73187947','rs9815484','rs1456193'}
spec=[('rs62264764',117639575,'A','G')]
for row in cand:
 v=row['variant37'];c,pos,ref,alt=v.split(':');assert c=='chr3';spec.append((row['rsid'],int(pos),alt,ref))
meta=list(csv.DictReader(panel.open(),delimiter='\t'));eur={r['sample'] for r in meta if r['super_pop']=='EUR'};assert len(eur)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);calls={};statuses=[]
for rs,pos,counted,other in spec:
 m=[x for x in vcf.fetch('3',pos-1,pos) if x.pos==pos and x.id==rs and counted in x.alleles and other in x.alleles]
 if len(m)!=1:statuses.append(dict(rsid=rs,pos=pos,status='nonunique_or_missing',matches=len(m)));continue
 x=m[0];a=x.alleles.index(counted);b=x.alleles.index(other);d={};exclude=collections.Counter()
 for sample in sorted(eur):
  if sample not in vcf.header.samples:exclude['absent']+=1;continue
  c=x.samples[sample];g=c.get('GT')
  if not g or len(g)!=2 or any(v is None for v in g) or not c.phased:exclude['missing_or_unphased']+=1;continue
  if any(v not in (a,b) for v in g):exclude['other_allele']+=1;continue
  d[sample]=tuple(int(v==a) for v in g)
 calls[rs]=d;statuses.append(dict(rsid=rs,pos=pos,status='matched',alleles=list(x.alleles),counted=counted,other=other,n=len(d),excluded=dict(exclude)))
assert 'rs62264764' in calls
def r2(a,b,c):
 den=a*(1-a)*b*(1-b)
 return (c-a*b)**2/den if den>0 else None
rows=[]
for rs,pos,counted,other in spec[1:]:
 if rs not in calls:continue
 shared=sorted(calls['rs62264764'].keys()&calls[rs].keys());n=2*len(shared);cnt=collections.Counter()
 for s in shared:
  for a,b in zip(calls['rs62264764'][s],calls[rs][s]):cnt[a,b]+=1
 f1=(cnt[1,0]+cnt[1,1])/n;f2=(cnt[0,1]+cnt[1,1])/n;v=r2(f1,f2,cnt[1,1]/n)
 rng=random.Random(20261028);per=[(sum(calls['rs62264764'][s]),sum(calls[rs][s]),sum(a*b for a,b in zip(calls['rs62264764'][s],calls[rs][s]))) for s in shared];boot=[]
 for _ in range(2000):
  d=[per[rng.randrange(len(per))] for j in per];z=r2(sum(x[0] for x in d)/n,sum(x[1] for x in d)/n,sum(x[2] for x in d)/n)
  if z is not None:boot.append(z)
 boot.sort();rows.append(dict(rsid=rs,pos_grch37=pos,counted=counted,n=len(shared),n00=cnt[0,0],n01=cnt[0,1],n10=cnt[1,0],n11=cnt[1,1],freq1=f1,freq2=f2,r2=v,bootstrap_valid=len(boot),lo=boot[int(.025*(len(boot)-1))] if boot else None,hi=boot[int(.975*(len(boot)-1))] if boot else None))
if rows:
 with (p/'ld.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
(p/'receipt.json').write_text(json.dumps(dict(url=url,whole_vcf_hashed=False,panel_sha256=sha(panel),source_sha256=sha(src),seed=20261028,bootstrap_replicates=2000,statuses=statuses),indent=2)+'\n');print('matched',len(rows),'of 3; no genotypes saved')
