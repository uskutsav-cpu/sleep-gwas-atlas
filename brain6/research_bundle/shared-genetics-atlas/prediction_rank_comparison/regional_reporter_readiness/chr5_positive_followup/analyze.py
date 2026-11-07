import collections,csv,hashlib,json,math,pathlib,random
import pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;atlas=p.parents[2];panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');chain=atlas/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz'
meta=list(csv.DictReader(panel.open(),delimiter='\t'));eur={r['sample'] for r in meta if r['super_pop']=='EUR'};assert len(eur)==503
spec=[('rs12659965',103904399,'A','G'),('rs2431108',103947968,'C','T'),('rs77960',103964585,'A','G')];url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);calls={};statuses=[]
for rs,pos,a,b in spec:
 matches=[x for x in vcf.fetch('5',pos-1,pos) if x.pos==pos and x.id==rs and a in x.alleles and b in x.alleles];assert len(matches)==1
 x=matches[0];ia=x.alleles.index(a);ib=x.alleles.index(b);d={};ex=collections.Counter()
 for s in sorted(eur):
  c=x.samples[s];g=c.get('GT')
  if not g or len(g)!=2 or any(z is None for z in g) or not c.phased:ex['missing_unphased']+=1;continue
  if any(z not in (ia,ib) for z in g):ex['other_allele']+=1;continue
  d[s]=tuple(int(z==ia) for z in g)
 calls[rs]=d;statuses.append(dict(rsid=rs,position37=pos,counted=a,other=b,reference_alleles=list(x.alleles),n=len(d),excluded=dict(ex)))
rows=[]
for rs,_,_,_ in spec[1:]:
 shared=sorted(calls[spec[0][0]].keys()&calls[rs].keys());n=2*len(shared);co=collections.Counter()
 for s in shared:
  for a,b in zip(calls[spec[0][0]][s],calls[rs][s]):co[a,b]+=1
 def r2(c):
  a,b,d=(sum(x[i] for x in c)/(2*len(c)) for i in range(3));den=a*(1-a)*b*(1-b);return (d-a*b)**2/den if den else None
 per=[(sum(calls[spec[0][0]][s]),sum(calls[rs][s]),sum(a*b for a,b in zip(calls[spec[0][0]][s],calls[rs][s]))) for s in shared];v=r2(per);rng=random.Random(20261028);boot=[]
 for _ in range(2000):
  draw=[per[rng.randrange(len(per))] for j in per];z=r2(draw)
  if z is not None:boot.append(z)
 boot.sort();rows.append(dict(lead=rs,n=len(shared),n00=co[0,0],n01=co[0,1],n10=co[1,0],n11=co[1,1],r2=v,lo=boot[int(.025*(len(boot)-1))],hi=boot[int(.975*(len(boot)-1))],bootstrap_valid=len(boot)))
with (p/'ld.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
mapping=LiftOver(str(chain)).convert_coordinate('chr5',103904398);assert len(mapping)==1 and mapping[0][0]=='chr5' and mapping[0][2]=='+';pos38=mapping[0][1]+1
src=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');base='placo_insomnia__adhd_chr5_103447968_104447968__';paths={'brain_insomnia':src/(base+'insomnia.tsv'),'brain_adhd':src/(base+'adhd.tsv'),'FI':atlas/'frailty_coloc_chr5/FI_extracted.tsv','R13':atlas/'frailty_coloc_chr5/R13_extracted.tsv'}
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();reports={}
for name,path in paths.items():
 with path.open() as f:data=list(csv.DictReader(f,delimiter='\t'))
 idkey='rsids' if name=='R13' else 'SNP';poskey='pos' if name=='R13' else ('BP' if name=='FI' else 'POS');position=pos38 if name=='R13' else 103904399
 matches=[r for r in data if r[idkey]=='rs12659965' or r[poskey]==str(position)]
 reports[name]=dict(path=str(path),sha256=sha(path),total_rows=len(data),rows=matches)
(p/'receipt.json').write_text(json.dumps(dict(reference_url=url,whole_vcf_hashed=False,reference_status=statuses,panel_sha256=sha(panel),chain_sha256=sha(chain),mapping38=pos38,source_reports=reports,bootstrap_seed=20261028),indent=2)+'\n');print('GRCh38',pos38,'matches',[(k,len(v['rows'])) for k,v in reports.items()]);print('no genotypes saved')
