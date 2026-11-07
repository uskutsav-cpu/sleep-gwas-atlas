import pathlib,csv,json,collections,hashlib,math,random,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;src=p.parent
paths={t:src/('Adrenal_Gland_ENSG00000185565_'+t+'_matched.tsv') for t in ('insomnia','adhd')}
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
a={k:rows(v) for k,v in paths.items()};q={t:min(x,key=lambda z:float(z['qtl_p'])) for t,x in a.items()};assert q['insomnia']['variant']==q['adhd']['variant']
selected=[('eqtl_lead',q['insomnia']),('insomnia_lead',min(a['insomnia'],key=lambda z:float(z['gwas_p_from_beta_se']))),('adhd_lead',min(a['adhd'],key=lambda z:float(z['gwas_p_from_beta_se'])))]
anchor=[x for x in a['insomnia'] if x['variant']=='chr3_117920728_G_A'];assert len(anchor)==1;selected.append(('alphagenome_candidate',anchor[0]));assert len({x['variant'] for _,x in selected})==4
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');chain=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
with panel.open() as f:eur={x['sample'] for x in csv.DictReader(f,delimiter='\t') if x['super_pop']=='EUR'}
assert len(eur)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';v=pysam.VariantFile(url);lo=LiftOver(str(chain));calls={};metadata={}
for role,x in selected:
 pos=int(x['pos38']);m=lo.convert_coordinate('chr3',pos-1);assert len(m)==1 and m[0][0]=='chr3' and m[0][2]=='+';pos37=m[0][1]+1;rs=x['gwas_rsid'];ref=x['ref'];alt=x['alt']
 candidates=[q for q in v.fetch('3',pos37-1,pos37) if q.pos==pos37 and q.id and rs in q.id.split(';') and ref in q.alleles and alt in q.alleles];assert len(candidates)==1,(role,rs,len(candidates));q=candidates[0];ai=q.alleles.index(alt);ri=q.alleles.index(ref);data={};excluded=collections.Counter()
 for sample in sorted(eur & set(v.header.samples)):
  s=q.samples[sample];gt=s.get('GT')
  if not gt or len(gt)!=2 or not s.phased or any(g not in (ri,ai) for g in gt):excluded['missing_unphased_or_other']+=1;continue
  data[sample]=tuple(int(g==ai) for g in gt)
 calls[role]=data;metadata[role]={'rsid':rs,'grch38_pos':pos,'grch37_pos':pos37,'ref':ref,'counted_allele':alt,'vcf_alleles':list(q.alleles),'people':len(data),'excluded':dict(excluded)}
result=[]
for i,(one,_) in enumerate(selected):
 for two,_ in selected[i+1:]:
  shared=sorted(calls[one].keys() & calls[two].keys());pairs=[(sum(calls[one][s]),sum(calls[two][s]),sum(a*b for a,b in zip(calls[one][s],calls[two][s]))) for s in shared]
  def corr(xs):
   n=2*len(xs);a=sum(x[0] for x in xs)/n;b=sum(x[1] for x in xs)/n;c=sum(x[2] for x in xs)/n;d=math.sqrt(a*(1-a)*b*(1-b));return (c-a*b)/d if d else None
  r=corr(pairs);assert r is not None;rng=random.Random(20261027);bs=[]
  for _ in range(2000):
   x=corr([pairs[rng.randrange(len(pairs))] for _ in pairs]);
   if x is not None:bs.append(x*x)
  assert len(bs)==2000;bs.sort();result.append({'role_1':one,'rsid_1':metadata[one]['rsid'],'allele_1':metadata[one]['counted_allele'],'role_2':two,'rsid_2':metadata[two]['rsid'],'allele_2':metadata[two]['counted_allele'],'people':len(shared),'r':r,'r2':r*r,'r2_boot_025':bs[49],'r2_boot_975':bs[1949]})
with (p/'pairwise_ld.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'input_sha256':{k:sha(v) for k,v in paths.items()},'chain_sha256':sha(chain),'panel_sha256':sha(panel),'reference_url':url,'variants':metadata,'selected_pvalues':{role:{'qtl_p':x['qtl_p'],'gwas_p_from_beta_se':x['gwas_p_from_beta_se']} for role,x in selected},'bootstrap_replicates':2000,'seed':20261027,'output_sha256':sha(p/'pairwise_ld.tsv'),'note':'post-hoc, EUR reference; remote whole VCF unverified'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))
