import csv,json,pathlib,math,hashlib,collections,random,pysam
p=pathlib.Path(__file__).resolve().parent
posterior_file=p.parent/'frailty_coloc/conditional_variant_posteriors.tsv'
with posterior_file.open() as f: posterior=list(csv.DictReader(f,delimiter='\t'))
chr7=sorted((x for x in posterior if x['chr']=='7'),key=lambda x:-float(x['conditional_h4']))
assert [x['rsid'] for x in chr7[:2]]==['rs10953765','rs1456031']
for x in chr7[:2]:
 assert (x['bp37'],x['effect_allele'],x['other_allele']) in [('114291435','A','G'),('114296102','T','C')]
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');rows=list(csv.DictReader(panel.open(),delimiter='\t'));eur={r['sample'] for r in rows if r['super_pop']=='EUR'};assert len(eur)==503
variants={'rs10953765':('7',114291435,'A','G'),'rs1456031':('7',114296102,'T','C'),'rs1476535':('7',114071035,'C','T'),'rs12536335':('7',114043159,'G','A')};pairs=[('rs10953765','rs1456031'),('rs10953765','rs1476535'),('rs10953765','rs12536335'),('rs1456031','rs1476535'),('rs1456031','rs12536335')];base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';calls={};meta={}
for rs,(chrom,pos,counted,other) in variants.items():
 vcf=pysam.VariantFile(base.format(chrom=chrom));records=[x for x in vcf.fetch(chrom,pos-1,pos) if x.pos==pos and x.id==rs];matches=[x for x in records if counted in x.alleles and other in x.alleles];assert len(matches)==1,(rs,len(records),len(matches));x=matches[0];a=x.alleles.index(counted);b=x.alleles.index(other);paired={};counts=collections.Counter()
 for s in sorted(eur & set(vcf.header.samples)):
  call=x.samples[s];g=call.get('GT')
  if not g or len(g)!=2 or any(z is None for z in g) or not call.phased:counts['missing_or_unphased']+=1;continue
  if any(z not in (a,b) for z in g):counts['other_allele_genotype']+=1;continue
  paired[s]=tuple(int(z==a) for z in g)
 calls[rs]=paired;meta[rs]={'chrom':chrom,'position_grch37':pos,'counted':counted,'other':other,'reference_alleles':list(x.alleles),'samples_eligible':len(paired),'excluded':dict(counts),'url':base.format(chrom=chrom)};print(rs,meta[rs],flush=True)
res=[]
for one,two in pairs:
 shared=sorted(set(calls[one])&set(calls[two])); x=[];y=[]
 for s in shared:x.extend(calls[one][s]);y.extend(calls[two][s])
 n=len(x);p1=sum(x)/n;p2=sum(y)/n;p11=sum(a*b for a,b in zip(x,y))/n;den=math.sqrt(p1*(1-p1)*p2*(1-p2));r=(p11-p1*p2)/den if den else None
 rng=random.Random(20261027);per=[]
 for s in shared:
  a=calls[one][s];b=calls[two][s];per.append((sum(a),sum(b),a[0]*b[0]+a[1]*b[1]))
 bootstrap=[]
 for _ in range(2000):
  draws=[per[rng.randrange(len(per))] for j in range(len(per))];a=sum(z[0] for z in draws)/n;b=sum(z[1] for z in draws)/n;c=sum(z[2] for z in draws)/n;dd=math.sqrt(a*(1-a)*b*(1-b));
  if dd:bootstrap.append(((c-a*b)/dd)**2)
 bootstrap.sort();assert len(bootstrap)==2000;res.append({'variant_1':one,'counted_allele_1':variants[one][2],'variant_2':two,'counted_allele_2':variants[two][2],'individuals':len(shared),'haplotypes':n,'allele1_frequency':p1,'allele2_frequency':p2,'p11':p11,'r':r,'r2':r*r if r is not None else None,'bootstrap_r2_025':bootstrap[49],'bootstrap_r2_975':bootstrap[1949],'bootstrap_replicates':2000})
with (p/'chr7_coloc_top_ld.tsv').open('w',newline='') as f:w=csv.DictWriter(f,res[0],delimiter='\t');w.writeheader();w.writerows(res)
(p/'chr7_coloc_top_receipt.json').write_text(json.dumps({'variants':meta,'pairs':res,'posterior_sha256':hashlib.sha256(posterior_file.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'panel_n_EUR':len(eur),'vcf_release':'1000 Genomes Phase 3 GRCh37, UCSC-hosted indexed reference; remote indexed access, no whole-file hash','method':'phased allele indicators; pairwise-complete samples with both called haplotypes and listed allele contrasts; percentile bootstrap resamples individuals 2000 times with seed 20261027'},indent=2)+'\n');print(json.dumps(res,indent=2))
