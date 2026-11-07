import csv,json,pathlib,math,hashlib,collections,random,pysam
p=pathlib.Path(__file__).resolve().parent
sel=p.parent/'candidate_positions.tsv';panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
rows=list(csv.DictReader(panel.open(),delimiter='\t'));eur={r['sample'] for r in rows if r['super_pop']=='EUR'};assert len(eur)==503
with sel.open() as f: c=[r for r in csv.DictReader(f,delimiter='\t') if r['CHR']=='7']
assert len(c)==7 and all(r['position_status']=='UNIQUE_RSID_ALLELE_MATCH' for r in c)
variants={r['SNP']:('7',int(r['BP_GRCH37']),r['A1'],r['A2']) for r in c};variants.update({'rs1476535':('7',114071035,'C','T'),'rs12536335':('7',114043159,'G','A')})
pairs=[]
for r in c:
 for anchor in ('rs2396766','rs1476535','rs12536335'):
  if r['SNP']!=anchor:pairs.append((r['SNP'],anchor))
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);calls={};meta={}
for rs,(chrom,pos,counted,other) in variants.items():
 rec=[x for x in vcf.fetch(chrom,pos-1,pos) if x.pos==pos and x.id==rs and counted in x.alleles and other in x.alleles]
 meta[rs]={'position_grch37':pos,'counted':counted,'other':other,'matching_reference_records':len(rec)}
 if len(rec)!=1:continue
 x=rec[0];a=x.alleles.index(counted);b=x.alleles.index(other);paired={};excluded=collections.Counter()
 for sample in sorted(eur&set(vcf.header.samples)):
  call=x.samples[sample];g=call.get('GT')
  if not g or len(g)!=2 or any(z is None for z in g) or not call.phased:excluded['missing_or_unphased']+=1;continue
  if any(z not in (a,b) for z in g):excluded['other_allele_genotype']+=1;continue
  paired[sample]=tuple(int(z==a) for z in g)
 calls[rs]=paired;meta[rs].update(reference_alleles=list(x.alleles),eligible=len(paired),excluded=dict(excluded))
results=[]
for first,second in pairs:
 if first not in calls or second not in calls:results.append({'variant_1':first,'variant_2':second,'status':'MISSING_REFERENCE'});continue
 samples=sorted(set(calls[first])&set(calls[second]));per=[(sum(calls[first][s]),sum(calls[second][s]),sum(a*b for a,b in zip(calls[first][s],calls[second][s]))) for s in samples];n=2*len(samples)
 def calc(pairs):
  a=sum(t[0] for t in pairs)/n;b=sum(t[1] for t in pairs)/n;c=sum(t[2] for t in pairs)/n;den=math.sqrt(a*(1-a)*b*(1-b));return (c-a*b)/den if den else None
 corr=calc(per);rng=random.Random(20261027);boot=[]
 for _ in range(1000):
  z=calc([per[rng.randrange(len(per))] for _ in per]);
  if z is not None:boot.append(z*z)
 boot.sort();results.append({'variant_1':first,'variant_2':second,'status':'OK' if corr is not None else 'MONOMORPHIC','counted_allele_1':variants[first][2],'counted_allele_2':variants[second][2],'n_individuals':len(samples),'r':corr,'r2':corr*corr if corr is not None else None,'bootstrap_r2_025':boot[24] if len(boot)==1000 else None,'bootstrap_r2_975':boot[974] if len(boot)==1000 else None,'bootstrap_n':len(boot)})
fields=['variant_1','variant_2','status','counted_allele_1','counted_allele_2','n_individuals','r','r2','bootstrap_r2_025','bootstrap_r2_975','bootstrap_n']
with (p/'pairs.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(results)
receipt={'source':url,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'selection_sha256':hashlib.sha256(sel.read_bytes()).hexdigest(),'variants':meta,'n_pairs':len(results),'n_ok':sum(r['status']=='OK' for r in results),'caveat':'Post-hoc reference LD is not shared causality or independent replication.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
