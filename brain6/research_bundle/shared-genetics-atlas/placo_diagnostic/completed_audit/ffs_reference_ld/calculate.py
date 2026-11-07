import collections,csv,hashlib,json,math,pathlib,sys,pysam
p=pathlib.Path(__file__).resolve().parent; selection=p.parent/'candidate_positions.tsv'
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
chrom=sys.argv[1]; assert chrom in {'5','11','16','17','22'}
with panel.open() as f: eur={x['sample'] for x in csv.DictReader(f,delimiter='\t') if x['super_pop']=='EUR'}
assert len(eur)==503
with selection.open() as f: allrows=list(csv.DictReader(f,delimiter='\t'))
rows=[x for x in allrows if x['CHR']==chrom and (x['pair']=='FFS' or chrom=='22' and x['pair']=='FI')]
ffs=[x for x in rows if x['pair']=='FFS']; assert len(ffs)=={'5':27,'11':9,'16':4,'17':14,'22':12}[chrom]
variants={x['SNP']:x for x in rows};assert len(variants)==len(rows)
url=f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);calls={};meta={}
for snp,x in variants.items():
 pos=int(x['BP_GRCH37']);a=x['A1'];b=x['A2']
 rec=[v for v in vcf.fetch(chrom,pos-1,pos) if v.pos==pos and v.id==snp and a in v.alleles and b in v.alleles]
 info={'pos':pos,'A1':a,'A2':b,'matches':len(rec)};meta[snp]=info
 if len(rec)!=1:continue
 v=rec[0];ai=v.alleles.index(a);bi=v.alleles.index(b);c={};excluded=collections.Counter()
 for sample in sorted(eur & set(vcf.header.samples)):
  call=v.samples[sample];g=call.get('GT')
  if not g or len(g)!=2 or any(z is None for z in g) or not call.phased:excluded['missing_or_unphased']+=1;continue
  if any(z not in (ai,bi) for z in g):excluded['other_allele']+=1;continue
  c[sample]=tuple(int(z==ai) for z in g)
 calls[snp]=c;info.update(alleles=list(v.alleles),eligible=len(c),excluded=dict(excluded))
result=[]
for i,x in enumerate(rows):
 for y in rows[i+1:]:
  if x['pair']!=y['pair'] and chrom!='22':continue
  if x['pair']==y['pair']=='FI':continue
  a,b=x['SNP'],y['SNP'];r={'variant_1':a,'variant_2':b,'pair_1':x['pair'],'pair_2':y['pair'],'A1_1':x['A1'],'A1_2':y['A1'],'n_individuals':0,'r':'','r2':'','status':'MISSING_REFERENCE'}
  if a in calls and b in calls:
   samples=set(calls[a])&set(calls[b]); n=2*len(samples);r['n_individuals']=len(samples)
   if n:
    af=sum(sum(calls[a][s]) for s in samples)/n;bf=sum(sum(calls[b][s]) for s in samples)/n;ab=sum(sum(k*l for k,l in zip(calls[a][s],calls[b][s])) for s in samples)/n;den=math.sqrt(af*(1-af)*bf*(1-bf))
    if den:
     z=(ab-af*bf)/den;r.update(r=z,r2=z*z,status='OK')
    else:r['status']='MONOMORPHIC'
  result.append(r)
with (p/f'chr{chrom}_pairs.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
anchor=min(ffs,key=lambda x:float(x['p_placo_plus']))['SNP']
receipt={'reference':url,'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'selected_sha256':hashlib.sha256(selection.read_bytes()).hexdigest(),'anchor':anchor,'variants':meta,'pair_counts':dict(collections.Counter(x['status'] for x in result)),'limitations':'Descriptive post-selected EUR reference LD only; not independent loci or replication.'}
(p/f'chr{chrom}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(chrom,'anchor',anchor,'n variants',len(variants),'reference matches',len(calls),'pair counts',receipt['pair_counts'])
