import csv,json,hashlib,pathlib,math,collections,pysam
p=pathlib.Path(__file__).resolve().parent;sel=p.parent.parent/'candidate_positions.tsv';panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
anchors={'rs2431108':(103947968,'T','C'),'rs1363103':(103917837,'T','C'),'rs13164188':(103914523,'T','C'),'rs7703746':(103907337,'G','A')}
with sel.open() as f:rows=[x for x in csv.DictReader(f,delimiter='\t') if x['pair']=='FFS' and x['CHR']=='5']
assert len(rows)==27 and len({x['SNP'] for x in rows})==27
variants={x['SNP']:(int(x['BP_GRCH37']),x['A1'],x['A2']) for x in rows}
for k,v in anchors.items():
 if k in variants: assert variants[k][0]==v[0] and set(variants[k][1:])==set(v[1:])
 else:variants[k]=v
with panel.open() as f: samples=list(csv.DictReader(f,delimiter='\t'))
pops={'EUR':{r['sample'] for r in samples if r['super_pop']=='EUR'},'FIN':{r['sample'] for r in samples if r['pop']=='FIN'}}
assert len(pops['EUR'])==503 and pops['FIN']<=pops['EUR']
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);calls={};metadata={}
for rs,(pos,a,b) in variants.items():
 hits=[x for x in vcf.fetch('5',pos-1,pos) if x.pos==pos and x.id==rs and a in x.alleles and b in x.alleles];meta={'position':pos,'counted':a,'other':b,'matches':len(hits)};metadata[rs]=meta
 if len(hits)!=1:continue
 rec=hits[0];ia=rec.alleles.index(a);ib=rec.alleles.index(b);d={};skip=collections.Counter()
 for sample in pops['EUR']&set(vcf.header.samples):
  call=rec.samples[sample];g=call.get('GT')
  if not g or len(g)!=2 or any(x is None for x in g) or not call.phased:skip['missing_or_unphased']+=1;continue
  if any(x not in (ia,ib) for x in g):skip['other_allele']+=1;continue
  d[sample]=tuple(int(x==ia) for x in g)
 calls[rs]=d;meta.update(reference_alleles=list(rec.alleles),eligible=len(d),excluded=dict(skip))
result=[]
for row in rows:
 for target in anchors:
  for pop in pops:
   a=row['SNP'];b=target;out={'FFS_marker':a,'anchor':b,'population':pop,'same_variant':a==b,'FFS_counted_allele':variants[a][1],'anchor_counted_allele':variants[b][1],'n_individuals':0,'r':'','r2':'','status':'MISSING_REFERENCE'}
   if a in calls and b in calls:
    common=pops[pop]&calls[a].keys()&calls[b].keys();n=2*len(common);out['n_individuals']=len(common)
    if n:
     x=sum(sum(calls[a][s]) for s in common)/n;y=sum(sum(calls[b][s]) for s in common)/n;z=sum(sum(i*j for i,j in zip(calls[a][s],calls[b][s])) for s in common)/n;den=math.sqrt(x*(1-x)*y*(1-y))
     if den:
      r=(z-x*y)/den;out.update(r=r,r2=r*r,status='OK')
     else:out['status']='MONOMORPHIC'
   result.append(out)
with (p/'pairs.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
receipt={'selection_sha256':hashlib.sha256(sel.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'vcf':url,'populations':{k:len(v) for k,v in pops.items()},'variants':metadata,'status_counts':dict(collections.Counter(x['status'] for x in result)),'limitation':'Post-selected LD in nested reference samples is not evidence of causal colocalization or independent replication.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('population sizes',receipt['populations'],'pair statuses',receipt['status_counts'],'records',len(result))
