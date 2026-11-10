import csv,pathlib,pysam,numpy as np,json,hashlib,itertools
score=pathlib.Path('outputs/unvalidated_22_shared_direction_score.csv')
rows=list(csv.DictReader(score.open())); panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'))
eur={x['sample'] for x in panel if x['super_pop']=='EUR'}
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
files={}; data=[]
for row in rows:
 c=row['chr_GRCh37_assumed']; p=int(row['position_GRCh37_assumed']); a=row['counted_allele']; b=row['other_allele']
 if c not in files: files[c]=pysam.VariantFile(base.format(chrom=c))
 f=files[c]; matches=[v for v in f.fetch(c,p-1,p) if v.id==row['variant'] and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={a,b}]; assert len(matches)==1
 v=matches[0]; ix=v.alleles.index(a); d={}
 for s in eur & set(f.header.samples):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(k is not None for k in gt): d[s]=sum(k==ix for k in gt)
 data.append(d)
 print(row['variant'],flush=True)
samples=sorted(set.intersection(*(set(d) for d in data))); assert len(samples)==503
X=np.array([[d[s] for d in data] for s in samples],float); full=X.sum(axis=1); groups={}
for i,r in enumerate(rows): groups.setdefault(r['chr_GRCh37_assumed'],[]).append(i)
# Same-chromosome components joined by pairs within 1 Mb; a descriptive grouping, not independent loci.
components=[]
for chrom,idx in groups.items():
 unseen=set(idx)
 while unseen:
  todo=[unseen.pop()]; component=[]
  while todo:
   i=todo.pop(); component.append(i)
   neighbors=[j for j in unseen if abs(int(rows[i]['position_GRCh37_assumed'])-int(rows[j]['position_GRCh37_assumed']))<=1000000]
   for j in neighbors: unseen.remove(j);todo.append(j)
  components.append(component)
def midrank(a):
 order=np.argsort(a,kind='stable'); ranks=np.empty(len(a),float); start=0
 while start<len(a):
  end=start+1
  while end<len(a) and a[order[end]]==a[order[start]]: end+=1
  ranks[order[start:end]]=(start+1+end)/2;start=end
 return ranks
results=[]
for comp in components:
 removed=X[:,comp].sum(axis=1); rest=full-removed
 results.append(dict(chromosome=rows[comp[0]]['chr_GRCh37_assumed'],variants=';'.join(rows[i]['variant'] for i in comp),positional_genes=';'.join(rows[i]['nearest_gene_body_positional'] for i in comp),n_markers=len(comp),removed_mean=float(removed.mean()),removed_sd=float(removed.std(ddof=1)),remaining_sd=float(rest.std(ddof=1)),correlation_with_full=float(np.corrcoef(full,rest)[0,1]),spearman_with_full=float(np.corrcoef(midrank(full),midrank(rest))[0,1])))
results.sort(key=lambda r:r['correlation_with_full'])
with open('outputs/shared22_leave_region_reference.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
assert sum(r['n_markers'] for r in results)==22 and len(results)==18
print('full SD',full.std(ddof=1),'components',len(results));print(*results[:5],sep='\n')
