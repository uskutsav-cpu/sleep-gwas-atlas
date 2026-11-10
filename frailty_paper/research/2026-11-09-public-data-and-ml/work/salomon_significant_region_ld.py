import csv,gzip,json,hashlib,collections
from pathlib import Path
import numpy as np,pysam
from pyliftover import LiftOver
out=Path('outputs');src=Path('work/salomon_neuronal_variant_readout.tsv.gz'); lo=LiftOver('work/hg38ToHg19.over.chain.gz')
names=['rs12536335','rs12671330','rs17312861','rs1942262','rs11152363'];candidates={}
for fn in ['integrated_15_candidate_review.csv','expanded_22_candidate_risk_gene_table.csv']:
 for r in csv.DictReader((out/fn).open()):
  if r['variant'] not in names or r['variant'] in candidates:continue
  candidates[r['variant']]={'chrom':r.get('GRCh37_chr',r.get('GRCh37_chr_assumed')),'pos37':int(r.get('GRCh37_pos',r.get('GRCh37_pos_assumed'))),'pos38':int(r['GRCh38_pos']),'alleles':{r['counted_allele'],r['other_allele']}}
assert len(candidates)==5
significant={}
with gzip.open(src,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['NGN2_variant_is_significant']!='True':continue
  fields=r['chr_pos_ref_alt'].split('-')
  if len(fields)!=4:continue
  c,p,a,b=fields
  try:p=int(p)
  except ValueError:continue
  matched=[name for name,t in candidates.items() if c==t['chrom'] and abs(t['pos38']-p)<=1000000]
  if matched:
   k=r['chr_pos_ref_alt'];s=significant.setdefault(k,{'chrom':c,'pos38':p,'alleles':{a,b},'candidates':set(),'adjusted_p':set(),'genes':set()});s['candidates'].update(matched);s['adjusted_p'].add(r['NGN2_variant_adj_p-value']);s['genes'].add(r['associated_gene_name'])
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));samples=sorted(r['sample'] for r in panel if r['super_pop']=='EUR');assert len(samples)==503;handles={}
def gt(c,p,alleles):
 if c not in handles:handles[c]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 f=handles[c];v=[r for r in f.fetch(c,p-1,p) if r.pos==p and len(r.alleles)==2 and set(r.alleles)==alleles]
 if len(v)!=1:return None,'no_unique_biallelic_match'
 arr=[]
 for s in samples:
  call=v[0].samples[s].get('GT')
  if call is None or len(call)!=2 or any(i is None for i in call):return None,'missing_calls'
  arr.append(sum(i==1 for i in call))
 arr=np.asarray(arr);return (None,'monomorphic_reference') if arr.std()==0 else (arr,v[0].id)
cache={};result=[]
for k,v in significant.items():
 m=lo.convert_coordinate('chr'+v['chrom'],v['pos38']-1);record={'reporter_variant':k,'source_genes':sorted(v['genes']),'adjusted_p_values':sorted(v['adjusted_p']),'reference_status':None,'candidates':[]}
 if len(m)!=1 or m[0][0]!='chr'+v['chrom'] or m[0][2]!='+':record['reference_status']='nonunique_or_reverse_liftover';result.append(record);continue
 yy,yid=gt(v['chrom'],int(m[0][1])+1,v['alleles']);record['reference_status']=yid
 if yy is not None:record['reference_alt_count']=int(yy.sum())
 for name in sorted(v['candidates']):
  t=candidates[name]
  if name not in cache:cache[name]=gt(t['chrom'],t['pos37'],t['alleles'])
  xx,xid=cache[name];x={'candidate':name,'distance_bp':v['pos38']-t['pos38'],'candidate_reference_status':xid}
  if xx is not None and yy is not None:x['r2_eur503']=float(np.corrcoef(xx,yy)[0,1]**2)
  record['candidates'].append(x)
 result.append(record)
receipt={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'panel_n':503,'window_bp':1000000,'candidate_count':len(candidates),'reporter_count':len(result),'results':result,'limitations':'Postselected significant reporter variants; frequency, LD, allele/coding, cell type and model differ. No participant outcomes, not causal or score validation. Small reference and remote VCF not full-checksummed.'}
(out/'stage6_salomon_significant_region_ld.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))
