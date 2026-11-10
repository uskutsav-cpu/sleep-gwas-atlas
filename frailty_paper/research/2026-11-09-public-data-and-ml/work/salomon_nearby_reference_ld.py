import csv,json,hashlib,collections
from pathlib import Path
import numpy as np,pysam
from pyliftover import LiftOver
root=Path('outputs'); source=root/'stage6_salomon_candidate_nearby_screen.json'; data=json.loads(source.read_text()); lo=LiftOver('work/hg38ToHg19.over.chain.gz')
cand={r['variant']:r for r in csv.DictReader((root/'integrated_15_candidate_review.csv').open())};cand.update({r['variant']:r for r in csv.DictReader((root/'expanded_22_candidate_risk_gene_table.csv').open()) if r['variant'] not in cand})
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur=sorted(x['sample'] for x in panel if x['super_pop']=='EUR');assert len(eur)==503
handles={}
def read(ch,pos,alleles):
 if ch not in handles:handles[ch]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{ch}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 f=handles[ch];rows=[v for v in f.fetch(ch,pos-1,pos) if v.pos==pos and len(v.alleles)==2 and set(v.alleles)==alleles]
 if len(rows)!=1:return None,'no_unique_biallelic_match'
 v=rows[0];alt=v.alleles[1];arr=[]
 for s in eur:
  gt=v.samples[s].get('GT')
  if gt is None or len(gt)!=2 or any(x is None for x in gt):return None,'missing_reference_call'
  arr.append(sum(x==1 for x in gt))
 if np.std(arr)==0:return None,'monomorphic_reference'
 return np.asarray(arr),v.id
results=[]
for name,near in data['matches'].items():
 r=cand[name];ch=r.get('GRCh37_chr',r.get('GRCh37_chr_assumed'));pos=int(r.get('GRCh37_pos',r.get('GRCh37_pos_assumed')));a={r['counted_allele'],r['other_allele']};x,xid=read(ch,pos,a)
 for v,meta in near.items():
  c,p,ref,alt=v.split('-');m=lo.convert_coordinate('chr'+c,int(p)-1)
  row={'candidate':name,'assay_variant_grch38':v,'distance_bp':meta['distance_bp'],'candidate_ref_record':xid,'status':None}
  if x is None:row['status']='candidate_'+xid
  elif len(m)!=1 or m[0][0]!='chr'+c or m[0][2]!='+':row['status']='nonunique_or_reverse_liftover'
  else:
   y,yid=read(c,int(m[0][1])+1,{ref,alt});row['reference_record']=yid
   if y is None:row['status']=yid
   else:row['status']='matched';row['r2_eur503']=float(np.corrcoef(x,y)[0,1]**2);row['assay_alt_count_1kg']=int(y.sum())
  results.append(row)
out={'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'candidate_inputs_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['integrated_15_candidate_review.csv','expanded_22_candidate_risk_gene_table.csv']},'chain_sha256':hashlib.sha256(Path('work/hg38ToHg19.over.chain.gz').read_bytes()).hexdigest(),'reference_panel_sha256':hashlib.sha256(Path('work/1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'sample_n':503,'results':results,'limitations':'Reference sample European, potentially related, not target donor LD; no confidence intervals, exact pair/positive liftover required, remote VCF not fully verified; reporter non-significance does not establish absence of molecular or disease effects.'}
(root/'stage6_salomon_nearby_reference_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(results,indent=2))
