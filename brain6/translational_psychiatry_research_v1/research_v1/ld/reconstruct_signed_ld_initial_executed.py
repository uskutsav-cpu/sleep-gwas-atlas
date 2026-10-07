#!/usr/bin/env python3
"""Actual independent-reference LD from exact indexed EBI VCF range bytes."""
import csv,gzip,hashlib,json,math,pathlib,re,struct,sys
from collections import Counter
import numpy as np
from scipy import stats
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[3]
RAW=ROOT/'work'/'ld_genotypes_research_v1'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for d in iter(lambda:f.read(1024*1024),b''):h.update(d)
 return h.hexdigest()
def dump(p,o):p.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')
def write(p,rows,fields=None):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def complete_bgzf(data):
 """Yield only complete BGZF members, validating their CRC via gzip."""
 offset=0
 while offset+18<=len(data):
  if data[offset:offset+4]!=b'\x1f\x8b\x08\x04':raise ValueError('Invalid BGZF header')
  xlen=struct.unpack_from('<H',data,offset+10)[0];extra=data[offset+12:offset+12+xlen];i=0;size=None
  while i+4<=len(extra):
   tag=extra[i:i+2];n=struct.unpack_from('<H',extra,i+2)[0]
   if tag==b'BC' and n==2:size=struct.unpack_from('<H',extra,i+4)[0]+1
   i+=4+n
  if size is None:raise ValueError('Missing BGZF BC subfield')
  if offset+size>len(data):break
  yield gzip.decompress(data[offset:offset+size]);offset+=size

def header_and_samples(path):
 text=b''.join(complete_bgzf(path.read_bytes())).decode()
 lines=text.splitlines();head=[x for x in lines if x.startswith('#')]
 cols=next(x for x in head if x.startswith('#CHROM')).split('\t')
 if len(cols[9:])!=len(set(cols[9:])):raise ValueError('Duplicate sample IDs')
 if not any('assembly=b37' in x for x in head):raise ValueError('Not verified GRCh37')
 return '\n'.join(head)+'\n',cols[9:]

def regional_rows(paths,chrom,start,stop):
 for p in paths:
  text=b''.join(complete_bgzf(p.read_bytes())).decode()
  lines=text.split('\n')[1:-1] # discard possible leading/trailing fragments
  for l in lines:
   x=l.split('\t')
   if len(x)<9:continue
   if x[0]!=str(chrom):continue
   pos=int(x[1])
   if start<=pos<=stop:yield x

def scalar_pearson(a,b):
 a=[float(v) for v in a];b=[float(v) for v in b];ma=sum(a)/len(a);mb=sum(b)/len(b)
 num=sum((x-ma)*(y-mb) for x,y in zip(a,b));den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
 return num/den

def genotype_reconstruction(obj,panel):
 c=obj['chr'];headerpath=next(RAW.glob(f'chr{c}_header_*.bgzf'));head,samples=header_and_samples(headerpath)
 if set(samples)!=set(panel):raise ValueError('VCF/panel sample mismatch')
 eur=[s for s in samples if panel[s]['super_pop']=='EUR'];indices=[samples.index(s)+9 for s in eur]
 (HERE/f'chr{c}_vcf_header.txt').write_text(head)
 write(HERE/f'chr{c}_sample_manifest.tsv',[{'order':i,'sample':s,'pop':panel[s]['pop'],'super_pop':panel[s]['super_pop']} for i,s in enumerate(eur)])
 paths=sorted(RAW.glob(f'chr{c}_genotypes_*.bgzf'),key=lambda p:int(p.name.split('_')[2]))
 counts=Counter();admitted=[];excluded=[];dosages=[];seen=set();rs_seen=set();maplead={}
 if c==5:
  for rs in ['rs2431108','rs77960']:
   data=json.loads((RAW/(rs+'_grch37.json')).read_text());m=[m for m in data['mappings'] if m['assembly_name']=='GRCh37' and m['seq_region_name']=='5' and m['strand']==1]
   if len(m)!=1:raise ValueError('Ambiguous lead mapping')
   maplead[m[0]['start']]=(rs,set(m[0]['allele_string'].split('/')))
 for row in regional_rows(paths,c,obj['start'],obj['stop']):
  chrom,pos,rs,ref,alt,qual,filt,info,fmt=row[:9];pos=int(pos);key=f'{c}:{pos}:{ref}:{alt}'
  if key in seen:counts['duplicate_retrieval_record']+=1;continue
  seen.add(key);counts['regional_unique_records']+=1;reason='';gt=None
  if ref not in 'ACGT' or len(ref)!=1 or alt not in 'ACGT' or len(alt)!=1 or ref==alt:reason='NOT_BIALLELIC_ACGT_SNP'
  elif filt!='PASS':reason='NON_PASS_FILTER'
  elif {ref,alt} in [{'A','T'},{'C','G'}]:reason='PALINDROMIC_EXCLUDED'
  elif fmt!='GT':reason='NON_GT_ONLY_FORMAT'
  else:
   tokens=[row[i] for i in indices]
   if not all(re.fullmatch(r'[01][|/][01]',x) for x in tokens):reason='NONCOMPLETE_DIPLOID_BIALLELIC_GT'
   else:
    gt=np.array([int(x[0])+int(x[2]) for x in tokens],dtype=np.int8);freq=float(gt.mean()/2)
    if min(freq,1-freq)<.01:reason='EUR_MAF_LT_0.01'
  if rs!='.' and rs in rs_seen:raise ValueError('Duplicate source rsID')
  if rs!='.':rs_seen.add(rs)
  mapped=''
  if pos in maplead and {ref,alt}.issubset(maplead[pos][1]):mapped=maplead[pos][0]
  if reason:excluded.append({'variant_key':key,'source_id':rs,'mapped_lead_rsid':mapped,'reason':reason});counts[reason]+=1;continue
  counts['admitted']+=1
  details={'variant_order':0,'variant_key':key,'chr':c,'pos':pos,'REF':ref,'ALT':alt,'source_id':rs,'mapped_lead_rsid':mapped,'EUR_ALT_AF':freq,'EUR_MAF':min(freq,1-freq),'n_reference':len(eur),'phased_GT':all('|' in x for x in tokens)}
  # INFO frequency is rounded and used only as an independent dosage parsing check.
  info_map=dict(x.split('=',1) for x in info.split(';') if '=' in x)
  details['INFO_EUR_AF']=info_map.get('EUR_AF','')
  details['INFO_EUR_AF_abs_error']=abs(freq-float(info_map['EUR_AF'])) if 'EUR_AF' in info_map else ''
  admitted.append(details);dosages.append(gt)
 order=sorted(range(len(admitted)),key=lambda i:(admitted[i]['pos'],admitted[i]['REF'],admitted[i]['ALT'],admitted[i]['source_id']))
 admitted=[admitted[i] for i in order];G=np.array([dosages[i] for i in order],dtype=np.int8).T
 for i,r in enumerate(admitted):r['variant_order']=i
 write(HERE/f'chr{c}_variants.tsv',admitted);write(HERE/f'chr{c}_excluded_variants.tsv',excluded,['variant_key','source_id','mapped_lead_rsid','reason'])
 np.savez_compressed(HERE/f'chr{c}_genotype_factors.npz',ALT_dosage=G,samples=np.array(eur),variant_keys=np.array([r['variant_key'] for r in admitted]),population=np.array([panel[s]['pop'] for s in eur]))
 centered=G.astype(np.float64)-G.mean(axis=0);sd=np.sqrt(np.sum(centered**2,axis=0));X=centered/sd
 n,p=G.shape
 if p*p*8>512*2**20:raise RuntimeError('512MiB matrix budget exceeded')
 R=X.T@X
 nonzero_spectrum=np.linalg.eigvalsh(X@X.T)
 maxe=float(nonzero_spectrum[-1]);mine=float(min(0,nonzero_spectrum[0])) if p>n else float(np.linalg.eigvalsh(R)[0]);rank=int(sum(nonzero_spectrum>1e-8*maxe))
 rng=np.random.default_rng(20261007);pairs=rng.integers(0,p,size=(100,2));inderr=max(abs(float(R[a,b])-scalar_pearson(G[:,a],G[:,b])) for a,b in pairs)
 scipyerr=max(abs(float(R[a,b])-float(stats.pearsonr(G[:,a],G[:,b]).statistic)) for a,b in pairs)
 diag=float(np.max(np.abs(np.diag(R)-1)));asym=float(np.max(np.abs(R-R.T)));maxr=float(np.max(np.abs(R)));finite=bool(np.isfinite(R).all())
 passed=finite and diag<=1e-12 and asym<=1e-12 and maxr<=1+1e-12 and mine>=-1e-8*maxe and inderr<=1e-12
 # Save exact float64 signed LD; no rounding, clipping or projection.
 np.savez_compressed(HERE/f'chr{c}_signed_LD.npz',R=R,variant_keys=np.array([r['variant_key'] for r in admitted]),effect_allele=np.array([r['ALT'] for r in admitted]))
 summary={'region':f'chr{c}:{obj["start"]}-{obj["stop"]}','n_reference':n,'n_variants':p,'matrix_status':'PASS_REFERENCE_NUMERIC_QC' if passed else 'FAIL_REFERENCE_NUMERIC_QC','min_Gram_eigenvalue':mine,'max_eigenvalue':maxe,'effective_rank_relative_1e-8':rank,'diagonal_max_error':diag,'asymmetry_max_error':asym,'max_abs_signed_r':maxr,'nonfinite_entries':int(np.size(R)-np.isfinite(R).sum()),'independent_scalar_max_error_100_pairs':inderr,'scipy_pearson_max_error_100_pairs':scipyerr,'max_INFO_AF_abs_error':max(float(r['INFO_EUR_AF_abs_error']) for r in admitted if r['INFO_EUR_AF_abs_error']!=''),'genotype_status':'SOURCE_VERIFIED_503_EUR_REFERENCE','GWAS_alignment_status':'NOT_YET_VERIFIED','fine_mapping_admission':'BLOCKED_DENSE_GWAS_CONTRACTS_AND_SUMMARY_LD_CONSISTENCY','source_index_sha256':obj['index_sha256'],'protocol_sha256':sha(HERE/'LD_PRE_OUTCOME_PROTOCOL.md'),'genotype_factor_sha256':sha(HERE/f'chr{c}_genotype_factors.npz'),'matrix_sha256':sha(HERE/f'chr{c}_signed_LD.npz')}
 dump(HERE/f'chr{c}_reconstruction_diagnostics.json',{'summary':summary,'variant_counts':dict(counts),'matrix_allele_orientation':'ALT','no_LD_values_modified':True,'rank_deficiency_expected_p_gt_n':p>n,'reference_source':'EBI_PHASE3_V5B_20130502_GRCH37','independent_reference_not_independent_GWAS_replication':True})
 if c==5:
  lead_index={r['mapped_lead_rsid']:i for i,r in enumerate(admitted) if r['mapped_lead_rsid']}
  if len(lead_index)!=2:raise RuntimeError('Two exact-allele chr5 lead variants not both admitted: '+str(lead_index))
  a,b=lead_index['rs2431108'],lead_index['rs77960'];res=[]
  for pop in ['EUR','CEU','FIN','GBR','IBS','TSI']:
   ix=[i for i,s in enumerate(eur) if pop=='EUR' or panel[s]['pop']==pop];ga=G[ix,a];gb=G[ix,b];r=scalar_pearson(ga,gb)
   res.append({'population':pop,'n_reference':len(ix),'variant1':'rs2431108','variant1_key':admitted[a]['variant_key'],'variant1_counted_allele':admitted[a]['ALT'],'variant1_ALT_AF':float(ga.mean()/2),'variant2':'rs77960','variant2_key':admitted[b]['variant_key'],'variant2_counted_allele':admitted[b]['ALT'],'variant2_ALT_AF':float(gb.mean()/2),'signed_r':r,'r_squared':r*r,'method':'PEARSON_ALT_DOSAGE','scope':'REFERENCE_LD_SENSITIVITY_NOT_CONDITIONAL_GWAS'})
  write(HERE/'chr5_prior_vs_brain6_lead_LD.tsv',res)
 print(json.dumps(summary),flush=True);return summary

def main():
 plan=json.loads((HERE/'genotype_acquisition_manifest.json').read_text())
 if plan['protocol_sha256']!=sha(HERE/'LD_PRE_OUTCOME_PROTOCOL.md'):raise RuntimeError('Protocol hash mismatch')
 receipts=json.loads((HERE/'genotype_acquisition_receipts.json').read_text())
 for r in receipts:
  if sha(pathlib.Path(r['path']))!=r['sha256']:raise RuntimeError('Raw checksum mismatch')
 panel={r['sample']:r for r in csv.DictReader((RAW/'integrated_call_samples_v3.20130502.ALL.panel').open(),delimiter='\t')}
 summaries=[genotype_reconstruction(o,panel) for o in plan['regions']]
 write(HERE.parent/'LD_VALIDATION.tsv',summaries)
 dump(HERE/'reconstruction_provenance.json',{'script_sha256':sha(pathlib.Path(__file__)),'acquisition_manifest_sha256':sha(HERE/'genotype_acquisition_manifest.json'),'acquisition_receipts_sha256':sha(HERE/'genotype_acquisition_receipts.json'),'metadata_receipts_sha256':sha(HERE/'metadata_acquisition_receipts.json'),'sample_pedigree_audit_sha256':sha(HERE/'sample_pedigree_audit.json'),'numpy':np.__version__,'python':sys.version,'results':summaries})
if __name__=='__main__':main()
