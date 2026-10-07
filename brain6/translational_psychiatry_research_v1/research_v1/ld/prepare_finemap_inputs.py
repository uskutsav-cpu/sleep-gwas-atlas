#!/usr/bin/env python3
"""Fixed gates for source-qualified chr5C fine mapping; raw inputs remain local."""
import csv,gzip,hashlib,json,math,pathlib
from collections import Counter,defaultdict
import numpy as np
from scipy import stats
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[3];DISC=HERE.parent/'discovery';RAW=ROOT/'work'/'ld_genotypes_research_v1';OUT=RAW/'C_finemap_inputs'
CASE_FRAC=109402/386533

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')
def write(p,rows,fields=None):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def orientation(a,b,ref,alt):
 if {a,b} in [{'A','T'},{'C','G'}]:return None
 if (a,b)==(alt,ref):return ('EXACT',1)
 if (a,b)==(ref,alt):return ('SWAP',-1)
 comp=str.maketrans('ACGT','TGCA')
 if (a.translate(comp),b.translate(comp))==(alt,ref):return ('COMPLEMENT',1)
 if (a.translate(comp),b.translate(comp))==(ref,alt):return ('COMPLEMENT_SWAP',-1)
 return None

def main():
 protocol=HERE/'FINEMAP_PRE_FIT_PROTOCOL.md';expected=(HERE/'FINEMAP_PRE_FIT_PROTOCOL.sha256').read_text().split()[0]
 if sha(protocol)!=expected:raise ValueError('Fine-map protocol changed after freeze')
 receipts={r['trait']:r for r in json.loads((DISC/'acquisition_receipts.json').read_text())}
 variants=list(csv.DictReader((HERE/'chr5_variants.tsv').open(),delimiter='\t'));bypos=defaultdict(list)
 for v in variants:bypos[int(v['pos'])].append(v)
 OUT.mkdir(exist_ok=True)
 counts=[];trait_data={};exclusions=[]
 for trait in ['insomnia','adhd']:
  p=DISC/f'{trait}_C.raw.tsv.gz'
  if sha(p)!=receipts[trait]['regional_sha256']['C']:raise ValueError('Regional input checksum mismatch')
  cnt=Counter();valid={};zs_from_p_deltas=[];min_source_p=1
  for x in csv.DictReader(gzip.open(p,'rt'),delimiter='\t'):
   cnt['source_rows']+=1;reason='';bp=int(x['BP']);a,b=x['A1'],x['A2'];o=None;v=None;min_source_p=min(min_source_p,float(x['P']))
   if a not in 'ACGT' or len(a)!=1 or b not in 'ACGT' or len(b)!=1 or a==b:reason='INVALID_OR_NONBIALLELIC_ALLELES'
   elif {a,b} in [{'A','T'},{'C','G'}]:reason='PALINDROMIC'
   else:
    try:
     beta=math.log(float(x['OR']));se=float(x['SE']);pval=float(x['P']);info=float(x['INFO'])
     if trait=='insomnia':n=float(x['N']);neff=n*4*CASE_FRAC*(1-CASE_FRAC);freq=float(x['MAF'])
     else:nca,nco=float(x['Nca']),float(x['Nco']);n=nca+nco;neff=4/(1/nca+1/nco);freq=float(x['FRQ_U_186843'])
     if not all(math.isfinite(t) for t in [beta,se,pval,info,n,neff,freq]) or se<=0 or n<=0 or neff<=0 or not(0<pval<=1) or not(0<freq<1):reason='INVALID_NUMERIC_FIELDS'
     elif info<.9:reason='INFO_LT_0.9'
     elif min(freq,1-freq)<.01:reason='GWAS_MAF_LT_0.01'
    except (ValueError,ZeroDivisionError):reason='INVALID_NUMERIC_FIELDS'
    if not reason:
     cnt['valid_source_nonpal_QC_rows']+=1
     matches=[(vv,oo) for vv in bypos[bp] if (oo:=orientation(a,b,vv['REF'],vv['ALT']))]
     if len(matches)>1:raise ValueError('Ambiguous variant reference identity')
     if not matches:reason='REFERENCE_COORDINATE_ALLELE_ABSENT'
     else:
      v,o=matches[0];altfreq=freq if o[1]==1 else 1-freq
      if abs(altfreq-float(v['EUR_ALT_AF']))>.15:reason='SOURCE_REFERENCE_AF_DIFF_GT_0.15'
   if reason:cnt[reason]+=1;exclusions.append({'trait':trait,'source_rsid':x['SNP'],'chr':x['CHR'],'pos':bp,'A1':a,'A2':b,'reason':reason});continue
   if v['variant_key'] in valid:raise ValueError('Duplicate harmonized GWAS key')
   valid[v['variant_key']]={'variant_key':v['variant_key'],'SNP':x['SNP'],'CHR':5,'BP':bp,'ALT':v['ALT'],'REF':v['REF'],'beta_ALT':beta*o[1],'SE':se,'Z_ALT':beta*o[1]/se,'P_source':pval,'N_analyzed':n,'N_eff':neff,'ALT_AF_source':altfreq,'ALT_AF_reference':float(v['EUR_ALT_AF']),'INFO':info,'orientation':o[0],'reference_index':int(v['variant_order'])}
   zs_from_p_deltas.append(abs(abs(beta/se)-float(stats.norm.isf(pval/2))))
  median=float(np.median([v['N_eff'] for v in valid.values()]));consistent={k:v for k,v in valid.items() if .8*median<=v['N_eff']<=1.2*median}
  for k,v in valid.items():
   if k not in consistent:exclusions.append({'trait':trait,'source_rsid':v['SNP'],'chr':5,'pos':v['BP'],'A1':'','A2':'','reason':'N_EFF_OUTSIDE_20PCT_MEDIAN'});cnt['N_EFF_OUTSIDE_20PCT_MEDIAN']+=1
  gates={'500_oriented_snps':len(consistent)>=500,'80pct_reference_coverage':len(consistent)/len(variants)>=.8,'80pct_valid_source_QC_coverage':len(consistent)/cnt['valid_source_nonpal_QC_rows']>=.8,'80pct_N_consistency':len(consistent)/len(valid)>=.8,'source_full_window_min_P_lt5e8':min_source_p<5e-8,'admitted_min_P_lt5e8':min(v['P_source'] for v in consistent.values())<5e-8}
  cnt['oriented_before_N_gate']=len(valid);cnt['after_N_gate']=len(consistent)
  summary={'trait':trait,'region':'chr5:103447968-104447968','source_sha256':receipts[trait]['source_sha256'],'regional_raw_sha256':sha(p),'source_rows':cnt['source_rows'],'reference_eligible':len(variants),'valid_source_QC_rows':cnt['valid_source_nonpal_QC_rows'],'oriented_before_N_gate':len(valid),'after_N_gate':len(consistent),'reference_coverage':len(consistent)/len(variants),'source_QC_coverage':len(consistent)/cnt['valid_source_nonpal_QC_rows'],'N_consistency_fraction':len(consistent)/len(valid),'N_eff_median':median,'N_eff_min':min(v['N_eff'] for v in consistent.values()),'N_eff_max':max(v['N_eff'] for v in consistent.values()),'min_P_source_full_window':min_source_p,'min_P_admitted':min(v['P_source'] for v in consistent.values()),'max_abs_Z_OR_SE_vs_sourceP_Z_discrepancy':max(zs_from_p_deltas),'median_abs_Z_OR_SE_vs_sourceP_Z_discrepancy':float(np.median(zs_from_p_deltas)),'input_gate_status':'PASS' if all(gates.values()) else 'HOLD_INPUT_QC','failed_gates':';'.join(k for k,v in gates.items() if not v),'source_N_semantics':'LITERAL_ANALYZED_N_X_GLOBAL_CASE_FRACTION_NEFF_APPROXIMATION' if trait=='insomnia' else 'PER_VARIANT_NCA_NCO_EFFECTIVE_N','source_N_approximation_limits':'Insomnia perSNP casefraction not supplied; published UKB casefraction109402/386533 assumed; binaryRSS is approximation. ADHD metaN pervariant may represent heterogeneous cohort weights.','protocol_sha256':expected}
  counts.append(summary);trait_data[trait]=consistent
  dump(HERE/f'{trait}_C_input_diagnostics.json',{'summary':summary,'counts':dict(cnt),'gates':gates})
 shared=sorted(set(trait_data['insomnia'])&set(trait_data['adhd']),key=lambda k:trait_data['insomnia'][k]['reference_index'])
 sharedgate=len(shared)>=500 and all(r['input_gate_status']=='PASS' for r in counts)
 write(HERE/'finemap_input_validation.tsv',counts)
 write(HERE/'finemap_input_exclusions.tsv',exclusions,['trait','source_rsid','chr','pos','A1','A2','reason'])
 manifest={'protocol_sha256':expected,'shared_snps':len(shared),'status':'READY_FOR_RSS_CONSISTENCY_DIAGNOSTICS' if sharedgate else 'HOLD_INPUT_QC','trait_checks':counts,'sources_same_discovery_as_historical_ABF':True,'independent_trait_replication':False,'local_only_input_files':{},'LD_full_sha256':sha(RAW/'chr5_signed_LD.npz')}
 if sharedgate:
  for t in ['insomnia','adhd']:
   rows=[trait_data[t][k] for k in shared];write(OUT/f'{t}_C_shared.tsv',rows);manifest['local_only_input_files'][t]={'path':str(OUT/f'{t}_C_shared.tsv'),'sha256':sha(OUT/f'{t}_C_shared.tsv'),'median_n_eff_shared':float(np.median([v['N_eff'] for v in rows]))}
  ix=[trait_data['insomnia'][k]['reference_index'] for k in shared];full=np.load(RAW/'chr5_signed_LD.npz')['R'];mat=full[np.ix_(ix,ix)];(OUT/'shared_signed_LD.float64.bin').write_bytes(mat.astype('<f8').tobytes(order='F'))
  manifest['local_only_input_files']['LD']={'path':str(OUT/'shared_signed_LD.float64.bin'),'sha256':sha(OUT/'shared_signed_LD.float64.bin'),'n_snps':len(shared),'byte_order':'little','matrix_order':'column_major'}
 dump(HERE/'finemap_input_manifest.json',manifest)
 print(json.dumps({'status':manifest['status'],'shared_snps':len(shared),'trait_diagnostics':counts}),flush=True)
if __name__=='__main__':main()
