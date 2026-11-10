import csv,gzip,hashlib,json,pathlib,collections,numpy as np
score=pathlib.Path('work/PGS005228_complete.txt.gz');fi=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz'); N=2379429
expected_md5='584c4839ba4d8caa299ff558ed161259';expected_fi='241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de'
def hfile(f,algo):
 h=hashlib.new(algo)
 with f.open('rb') as z:
  for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert hfile(score,'md5')==expected_md5 and hfile(fi,'sha256')==expected_fi
weights=np.zeros(N,dtype=np.float64); fi_aligned=np.full(N,np.nan,dtype=np.float64); fi_p=np.full(N,np.nan,dtype=np.float64)
ids=np.zeros(N,dtype=np.uint64); chrs=np.zeros(N,dtype=np.uint8); pos=np.zeros(N,dtype=np.uint32); allele=np.zeros(N,dtype=np.uint8);state=np.zeros(N,dtype=np.uint8)
count=0;non_numeric=0;invalid=0; chr6_mhc=np.zeros(N,dtype=bool); sample_indices={}
with gzip.open(score,'rt',newline='') as f:
 for row in csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t'):
  i=count;count+=1;assert i<N
  weights[i]=float(row['effect_weight'])
  name=row['rsID'];ch=row['chr_name'];bp=row['chr_position'];a=row['effect_allele']
  key=name+'|'+ch+'|'+bp
  if int.from_bytes(hashlib.sha256(key.encode()).digest()[:2],'big')<64:sample_indices[name]=i
  if not name.startswith('rs') or not name[2:].isdigit() or int(name[2:])==0:
   state[i]=6;non_numeric+=1;continue
  ids[i]=int(name[2:])
  if not ch.isdigit() or not (1<=int(ch)<=22) or not bp.isdigit() or not 0<int(bp)<2**32 or len(a)!=1 or a not in 'ACGT':
   state[i]=6;invalid+=1;continue
  chrs[i]=int(ch);pos[i]=int(bp);allele[i]=ord(a);chr6_mhc[i]=(ch=='6' and 25000000<=int(bp)<=34000000)
assert count==N
order=np.argsort(ids,kind='stable');sortids=ids[order];duplicate=np.zeros(N,dtype=bool);duplicate[1:]|=sortids[1:]==sortids[:-1];duplicate[:-1]|=sortids[:-1]==sortids[1:];state[order[duplicate & (sortids!=0)]]=5
print('Loaded',count,'score rows; nonnumeric IDs',non_numeric,'invalid format',invalid,'duplicate-ID rows',int(np.count_nonzero(state==5)),flush=True)
seen=0; chunk=[]
def process(chunk):
 global seen
 n=len(chunk);fid=np.zeros(n,dtype=np.uint64);fchr=np.zeros(n,dtype=np.uint8);fpos=np.zeros(n,dtype=np.uint32);a1=np.zeros(n,dtype=np.uint8);a2=np.zeros(n,dtype=np.uint8); beta=np.full(n,np.nan,dtype=np.float64); pv=np.full(n,np.nan,dtype=np.float64)
 for i,x in enumerate(chunk):
  sn=x['SNP'];ch=x['CHR'];bp=x['BP'];aa=x['A1'];bb=x['A2']
  if sn.startswith('rs') and sn[2:].isdigit() and ch.isdigit() and bp.isdigit() and len(aa)==len(bb)==1:
   fid[i]=int(sn[2:]);fchr[i]=int(ch);fpos[i]=int(bp);a1[i]=ord(aa);a2[i]=ord(bb);beta[i]=float(x['BETA']);pv[i]=float(x['P'])
 ix=np.searchsorted(sortids,fid);inrange=ix<N;safe=np.minimum(ix,N-1); hit=inrange & (fid!=0) & (sortids[safe]==fid);findex=np.where(hit)[0]; si=order[ix[findex]]
 eligible=state[si]!=5;findex=findex[eligible];si=si[eligible];eligible=state[si]!=6;findex=findex[eligible];si=si[eligible]
 if len(si)!=len(np.unique(si)):raise RuntimeError('FI source duplicated rsID within a chunk')
 if np.any(state[si]!=0):raise RuntimeError('FI source duplicated rsID across chunks')
 coord=(chrs[si]==fchr[findex]) & (pos[si]==fpos[findex]);state[si]=4
 q=si[coord];j=findex[coord];state[q]=3
 direct=allele[q]==a1[j];swap=allele[q]==a2[j];state[q[direct]]=1;state[q[swap]]=2
 fi_aligned[q[direct]]=beta[j[direct]];fi_aligned[q[swap]]=-beta[j[swap]];fi_p[q]=pv[j]
 seen+=len(chunk)
with gzip.open(fi,'rt',newline='') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  chunk.append(row)
  if len(chunk)==50000:process(chunk);chunk=[]
 if chunk:process(chunk)
labels={0:'not_found_by_rsid_in_filtered_FI',1:'coordinate_match_effect_allele_matches_FI_A1',2:'coordinate_match_effect_allele_matches_FI_A2',3:'coordinate_match_effect_allele_not_in_FI_pair',4:'rsid_match_coordinate_mismatch',5:'score_duplicate_rsid_unresolved',6:'invalid_or_non_autosomal_score_record'}
counts={labels[k]:int(np.count_nonzero(state==k)) for k in labels};assert sum(counts.values())==N
samples=list(csv.DictReader(open('outputs/fi_pgs_allele_sample_audit.csv')));translated={0:'not_found_in_filtered_FI',1:'exact_coordinate_effect_equals_FI_A1',2:'exact_coordinate_effect_equals_FI_A2',3:'exact_coordinate_effect_not_in_FI_pair',4:'rsid_found_coordinate_mismatch',5:'multiple_rsid_matches_unresolved'}
compared=0; errors=[]
for x in samples:
 idx=sample_indices.get(x['PGS_variant']);assert idx is not None
 if int(state[idx]) in translated:
  compared+=1
  if translated[int(state[idx])]!=x['status']:errors.append(x['PGS_variant'])
assert not errors and compared==len(samples)
summary={'purpose':'Full PGS005228 counted-allele consistency check against filtered harmonized Atkins FI GWAS; no PGS modification, no prediction or independent validation','score_rows':N,'FI_rows':seen,'score_md5':expected_md5,'FI_sha256':expected_fi,'status_counts':counts,'score_chr6_MHC_rows':int(chr6_mhc.sum()),'unmatched_chr6_MHC':int(np.count_nonzero(chr6_mhc & (state==0))),'independent_previous_sample_rows_agree':compared,'limitations':'Match uses rsID, hg19 position and allele strings only; FFS not evaluated; genomic variant representation, effect weights, all sources, target genotypes and predictive utility remain unresolved. Unmatched is not a failed risk or zero effect.'}
matched=(state==1)|(state==2)
strata={'all_matched':matched,'source_p_below_5e-8':matched & (fi_p<5e-8),'source_p_5e-8_to_0.05':matched & (fi_p>=5e-8) & (fi_p<0.05),'source_p_at_least_0.05':matched & (fi_p>=0.05)}
stats={}
for label,m in strata.items():
 w=weights[m];b=fi_aligned[m]; valid=np.isfinite(w)&np.isfinite(b);w=w[valid];b=b[valid];nz=(w!=0)&(b!=0)
 stats[label]={'rows':int(m.sum()),'finite_pairs':len(w),'zero_score_weights':int(np.count_nonzero(w==0)),'zero_source_betas':int(np.count_nonzero(b==0)),'same_sign_nonzero':int(np.count_nonzero(w[nz]*b[nz]>0)),'opposite_sign_nonzero':int(np.count_nonzero(w[nz]*b[nz]<0)),'pearson_weight_vs_aligned_marginal_beta':float(np.corrcoef(w,b)[0,1]) if len(w)>1 and np.std(w)>0 and np.std(b)>0 else None}
sampleidx=np.array([sample_indices[x['PGS_variant']] for x in samples if x['status']=='exact_coordinate_effect_equals_FI_A2'],dtype=np.int64)
sm=sampleidx[(fi_p[sampleidx]<1)&(weights[sampleidx]!=0)&(fi_aligned[sampleidx]!=0)]
assert len(sm)==2286 and np.count_nonzero(weights[sm]*fi_aligned[sm]>0)==2104
summary['weight_orientation_audit']={'strata':stats,'sample_check':{'nonzero_pairs_p_lt_1':len(sm),'same_sign':int(np.count_nonzero(weights[sm]*fi_aligned[sm]>0))},'interpretation':'Descriptive within same underlying FI source; LD and SBayesR shrinkage can change signs and magnitudes; no independent predictive validation or accuracy claim.'}
pathlib.Path('outputs/fi_pgs_weight_orientation_audit.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
