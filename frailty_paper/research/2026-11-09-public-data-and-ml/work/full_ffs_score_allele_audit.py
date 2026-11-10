import csv,gzip,hashlib,json,pathlib,numpy as np
score=pathlib.Path('work/PGS005229_complete.txt.gz'); ffs=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv');N=2761868
md5='168e5c11469c0b27668a4c4d733276b1';sha='912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55'
def hashfile(f,algo):
 h=hashlib.new(algo)
 with f.open('rb') as z:
  for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert hashfile(score,'md5')==md5 and hashfile(ffs,'sha256')==sha
ids=np.zeros(N,dtype=np.uint64);chrs=np.zeros(N,dtype=np.uint8);pos=np.zeros(N,dtype=np.uint32);allele=np.zeros(N,dtype=np.uint8);state=np.zeros(N,dtype=np.uint8);samples={};n=0
with gzip.open(score,'rt',newline='') as f:
 for row in csv.DictReader((s for s in f if not s.startswith('#')),delimiter='\t'):
  i=n;n+=1;assert i<N
  id=row['rsID'];ch=row['chr_name'];bp=row['chr_position'];a=row['effect_allele']
  if id.startswith('rs') and id[2:].isdigit() and int(id[2:])>0:ids[i]=int(id[2:])
  else:state[i]=6;continue
  if not ch.isdigit() or not 1<=int(ch)<=22 or not bp.isdigit() or not 0<int(bp)<2**32 or len(a)!=1 or a not in 'ACGT':state[i]=6;continue
  chrs[i]=int(ch);pos[i]=int(bp);allele[i]=ord(a)
  if id in ('rs113851554','rs8067974','rs11152363','rs1942262','rs4886140','rs35881094','rs7734045','rs10865954'):samples[id]=i
assert n==N
order=np.argsort(ids,kind='stable');si=ids[order];dup=np.zeros(N,dtype=bool);dup[1:]|=si[1:]==si[:-1];dup[:-1]|=si[:-1]==si[1:];state[order[dup & (si!=0)]]=5
print('Loaded',n,'score records; invalid format',int(np.count_nonzero(state==6)),'duplicate-ID rows',int(np.count_nonzero(state==5)),flush=True)
seen=0;chunk=[]
def process(rows):
 global seen
 k=len(rows);fid=np.zeros(k,dtype=np.uint64);fc=np.zeros(k,dtype=np.uint8);fp=np.zeros(k,dtype=np.uint32);a1=np.zeros(k,dtype=np.uint8);a2=np.zeros(k,dtype=np.uint8);invalid=np.zeros(k,dtype=bool)
 for i,x in enumerate(rows):
  sn=x['SNP'];ch=x['CHROM'];p=x['POS'];a=x['A1'];b=x['A2']
  if sn.startswith('rs') and sn[2:].isdigit() and ch.isdigit() and p.isdigit() and 1<=int(ch)<=22:
   fid[i]=int(sn[2:]);fc[i]=int(ch);fp[i]=int(p)
   if len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b:a1[i]=ord(a);a2[i]=ord(b)
   else:invalid[i]=True
 ix=np.searchsorted(si,fid);hit=(ix<N)&(fid!=0)&(si[np.minimum(ix,N-1)]==fid);j=np.where(hit)[0];q=order[ix[j]];keep=(state[q]!=5)&(state[q]!=6);j=j[keep];q=q[keep]
 # Multiple input rows with the same ID/position are handled one-by-one to preserve ambiguity.
 for index,t in zip(j,q):
  t=int(t);v=int(state[t]);coordinate=(chrs[t]==fc[index] and pos[t]==fp[index])
  if not coordinate:
   if v==0:state[t]=4
   continue
  if v in (1,2,3,7,8):state[t]=7;continue
  if invalid[index]:state[t]=8
  elif allele[t]==a1[index]:state[t]=1
  elif allele[t]==a2[index]:state[t]=2
  else:state[t]=3
 seen+=k
with ffs.open('rt',newline='') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  chunk.append(row)
  if len(chunk)==50000:process(chunk);chunk=[]
 if chunk:process(chunk)
labels={0:'not_found_by_rsid',1:'exact_coord_effect_matches_FFS_A1',2:'exact_coord_effect_matches_FFS_A2',3:'exact_coord_effect_allele_not_in_FFS_pair',4:'rsid_found_no_coordinate_match',5:'score_duplicate_rsid_unresolved',6:'score_missing_coord_invalid_or_non_autosomal',7:'multiple_FFS_rows_same_rsid_coordinate',8:'FFS_exact_row_not_biallelic_SNV'}
counts={v:int(np.count_nonzero(state==k)) for k,v in labels.items()};assert sum(counts.values())==N
comparison=[]
existing={x['variant']:x for x in csv.DictReader(open('outputs/published_fi_ffs_pgs_candidate_lookup.csv')) if x['score_id']=='PGS005229' and x['matches']=='1'}
for id,i in samples.items():
 assert id in existing and existing[id]['coordinate_match']=='True'
 comparison.append({'variant':id,'status':labels[int(state[i])]})
assert len(comparison)==8
summary={'purpose':'Full PGS005229 counted-allele and coordinate comparison with Ye FFS source; no source modification, score application or independent validation','score_rows':N,'FFS_rows':seen,'score_md5':md5,'FFS_source_sha256':sha,'status_counts':counts,'eight_existing_candidate_checks':comparison,'limitations':'FFS A1/BETA convention is a documented format inference, source uses non-SNV variants; scores with missing coordinate and duplicate IDs remain unresolved. Exact rsID/coordinate and allele match does not establish weights, phase, LD, sample independence or predictive validity.'}
pathlib.Path('outputs/ffs_pgs_full_counted_allele_audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
