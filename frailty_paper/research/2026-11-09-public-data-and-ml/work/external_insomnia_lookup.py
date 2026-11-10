import csv,hashlib,json,pathlib,collections
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered')
inputs={'R13':(base/'R13_INSOMNIA_HM3.tsv','47ab07104aa27fd641711d7df04d64b9c90a3bf6986af1695c3f916fabd9f3c9'),'MVP':(base/'MVP_INSOMNIA_HM3.tsv','1a96b8b62d79e5d9b0854a4ea78dac7f06cdd3ba15e7ed0e9c4f4acac16c8a60')}
candidates=list(csv.DictReader(open('outputs/published_candidate_marginal_lookup.csv')));ids={x['candidate'] for x in candidates}; found={}; receipt={}
for name,(path,expected) in inputs.items():
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 assert h.hexdigest()==expected,(name,h.hexdigest())
 acc={};n=0
 with path.open(newline='') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   n+=1
   if x['SNP'] in ids:acc.setdefault(x['SNP'],[]).append(x)
 found[name]=acc;receipt[name]={'sha256':expected,'rows':n,'candidate_ids':len(acc),'duplicate_ids':sum(len(v)>1 for v in acc.values())}
comp={'A':'T','T':'A','C':'G','G':'C'}
rows=[]
for x in candidates:
 out={k:x[k] for k in ('candidate','song_rank_cpassoc','song_chr','song_pos_reported','insomnia_A1','insomnia_A2','insomnia_BETA','insomnia_P','FI_BETA','FI_P','FFS_BETA','FFS_P')}
 for n in inputs:
  ls=found[n].get(x['candidate'],[]);out[n+'_matches']=len(ls)
  for col in ('CHR','A1','A2','BETA','SE','P','N','EAF'):out[n+'_'+col]=ls[0].get(col,'') if len(ls)==1 else ''
  out[n+'_relation']='not_assessed';out[n+'_beta_aligned_to_Jansen_A1']='';out[n+'_sign_concordance']='not_assessed'
  if len(ls)!=1 or not x['insomnia_A1']:continue
  v=ls[0]; a=(x['insomnia_A1'],x['insomnia_A2']); b=(v['A1'],v['A2']);out[n+'_relation']='incompatible'
  if str(v['CHR'])!=str(x['insomnia_CHR']) or set(a) in ({'A','T'},{'C','G'}) or any(t not in comp for t in a+b):continue
  bc=(comp[b[0]],comp[b[1]])
  mode='direct' if a==b else 'swap' if a==b[::-1] else 'complement' if a==bc else 'complement_swap' if a==bc[::-1] else 'incompatible';out[n+'_relation']=mode
  if mode=='incompatible':continue
  val=float(v['BETA'])*(-1 if mode in ('swap','complement_swap') else 1);out[n+'_beta_aligned_to_Jansen_A1']=val
  product=val*float(x['insomnia_BETA']);out[n+'_sign_concordance']='same' if product>0 else 'opposite' if product<0 else 'zero'
 out['status']='Post-selection descriptive sleep-source lookup only; unknown participant overlap; rsID/chr/allele match without coordinates; not independent frailty validation or shared-risk score'
 rows.append(out)
with open('outputs/published_candidates_external_insomnia_lookup.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
pathlib.Path('work/external_lookup_receipt.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
for name in inputs:
 counts=collections.Counter(x[name+'_sign_concordance'] for x in rows)
 top=[x for x in rows if int(x['song_rank_cpassoc'])<=30 and x['insomnia_A1']]
 print(name,'direction count',dict(counts),'top prototype',dict(collections.Counter(x[name+'_sign_concordance'] for x in top)))
