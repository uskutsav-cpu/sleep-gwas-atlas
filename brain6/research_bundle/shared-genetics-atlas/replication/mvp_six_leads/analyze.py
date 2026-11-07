import csv,gzip,hashlib,json,math,pathlib,collections
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];six=root/'replication/six_lead_joint/six_leads.tsv';previous=root/'replication/six_lead_joint/receipt.json';brain=root/'evidence/BRAIN6_FINAL_CANDIDATES.tsv';prov=root.parent/'brain6-alphagenome-pilot/results/provenance.json';src=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');mvp=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90475826.tsv.gz');mvp_receipt=root/'placo_diagnostic/cross_source_followup/receipt.json'
def read(q):
 with q.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
expect=json.loads(mvp_receipt.read_text())['source_digests']['MVP_sha256'];actual=sha(mvp);assert actual==expect
rows=read(six);cands={r['lead_variants']:r for r in read(brain) if r['pair_id']=='insomnia__adhd'};assert len(rows)==len(cands)==6 and {r['variant'] for r in rows}==set(cands)
original={};inputs=json.loads(prov.read_text())['input_sha256']
for x in rows:
 d={}
 for t in ('insomnia','adhd'):
  name=cands[x['variant']]['candidate_locus_id']+'__'+t+'.tsv';file=src/name;assert sha(file)==inputs[name];v=[r for r in read(file) if r['SNP']==x['variant']];assert len(v)==1;z=v[0];assert {z['A1'],z['A2']}=={x['effect_allele'],x['other_allele']};d[t]=float(z['BETA'])*(1 if z['A1']==x['effect_allele'] else -1)
 original[x['variant']]=d
positions={(r['chr'],r['position38']) for r in rows};found=collections.defaultdict(list);n=0
with gzip.open(mvp,'rt') as f:
 header=next(f).rstrip('\n').split('\t');ci=header.index('chromosome');pi=header.index('base_pair_location')
 for line in f:
  n+=1;v=line.rstrip('\n').split('\t')
  if (v[ci],v[pi]) in positions:found[(v[ci],v[pi])].append(dict(zip(header,v)))
prior={r['variant']:r for r in json.loads(previous.read_text())['query_records']};out=[];audit=[]
for x in rows:
 snp=x['variant'];ea=x['effect_allele'];oa=x['other_allele'];q=[r for r in found[(x['chr'],x['position38'])] if r['rsid']==snp and {r['effect_allele'],r['other_allele']}=={ea,oa}];raw=prior[snp]['symptom_raw_rows'];s=[v.split('\t') for v in raw if len(v.split('\t'))==12 and v.split('\t')[8]==snp and {v.split('\t')[2],v.split('\t')[3]}=={ea,oa}]
 z={'variant':snp,'chr':x['chr'],'pos38':x['position38'],'effect_allele':ea,'status':'NOT_ESTIMABLE','mvp_logor_aligned':'','mvp_p':'','mvp_same_direction':'','symptom_beta_aligned':'','symptom_p':'','symptom_same_direction':'','joint_max_p':'','joint_times_six':''}
 if len(q)==len(s)==1:
  try:
   r=q[0];v=s[0];odds,mp,f,beta,sp=map(float,(r['odds_ratio'],r['p_value'],r['effect_allele_frequency'],v[4],v[7]));assert all(map(math.isfinite,(odds,mp,f,beta,sp))) and odds>0 and 0<=mp<=1 and 0<=sp<=1 and 0<=f<=1
   lb=math.log(odds)*(1 if r['effect_allele']==ea else -1);sb=beta*(1 if v[2]==ea else -1);m=max(mp,sp)
   z.update(status='MATCHED',mvp_logor_aligned=lb,mvp_p=mp,mvp_same_direction='YES' if lb*original[snp]['insomnia']>0 else 'NO_OR_ZERO',symptom_beta_aligned=sb,symptom_p=sp,symptom_same_direction='YES' if sb*original[snp]['adhd']>0 else 'NO_OR_ZERO',joint_max_p=m,joint_times_six=min(1,6*m))
  except (ValueError,AssertionError):z['status']='INVALID_STATISTICS'
 elif len(q)>1 or len(s)>1:z['status']='AMBIGUOUS'
 out.append(z);audit.append({'variant':snp,'mvp_rows_at_position':len(found[(x['chr'],x['position38'])]),'mvp_exact_matches':len(q),'symptom_exact_matches':len(s),'mvp_exact_rows':q})
with (p/'six_leads.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=out[0].keys(),delimiter='\t');w.writeheader();w.writerows(out)
rec={'mvp_path':str(mvp),'mvp_sha256':actual,'mvp_rows':n,'source_prior_receipt_sha256':sha(mvp_receipt),'six_lead_input_sha256':sha(six),'six_lead_query_receipt_sha256':sha(previous),'frozen_candidate_sha256':sha(brain),'audit':audit,'note':'Complete MVP file hash and stream checked; original per-variant PGC file data not redistributed.'};(p/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(out,indent=2))
