import csv,gzip,hashlib,json,math,pathlib,collections,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];brain=root/'evidence/BRAIN6_FINAL_CANDIDATES.tsv';prov=root.parent/'brain6-alphagenome-pilot/results/provenance.json';src=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');chain=root/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz';source_info=json.loads((root/'source_provenance/raw_source_hashes.json').read_text())['FinnGen_R13_insomnia_GRCh38'];raw=pathlib.Path(source_info['path']);symmeta=root/'replication/childhood_symptoms/metadata.yaml';url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90568001-GCST90569000/GCST90568440/harmonised/GCST90568440.h.tsv.gz'
def read(q):
 with q.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
record=json.loads(prov.read_text())['input_sha256'];cands=[r for r in read(brain) if r['pair_id']=='insomnia__adhd'];assert len(cands)==6 and len({r['lead_variants'] for r in cands})==6
lo=LiftOver(str(chain));items=[]
for c in cands:
 s=c['lead_variants'];d={}
 for t in ('insomnia','adhd'):
  name=c['candidate_locus_id']+'__'+t+'.tsv';f=src/name;assert sha(f)==record[name];rows=[r for r in read(f) if r['SNP']==s];assert len(rows)==1;d[t]=rows[0]
 a,b=d['insomnia'],d['adhd'];assert a['CHR']==b['CHR'] and a['POS']==b['POS'] and {a['A1'],a['A2']}=={b['A1'],b['A2']};alleles={a['A1'],a['A2']};assert len(alleles)==2 and alleles.issubset(set('ACGT')) and alleles not in ({'A','T'},{'C','G'})
 lift=lo.convert_coordinate('chr'+a['CHR'],int(a['POS'])-1);assert len(lift)==1 and lift[0][0]=='chr'+a['CHR'] and lift[0][2]=='+'
 items.append({'snp':s,'chr':a['CHR'],'pos37':int(a['POS']),'pos38':lift[0][1]+1,'A1':a['A1'],'A2':a['A2'],'beta_sleep':float(a['BETA']),'beta_adhd':float(b['BETA'])*(1 if b['A1']==a['A1'] else -1)})
assert sha(raw)==source_info['sha256'];targets={(i['chr'],str(i['pos38'])) for i in items};matches=collections.defaultdict(list)
with gzip.open(raw,'rt') as f:
 head=next(f).rstrip('\n').split('\t');assert all(x in head for x in ('#chrom','pos','ref','alt','rsids','pval','beta','sebeta'));ci=head.index('#chrom');pi=head.index('pos')
 n=0
 for line in f:
  n+=1;z=line.rstrip('\n').split('\t')
  if (z[ci],z[pi]) in targets:matches[(z[ci],z[pi])].append(dict(zip(head,z)))
assert n==21327010
sym=pysam.TabixFile(url);cols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','hm_coordinate_conversion','hm_code','variant_id'];out=[];queries=[]
for a in items:
 q=[z for z in matches[(a['chr'],str(a['pos38']))] if z['rsids']==a['snp'] and {z['ref'],z['alt']}=={a['A1'],a['A2']}]
 lines=list(sym.fetch(a['chr'],a['pos38']-1,a['pos38']));parsed=[]
 for l in lines:
  vals=l.split('\t');assert len(vals)==len(cols);parsed.append(dict(zip(cols,vals)))
 t=[z for z in parsed if z['rsid']==a['snp'] and z['base_pair_location']==str(a['pos38']) and {z['effect_allele'],z['other_allele']}=={a['A1'],a['A2']}]
 row={'variant':a['snp'],'chr':a['chr'],'position37':a['pos37'],'position38':a['pos38'],'effect_allele':a['A1'],'other_allele':a['A2'],'status':'NOT_ESTIMABLE','sleep_beta_aligned':'','sleep_se':'','sleep_p':'','sleep_same_direction':'','symptom_beta_aligned':'','symptom_se':'','symptom_p':'','symptom_same_direction':'','max_p':'','max_p_times_six':''}
 if len(q)==len(t)==1:
  try:
   sl=q[0];sy=t[0];sb,ss,sp=map(float,(sl['beta'],sl['sebeta'],sl['pval']));ab,ase,ap,af=map(float,(sy['beta'],sy['standard_error'],sy['p_value'],sy['effect_allele_frequency']))
   assert all(map(math.isfinite,(sb,ss,sp,ab,ase,ap,af))) and ss>0 and ase>0 and 0<=sp<=1 and 0<=ap<=1 and 0<=af<=1
   sb*=1 if sl['alt']==a['A1'] else -1;ab*=1 if sy['effect_allele']==a['A1'] else -1;m=max(sp,ap)
   row.update(status='MATCHED',sleep_beta_aligned=sb,sleep_se=ss,sleep_p=sp,sleep_same_direction='YES' if sb*a['beta_sleep']>0 else 'NO_OR_ZERO',symptom_beta_aligned=ab,symptom_se=ase,symptom_p=ap,symptom_same_direction='YES' if ab*a['beta_adhd']>0 else 'NO_OR_ZERO',max_p=m,max_p_times_six=min(1,6*m))
  except (ValueError,AssertionError):row['status']='INVALID_STATISTICS'
 elif len(q)>1 or len(t)>1:row['status']='AMBIGUOUS'
 out.append(row);queries.append({'variant':a['snp'],'sleep_matches':len(q),'symptom_matches':len(t),'symptom_raw_rows':lines})
with (p/'six_leads.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=out[0].keys(),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'brain_candidates_sha256':sha(brain),'brain_provenance_sha256':sha(prov),'chain_sha256':sha(chain),'finngen_source_sha256':sha(raw),'finngen_total_rows':n,'symptom_url':url,'symptom_metadata_sha256':sha(symmeta),'query_records':queries,'note':'Whole Finnish file SHA matches previously saved source; symptom remote file and index not hashed; no per-person overlap verification.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(out,indent=2))
