import pathlib,hashlib,gzip,csv,math,collections,json,openpyxl
p=pathlib.Path(__file__).resolve().parent
xlsx=p/'published_supplement.xlsx'; source=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
hash_source=sha(source);assert hash_source=='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049'
w=openpyxl.load_workbook(xlsx,read_only=True,data_only=True);s=w['ST2'];r=iter(s.values);next(r);headers=next(r);lead=[]
for number,row in enumerate(r,3):
 if row[0] is None:continue
 d=dict(zip(headers,row))
 if d['indep_lead_snp'] is True:d['source_excel_row']=number;lead.append(d)
assert len(lead)==42
ids={x['SNP'] for x in lead if isinstance(x['SNP'],str) and x['SNP'].startswith('rs')};found=collections.defaultdict(list)
with gzip.open(source,'rt') as f:
 header=f.readline().rstrip('\n').split('\t');idx={x:i for i,x in enumerate(header)};n=0
 for line in f:
  n+=1;cols=line.rstrip('\n').split('\t');rs=cols[idx['rsids']]
  if rs in ids:found[rs].append(cols)
assert n==21327010
out=[]
for x in lead:
 a=str(x['finngen_alt']).upper();b=str(x['finngen_ref']).upper();rs=x['SNP'];v={'excel_row':x['source_excel_row'],'rsid':rs,'chr':x['chr'],'pos_reported':x['pos'],'gene_annotation':x['gene'],'hfrs_alt':a,'hfrs_ref':b,'hfrs_beta':x['finngen_beta'],'hfrs_p':x['finngen_p'],'hfrs_af_alt':x['finngen_af_alt'],'ukb_alt':x['ukb_alt'],'ukb_ref':x['ukb_ref'],'ukb_beta':x['ukb_beta'],'ukb_p':x['ukb_p'],'r13_id_records':len(found[rs])}
 if not isinstance(rs,str) or not rs.startswith('rs'):v['status']='NO_RSID'
 elif len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:v['status']='NON_BIALLELIC_SNV'
 elif frozenset((a,b)) in (frozenset('AT'),frozenset('CG')):v['status']='PALINDROMIC_EXCLUDED'
 else:
  matches=[q for q in found[rs] if q[idx['#chrom']]==str(x['chr']) and q[idx['pos']]==str(x['pos']) and {q[idx['ref']].upper(),q[idx['alt']].upper()}=={a,b}]
  v['exact_matches']=len(matches)
  if len(matches)!=1:v['status']='NO_UNIQUE_MATCH'
  else:
   q=matches[0];prob=float(q[idx['pval']]);bet=float(q[idx['beta']]);freq=float(q[idx['af_alt']]);assert 0<=prob<=1 and math.isfinite(bet) and 0<=freq<=1
   sign=1 if q[idx['alt']].upper()==a else -1;v.update(status='MATCH',insomnia_p=prob,insomnia_beta_hfrs_alt=sign*bet,insomnia_af_hfrs_alt=freq if sign==1 else 1-freq,benchmark_p42=min(1,42*prob),same_sign=sign*bet*float(x['finngen_beta'])>0)
 out.append(v)
fields=['excel_row','rsid','chr','pos_reported','gene_annotation','hfrs_alt','hfrs_ref','hfrs_beta','hfrs_p','hfrs_af_alt','ukb_alt','ukb_ref','ukb_beta','ukb_p','r13_id_records','exact_matches','status','insomnia_p','insomnia_beta_hfrs_alt','insomnia_af_hfrs_alt','benchmark_p42','same_sign']
with (p/'lead_followup.tsv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(out)
receipt={'workbook_sha256':sha(xlsx),'insomnia_source_sha256':hash_source,'insomnia_source_rows':n,'author_lead_rows':len(out),'status_counts':dict(collections.Counter(x['status'] for x in out)),'threshold_42_below_005':sum(x.get('benchmark_p42',1)<.05 for x in out),'output_sha256':sha(p/'lead_followup.tsv')}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
