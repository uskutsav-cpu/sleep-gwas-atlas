import csv, pathlib, gzip, collections, math, hashlib, json
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent
inp=p.parent/'completed_audit/candidate_positions.tsv'
chain=p/'hg19ToHg38.over.chain.gz'
raw=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
sources={'HFRS':raw/'GCST90707291.tsv','MVP':raw/'GCST90475826.tsv.gz'}
expected={'HFRS':'38d9434f50a4493def9d2c9d1cc1b4b3','MVP':'46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835'}
def digest(path,algo):
 h=hashlib.new(algo)
 with path.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
for name,path in sources.items():
 got=digest(path,'md5' if name=='HFRS' else 'sha256');assert got==expected[name],(name,got);print('Verified',name,flush=True)
with inp.open() as f: candidates=list(csv.DictReader(f,delimiter='\t'))
assert len(candidates)==82 and len({x['SNP'] for x in candidates})==82
lo=LiftOver(str(chain));byid={x['SNP']:x for x in candidates};found={name:collections.defaultdict(list) for name in sources}
for name,path in sources.items():
 op=gzip.open if path.suffix=='.gz' else open
 with op(path,'rt') as f:
  reader=csv.DictReader(f,delimiter='\t'); count=0
  for x in reader:
   count+=1; rs=x['rs_id' if name=='HFRS' else 'rsid']
   if rs in byid:found[name][rs].append(x)
 print(name,'rows',count,'ids',len(found[name]),flush=True)
rows=[];statuses=collections.Counter()
for x in candidates:
 rs=x['SNP']; a=x['A1'].upper();b=x['A2'].upper(); c=x['CHR'];pos=x['BP_GRCH37']; row={'pair':x['pair'],'rsid':rs,'chr37':c,'pos37':pos,'a1':a,'a2':b,'original_placo_p':x['p_placo_plus']}
 if not c.isdigit() or not pos.isdigit():row['status']='MISSING_POSITION';rows.append(row);continue
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b or frozenset((a,b)) in (frozenset('AT'),frozenset('CG')):row['status']='ALLELE_EXCLUDED';rows.append(row);continue
 lift=lo.convert_coordinate('chr'+c,int(pos)-1)
 if len(lift)!=1 or lift[0][0]!='chr'+c or lift[0][2]!='+':row['status']='LIFTOVER_EXCLUDED';rows.append(row);continue
 pos38=lift[0][1]+1;row['pos38']=pos38;invalid=[]
 for name in sources:
  source_rows=found[name][rs]; matched=[y for y in source_rows if y['chromosome']==c and y['base_pair_location']==str(pos38) and {y['effect_allele'].upper(),y['other_allele'].upper()}=={a,b}]
  row[name+'_id_rows']=len(source_rows);row[name+'_match_rows']=len(matched)
  if len(matched)!=1:invalid.append(name+'_NO_UNIQUE_MATCH');continue
  y=matched[0]
  try:
   prob=float(y['p_value']);effect=float(y['beta'] if name=='HFRS' else y['odds_ratio']);freq=float(y['effect_allele_frequency'])
   if not(0<=prob<=1 and math.isfinite(effect) and 0<=freq<=1) or (name=='MVP' and effect<=0):raise ValueError()
  except ValueError:invalid.append(name+'_INVALID_NUMERIC');continue
  aligned=(effect if name=='HFRS' else math.log(effect))*(1 if y['effect_allele'].upper()==a else -1)
  row[name+'_p']=prob;row[name+'_a1_beta_or_logodds']=aligned;row[name+'_a1_frequency']=freq if y['effect_allele'].upper()==a else 1-freq
 if invalid:row['status']=';'.join(invalid)
 else:
  row['status']='ELIGIBLE';row['max_p']=max(row['MVP_p'],row['HFRS_p']);row['max_p_times_82']=min(1,row['max_p']*82);row['same_sign']=(row['MVP_a1_beta_or_logodds']*row['HFRS_a1_beta_or_logodds']>0)
 rows.append(row)
for row in rows:statuses[row['status']]+=1
fields=['pair','rsid','chr37','pos37','pos38','a1','a2','original_placo_p','status','MVP_id_rows','MVP_match_rows','MVP_p','MVP_a1_beta_or_logodds','MVP_a1_frequency','HFRS_id_rows','HFRS_match_rows','HFRS_p','HFRS_a1_beta_or_logodds','HFRS_a1_frequency','max_p','max_p_times_82','same_sign']
with (p/'followup.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
summary={'candidate_sha256':digest(inp,'sha256'),'chain_sha256':digest(chain,'sha256'),'source_digests':{'HFRS_md5':expected['HFRS'],'MVP_sha256':expected['MVP']},'n':len(rows),'status_counts':dict(statuses),'eligible':sum(x['status']=='ELIGIBLE' for x in rows),'benchmark_below_005':sum(x.get('max_p_times_82',1)<.05 for x in rows),'followup_sha256':digest(p/'followup.tsv','sha256')}
(p/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
