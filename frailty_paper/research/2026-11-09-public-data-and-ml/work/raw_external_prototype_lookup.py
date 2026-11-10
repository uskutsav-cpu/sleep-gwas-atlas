import csv,gzip,hashlib,collections,json,re,pathlib,time
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
inputs={'MVP':(base/'GCST90475826.tsv.gz','46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835'),'R13':(base/'finngen_R13_F5_INSOMNIA.gz','353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049')}
markers=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')));ids={x['variant'] for x in markers}; allrows=[]; receipt={}
for name,(path,expect) in inputs.items():
 start=time.time();h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 assert h.hexdigest()==expect,(name,h.hexdigest())
 found=collections.defaultdict(list);count=0
 with gzip.open(path,'rt',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   count+=1
   aliases=[r['rsid']] if name=='MVP' else re.split(r'[,;\s]+',r['rsids'])
   for v in set(aliases)&ids:found[v].append(r)
 for x in markers:
  for r in found[x['variant']]:
   if name=='MVP':v={'chr':r['chromosome'],'position':r['base_pair_location'],'A1':r['effect_allele'],'A2':r['other_allele'],'effect_value':r['odds_ratio'],'effect_type':'odds_ratio','SE':r['standard_error'],'p':r['p_value'],'frequency':r['effect_allele_frequency'],'N':r['n']}
   else:v={'chr':r['#chrom'],'position':r['pos'],'A1':r['alt'],'A2':r['ref'],'effect_value':r['beta'],'effect_type':'beta_log_odds_as_reported','SE':r['sebeta'],'p':r['pval'],'frequency':r['af_alt'],'N':''}
   allrows.append({'source':name,'variant':x['variant'],'source_id_match_count':len(found[x['variant']]),'source_build':'GRCh38' if name=='R13' else 'source release build to verify',**v,'candidate_GRCh37_chr':x['chr'],'candidate_GRCh37_pos':x['position_GRCh37_inferred'],'candidate_effect_allele':x['counted_allele'],'candidate_other_allele':x['other_allele'],'status':'Raw exact rsID match only; unresolved build, alleles, duplicates and eligibility; no sign comparison or risk score'})
 receipt[name]={'compressed_sha256':expect,'rows_read':count,'unique_ids_found':len(found),'ids_with_multiple_source_rows':{k:len(v) for k,v in found.items() if len(v)>1},'unfound':sorted(ids-set(found)),'elapsed_seconds':round(time.time()-start,1)}
 print(name,json.dumps(receipt[name]),flush=True)
with open('work/raw_external_prototype_matches.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
pathlib.Path('work/raw_external_prototype_receipt.json').write_text(json.dumps(receipt,indent=2))
