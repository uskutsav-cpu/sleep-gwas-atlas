import gzip,csv,hashlib,json,pathlib,collections
p=pathlib.Path('work');ids=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv'))); idsmap={x['variant']:x for x in ids};stats={};matches=[]
expected={'PGS005228':'584c4839ba4d8caa299ff558ed161259','PGS005229':'168e5c11469c0b27668a4c4d733276b1'}
for name,md5 in expected.items():
 f=p/(name+'_complete.txt.gz');h=hashlib.md5(); h2=hashlib.sha256()
 with f.open('rb') as b:
  for chunk in iter(lambda:b.read(8*1024*1024),b''):h.update(chunk);h2.update(chunk)
 assert h.hexdigest()==md5
 c=collections.Counter();found=collections.defaultdict(list);header=[]
 with gzip.open(f,'rt',newline='') as g:
  lines=(s for s in g if not s.startswith('#')); reader=csv.DictReader(lines,delimiter='\t'); cols=reader.fieldnames
  for row in reader:
   c['rows']+=1;ea=row.get('effect_allele','');oa=row.get('other_allele','');
   if ea and oa and ea==oa:c['same_alleles']+=1
   if not ea:c['empty_effect_allele']+=1
   if not oa:c['empty_other_allele']+=1
   if row['rsID'] in idsmap:found[row['rsID']].append(row)
 for v,x in idsmap.items():
  t=found.get(v,[])
  if not t:matches.append({'score_id':name,'variant':v,'matches':0,'record_status':'not_found; not interpreted as zero weight','source_url':'https://www.pgscatalog.org/score/'+name+'/'});continue
  for s in t:
   coord=(s.get('chr_name')==x['chr'] and s.get('chr_position')==x['position_GRCh37_inferred'])
   matches.append({'score_id':name,'variant':v,'matches':len(t),'score_chr':s.get('chr_name'),'score_position':s.get('chr_position'),'prototype_chr':x['chr'],'prototype_position':x['position_GRCh37_inferred'],'coordinate_match':coord,'score_effect_allele':s.get('effect_allele'),'score_other_allele':s.get('other_allele'),'score_effect_weight':s.get('effect_weight'),'prototype_effect_allele':x['counted_allele'],'prototype_other_allele':x['other_allele'],'record_status':'source other allele identical to effect allele; harmonization unresolved' if s.get('effect_allele') and s.get('effect_allele')==s.get('other_allele') else 'matched rsID only; coordinate/allele/source QC required before use','source_url':'https://www.pgscatalog.org/score/'+name+'/'})
 stats[name]={'source_url':'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/'+name+'/ScoringFiles/'+name+'.txt.gz','source_md5':h.hexdigest(),'source_sha256':h2.hexdigest(),'header':cols,'counts':dict(c),'prototype_variants_with_rsid_match':len(found),'prototype_variants_with_duplicate_rsid_rows':sum(len(v)>1 for v in found.values())}
 print(name,stats[name],flush=True)
keys=['score_id','variant','matches','score_chr','score_position','prototype_chr','prototype_position','coordinate_match','score_effect_allele','score_other_allele','score_effect_weight','prototype_effect_allele','prototype_other_allele','record_status','source_url']
with open('outputs/published_fi_ffs_pgs_candidate_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(matches)
pathlib.Path('outputs/published_fi_ffs_pgs_source_audit.json').write_text(json.dumps(stats,indent=2))
