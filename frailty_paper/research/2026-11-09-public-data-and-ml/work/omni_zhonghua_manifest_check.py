import csv,gzip,hashlib,json
from pathlib import Path
p=Path('outputs');f=Path('work/omni_zhonghua_v11_manifest.csv.gz');groups={};score_pairs={}
for label,fn in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((p/fn).open()):
  c=r.get('chr',r.get('chr_GRCh37_assumed'));bp=r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed'));groups.setdefault(r['variant'],[]).append((label,c,bp));score_pairs[(label,r['variant'])]=''.join(sorted((r['counted_allele'],r['other_allele'])))
assert len(groups)==33
byp={};byn={};n=0
with gzip.open(f,'rt',encoding='utf-8-sig',newline='') as h:
 for line in h:
  if line.strip()=='[Assay]':break
 for r in csv.DictReader(h):
  if not r or not r.get('Name') or r['Name']=='[Controls]':break
  n+=1
  if r['Name'] in groups:byn.setdefault(r['Name'],[]).append(r)
  if any((r['Chr'],r['MapInfo'])==(c,bp) for entries in groups.values() for _,c,bp in entries):byp.setdefault((r['Chr'],r['MapInfo']),[]).append(r)
out=[]
complement=str.maketrans('ACGT','TGCA')
def plus_pair(row):
 raw=row['SNP']; strand=row['RefStrand']
 if len(raw)!=5 or raw[0]!='[' or raw[2]!='/' or raw[4]!=']' or raw[1] not in 'ACGT' or raw[3] not in 'ACGT' or strand not in ('+','-'):return None
 pair=raw[1]+raw[3];return ''.join(sorted(pair if strand=='+' else pair.translate(complement)))
for rs,entries in groups.items():
 for label,c,bp in entries:
  hits=byn.get(rs,[]);same=[x for x in hits if x['Chr']==c and x['MapInfo']==bp];co=byp.get((c,bp),[])
  coordpairs=[plus_pair(x) for x in co]
  out.append({'list':label,'variant':rs,'GRCh37_chr':c,'GRCh37_pos':bp,'manifest_exact_name_rows':len(hits),'manifest_name_and_coordinate_rows':len(same),'any_coordinate_rows':len(co),'manifest_coordinate_names':';'.join(x['Name'] for x in co),'manifest_coordinate_plus_pairs':';'.join(x or 'unresolved' for x in coordpairs),'score_pair':score_pairs[(label,rs)],'plus_pair_matches_score':score_pairs[(label,rs)] in coordpairs,'manifest_snp_raw':';'.join(x['SNP'] for x in same),'manifest_ilmnstrand':';'.join(x['IlmnStrand'] for x in same),'manifest_refstrand':';'.join(x['RefStrand'] for x in same),'caveat':'Plus pair inferred using SNP and RefStrand; array version, genotypes, build and coding are not independently verified.'})
with (p/'stage6_clhls_omni_v11_coverage.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
res={'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL19864','file_url':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19864/suppl/GPL19864_HumanOmniZhongHua-8-v1-1-B.csv.gz','local_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'manifest_assay_rows_scanned':n,'results':{tag:{'candidates':sum(r['list']==tag for r in out),'name_and_coordinate_matches':sum(r['list']==tag and r['manifest_name_and_coordinate_rows']>=1 for r in out),'coordinate_only_or_other_name':sum(r['list']==tag and r['any_coordinate_rows']>=1 and r['manifest_name_and_coordinate_rows']==0 for r in out),'coordinate_plus_pair_matches':sum(r['list']==tag and r['plus_pair_matches_score'] for r in out)} for tag in ['score15','score22']},'caveat':'v1.1 representative manifest only; plus-strand pairs inferred from RefStrand, chip revision, performance, target-sample genotype and imputation quality are unverified.'}
(p/'stage6_clhls_omni_v11_coverage.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
