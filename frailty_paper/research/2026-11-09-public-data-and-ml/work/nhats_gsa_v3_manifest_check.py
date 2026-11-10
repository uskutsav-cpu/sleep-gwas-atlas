import csv,zipfile,io,hashlib,json
from pathlib import Path
p=Path('outputs');f=Path('work/gsa_v3_grch37_manifest.zip');groups={};score_pairs={}
for label,fn in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((p/fn).open()):
  c=r.get('chr',r.get('chr_GRCh37_assumed'));bp=r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed'));groups.setdefault(r['variant'],[]).append((label,c,bp));score_pairs[(label,r['variant'])]=''.join(sorted((r['counted_allele'],r['other_allele'])))
assert len(groups)==33
byp={};byn={};n=0;positive_control=0
with zipfile.ZipFile(f) as z, z.open('GSA-24v3-0_A1.csv') as raw, io.TextIOWrapper(raw,encoding='utf-8-sig',newline='') as h:
 for line in h:
  if line.strip()=='[Assay]':break
 for r in csv.DictReader(h):
  if not r or not r.get('Name') or r['Name']=='[Controls]':break
  n+=1
  if r['Name']=='rs7412' and r['Chr']=='19' and r['MapInfo']=='45412079':positive_control+=1
  keyname=r['Name']; keyname=keyname[4:] if keyname.startswith('GSA-rs') else keyname
  if keyname in groups:byn.setdefault(keyname,[]).append(r)
  if any((r['Chr'],r['MapInfo'])==(c,bp) for entries in groups.values() for _,c,bp in entries):byp.setdefault((r['Chr'],r['MapInfo']),[]).append(r)
assert positive_control>=1
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
with (p/'stage6_nhats_gsa_v3_coverage.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
res={'source':'https://support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html','file_url':'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A1-manifest-file-csv.zip','local_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'manifest_assay_rows_scanned':n,'positive_control_rs7412_19_45412079_rows':positive_control,'results':{tag:{'candidates':sum(r['list']==tag for r in out),'name_and_coordinate_matches':sum(r['list']==tag and r['manifest_name_and_coordinate_rows']>=1 for r in out),'coordinate_only_or_other_name':sum(r['list']==tag and r['any_coordinate_rows']>=1 and r['manifest_name_and_coordinate_rows']==0 for r in out),'coordinate_plus_pair_matches':sum(r['list']==tag and r['plus_pair_matches_score'] for r in out)} for tag in ['score15','score22']},'caveat':'GSA-24v3-0_A1 manifest only; plus-strand pairs inferred from RefStrand, chip revision, performance, target-sample genotype and imputation quality are unverified.'}
(p/'stage6_nhats_gsa_v3_coverage.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
