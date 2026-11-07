import csv,math,pathlib,collections,hashlib,json,statistics
out=pathlib.Path('outputs/shared-genetics-atlas/frailty_coloc_chr5/ffs_readiness');matched=out.parent/'matched.tsv';source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
with matched.open() as f:m={r['rsid']:r for r in csv.DictReader(f,delimiter='\t')}
assert len(m)==2337;found=collections.defaultdict(list);cnt=collections.Counter()
with source.open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  cnt['source_rows']+=1;rs=row['SNP']
  if rs not in m:continue
  cnt['matching_rsid_rows']+=1;x=m[rs]
  if row['CHROM']!='5' or int(row['POS'])!=int(x['bp37']):cnt['position_mismatch']+=1;continue
  if {row['A1'].upper(),row['A2'].upper()}!={x['effect_allele'],x['other_allele']}:cnt['allele_mismatch']+=1;continue
  try:
   b=float(row['BETA']);se=float(row['SE']);p=float(row['P']);assert all(map(math.isfinite,(b,se,p))) and se>0 and 0<p<=1
  except (ValueError,AssertionError):cnt['invalid']+=1;continue
  wp=math.erfc(abs(b/se)/math.sqrt(2));
  if wp==0:cnt['wald_underflow']+=1;continue
  found[rs].append({'rsid':rs,'bp37':row['POS'],'FFS_source_A1':row['A1'],'FFS_source_A2':row['A2'],'FFS_source_beta':b,'FFS_source_se':se,'FFS_source_p':p,'FFS_wald_p':wp,'abs_log10_p_difference':abs(math.log10(p)-math.log10(wp)),'FI_p':x['fi_p'],'R13_insomnia_p':x['r13_p']})
rows=[]
for rs in m:
 if len(found[rs])==1:rows.append(found[rs][0])
 elif len(found[rs])>1:cnt['multiple_valid_matches']+=1
 else:cnt['no_valid_match']+=1
cnt['unique_allele_position_matched']=len(rows);diff=sorted(x['abs_log10_p_difference'] for x in rows);assert diff
with (out/'matched_source_stats.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
h=hashlib.sha256()
with source.open('rb') as f:
 for block in iter(lambda:f.read(1<<20),b''):h.update(block)
meta={'ffs_source_sha256':h.hexdigest(),'input_fi_r13_matched_sha256':hashlib.sha256(matched.read_bytes()).hexdigest(),'ffs_source_path':str(source),'ffs_source_bytes':source.stat().st_size,'counts':dict(cnt),'absolute_log10_p_difference':{'median':statistics.median(diff),'p95':diff[int(.95*(len(diff)-1))],'max':max(diff)},'source_has_allele_frequency':False,'FFS_coloc_status':'NOT_ESTIMATED_UNVERIFIED_SOURCE_SCALE_AND_ALLELES'}
(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
