import csv,json,pathlib,hashlib,math
p=pathlib.Path(__file__).resolve().parent;study=p/'cohort_rows.tsv';atlas=pathlib.Path('outputs/shared-genetics-atlas');files=[atlas/'frailty_coloc/FI_extracted.tsv',atlas/'frailty_coloc_chr5/FI_extracted.tsv'];source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz');readme=source.with_name(source.name.replace('.txt.gz','_README.txt'))
with study.open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==4 and len({x['SNP'] for x in r})==4
meta={}
for f in files:
 with f.open() as z:
  for x in csv.DictReader(z,delimiter='\t'):
   if x['SNP'] in {q['SNP'] for q in r}:assert x['SNP'] not in meta;meta[x['SNP']]=x
assert len(meta)==4
rows=[]
for a in r:
 b=meta[a['SNP']];assert a['CHR']==b['CHR'] and a['BP']==b['BP'] and {a['EA'],a['OA']}=={b['A1'],b['A2']};sign=1 if a['EA']==b['A1'] else -1
 bu=sign*float(a['BETA_UKB1']);su=float(a['SE_UKB1']);raw=sign*float(a['BETA_TG']);st=float(a['SE_TG']);freq=float(a['EAF_TG']);u=float(a['EAF_UKB']);meta_beta=float(b['BETA'])
 for mode in ('README_same_EA','hypothesis_TwinGene_opposite_EA'):
  flip=1 if mode=='README_same_EA' else -1;bt=flip*raw;meta_ivw=(bu/su**2+bt/st**2)/(1/su**2+1/st**2);qz=(bu-bt)/math.sqrt(su**2+st**2);pq=math.erfc(abs(qz)/math.sqrt(2));tg_eaf=freq if flip==1 else 1-freq
  rows.append({'SNP':a['SNP'],'chr':a['CHR'],'bp37':a['BP'],'meta_A1':b['A1'],'meta_A2':b['A2'],'raw_source_EA':a['EA'],'hypothesis':mode,'UKB_beta_meta_A1':bu,'UKB_SE':su,'TG_beta_meta_A1':bt,'TG_SE':st,'meta_reported_beta':meta_beta,'approx_IVW_beta':meta_ivw,'absolute_delta_from_reported':abs(meta_ivw-meta_beta),'cohort_heterogeneity_p_approx':pq,'EAF_UKB_raw_EA':u,'EAF_TG_assumed_raw_EA':tg_eaf,'absolute_EAF_difference':abs(u-tg_eaf)})
with (p/'orientation_sensitivity.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while x:=f.read(8*1024*1024):h.update(x)
 return h.hexdigest()
receipt={'raw_cohort_source':str(source),'raw_cohort_sha256':sha(source),'README_path':str(readme),'README_sha256':sha(readme),'selected_rows_sha256':sha(study),'meta_extracts_sha256':{str(x):sha(x) for x in files},'source_rows':'four exact rsIDs selected from source; verify_source.py streams full source to check uniqueness and saved rows','caveat':'Approximate IVW ignores paper genomic control; allele adapter not author-confirmed.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('saved rows',len(rows))
