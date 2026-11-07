import csv, gzip, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
p=Path(__file__).resolve().parent
source=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
parent=Path('/Volumes/Brain6Work/chr7-three-trait-input.tsv')
out=Path('/Volumes/Brain6Work/chr7-three-trait-ukb-sensitivity.tsv')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
with parent.open() as f: old=list(csv.DictReader(f,delimiter='\t'))
assert len(old)==1055 and len({r['rsid'] for r in old})==len(old)
target={r['rsid']:r for r in old}; found=defaultdict(list); n=0
with gzip.open(source,'rt') as f:
 rd=csv.DictReader(f,delimiter='\t')
 for row in rd:
  n+=1
  if row['SNP'] in target:found[row['SNP']].append(row)
counts=Counter();res=[]
for a in old:
 candidates=found[a['rsid']]
 if len(candidates)!=1:counts['missing_or_duplicate']+=1;continue
 b=candidates[0]
 if b['CHR']!='7' or b['BP']!=a['pos37']:counts['position']+=1;continue
 alleles={b['EA'],b['OA']}
 if alleles!={a['effect_allele'],a['other_allele']} or alleles in ({'A','T'},{'C','G'}):counts['allele']+=1;continue
 try:
  beta=float(b['BETA_UKB1']);se=float(b['SE_UKB1'])
  if not(math.isfinite(beta) and math.isfinite(se) and se>0):raise ValueError()
 except ValueError:counts['statistics']+=1;continue
 res.append(dict(a,ukb_fi_beta=beta*(1 if b['EA']==a['effect_allele'] else -1),ukb_fi_se=se))
assert len(res)>500 and sum(counts.values())+len(res)==len(old)
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=res[0].keys(),delimiter='\t');w.writeheader();w.writerows(res)
receipt={'source':str(source),'source_sha256':sha(source),'source_rows_read':n,'parent':str(parent),'parent_sha256':sha(parent),'intersection':len(res),'excluded':dict(counts),'matched_input':str(out),'matched_input_sha256':sha(out),'source_note':'UKB1 effect uses source EA as documented; combined data remain on SSD under source terms.'}
(p/'matching_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
