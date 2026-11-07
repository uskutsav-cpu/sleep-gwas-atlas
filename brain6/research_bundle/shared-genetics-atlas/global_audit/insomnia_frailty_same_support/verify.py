import csv,gzip,hashlib,json,math,pathlib,re,collections
p=pathlib.Path(__file__).resolve().parent;data=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity');ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');root=pathlib.Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas');support=json.loads((data/'support.json').read_text())
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
files=[data/f'{x}.sumstats.gz' for x in ('FI','FFS','HFRS','S')];hashes={str(f):sha(f) for f in files}
for trait in ('FI','FFS','HFRS'):assert hashes[str(data/f'{trait}.sumstats.gz')]==support['output_sha256'][trait]
from contextlib import ExitStack
with ExitStack() as stack:
 readers=[csv.DictReader(stack.enter_context(gzip.open(f,'rt')),delimiter='\t') for f in files];n=0
 for rows in zip(*readers,strict=True):
  keys=[(x['SNP'],x['A1'],x['A2']) for x in rows];assert all(k==keys[0] for k in keys), (n,keys);n+=1
assert n==support['common_snps']==1134953
out=[]
for t in ('FI','FFS','HFRS'):
 log=p/f'{t}_insomnia.log';text=log.read_text();m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',text);assert m and '1134953 SNPs with valid alleles.' in text and 'Analysis finished at' in text
 rg,se,z,prob=map(float,m.groups());assert 0<prob<=1
 if t in ('FI','FFS'):
  original=json.loads((p.parents[1]/'frailty_sensitivity'/f'{t}_R13_ukbcommon.exact.json').read_text())[0]
  assert abs(rg-original['rg'])<.00005 and abs(se-original['rg_se'])<.00005 and math.isclose(prob,original['p'],rel_tol=.00002) and original['n_snps']==n
 out.append({'frailty_trait':t,'sleep_trait':'FinnGen R13 registry insomnia','support_snps':n,'rg':rg,'se':se,'z':z,'p':prob,'bonferroni_p_3':min(1,prob*3)})
with (p/'correlations.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
required=files+[data/'support.json',p/'run.sh',root/'ldsc/ldsc.py']+[p/f'{t}_insomnia.log' for t in ('FI','FFS','HFRS')]+[ref/f'{i}.l2.{k}' for i in range(1,23) for k in ('ldscore.gz','M_5_50')]+sorted((root/'ldsc/ldscore').glob('*.py'))
receipt={'input_support_rows':n,'all_input_snp_and_alleles_ordered_identically':True,'support_source_hashes_FI_FFS_HFRS':'matched','FI_FFS_native_comparison':'printed rg, se, p and counts agree','HFRS_insomnia':'new same-data sensitivity, overlapping FinnGen cohorts; not independent','multiple_testing_family':'three correlations, exploratory Bonferroni; not original project families','source_sha256':{str(f):sha(f) for f in required}}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('verified matched ordered SNP/allele keys',n,'and prior FI/FFS results')
