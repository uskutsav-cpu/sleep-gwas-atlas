import csv,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent; atlas=p.parents[2]; root=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas');munged=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged');ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1048576),b''):h.update(b)
 return h.hexdigest()
log=(p/'insomnia_adhd.log').read_text();m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',log);assert m and 'Analysis finished at' in log
rg,se,z,pv=map(float,m.groups());assert rg==.3817 and se==.0521
old=root/'results/track_b/02_independent_global_replication.tsv'
with old.open() as f:r=[x for x in csv.DictReader(f,delimiter='\t') if x['pair_id']=='B'];assert len(r)==1
previous=json.loads((atlas/'global_audit/insomnia_adhd_rerun/receipt.json').read_text()); paths=[munged/'insomnia.sumstats.gz',munged/'track_b_finngen_r13_F5_ADHD.sumstats.gz',root/'ldsc/ldsc.py',old,p/'run.sh',p/'insomnia_adhd.log',p/'console.log']
paths += [ref/f'{i}.l2.{suffix}' for i in range(1,23) for suffix in ('ldscore.gz','M_5_50')];paths += sorted((root/'ldsc/ldscore').glob('*.py'))
hashes={str(x):sha(x) for x in paths};shared=[str(x) for x in paths if str(x) in previous['source_sha256']];assert shared and all(hashes[x]==previous['source_sha256'][x] for x in shared)
rec={'rg':rg,'se':se,'z_printed':z,'p_printed':pv,'previous_aggregate':{'rg':r[0]['replication_rg'],'se':r[0]['replication_SE'],'p':r[0]['replication_P'],'h2':r[0]['h2'],'h2_se':r[0]['h2_SE']},'shared_hashes_matching_prior_ldsc_rerun':len(shared),'source_sha256':hashes,'notes':'Printed rg and SE match historical aggregate; p and phenotype-2 h2/SE differ. Historical original log and raw ADHD input remain unavailable.'}
(p/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print('verified rg/se; previous shared-file hashes agree:',len(shared));print(json.dumps(rec['previous_aggregate']))
