import csv,hashlib,json,pathlib,re,subprocess
p=pathlib.Path(__file__).resolve().parent;root=pathlib.Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas');munged=pathlib.Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged');ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');log=p/'insomnia_adhd.log'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
source=p.parent/'atlas_396_source.tsv'
with source.open() as f:rows=[x for x in csv.DictReader(f,delimiter='\t') if x['pair_id']=='insomnia__adhd']
assert len(rows)==1;a=rows[0];txt=log.read_text()
m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',txt);assert m, 'missing correlation section'
rg,se,z,pvalue=map(float,m.groups());assert rg==float(a['global_rg']) and se==float(a['global_rg_se']) and pvalue==float(a['global_rg_p']);assert 'Analysis finished at' in txt
required=[munged/'insomnia.sumstats.gz',munged/'adhd.sumstats.gz',root/'ldsc/ldsc.py',source,log,p/'run.sh']
required += [ref/f'{i}.l2.{suffix}' for i in range(1,23) for suffix in ('ldscore.gz','M_5_50')]
required += sorted((root/'ldsc/ldscore').glob('*.py'))
result={'rg':rg,'se':se,'z_printed':z,'p_printed':pvalue,'atlas_q_396':float(a['global_rg_fdr_all_396']),'matches_atlas_reported_precision':True,'input_frozen_hashes':{'insomnia':'b4f66dde69d0c793f77b674dca45526659e69e6fbb31a5a4f3e57870077e6fc4','adhd':'dadb8adcabcb1e8e2d59c8f85fbc4fdd48d38753668bf0a377e938e4ad065b4c'},'source_sha256':{str(x):sha(x) for x in required},'same_data_reexecution':True,'original_atlas_run_provenance_verified':False,'ldsc_runtime':'archived conda-forge Python 3.9.23, numpy 1.21.5, pandas 1.3.3, scipy 1.7.3','warning':'The log reports a compression runtime warning; both input files were read and the run completed. No independent cohort or local-QC resolution.'}
for k,v in result['input_frozen_hashes'].items():assert result['source_sha256'][str(munged/f'{k}.sumstats.gz')]==v
(p/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('verified same printed rg, se and p; reference and code files hashed:',len(required))
