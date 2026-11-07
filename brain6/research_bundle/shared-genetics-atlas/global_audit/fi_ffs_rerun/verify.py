import hashlib,json,pathlib,re
p=pathlib.Path(__file__).resolve().parent;root=pathlib.Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas');inp=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-tp-v3-20261009/munged_full_native');ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');log=p/'fi_ffs.log';original=p.parent/'fi_ffs_native_exact.json'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
d=json.loads(original.read_text());assert len(d)==1;a=d[0];txt=log.read_text();m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',txt);assert m
rg,se,z,pvalue=map(float,m.groups());assert abs(rg-a['rg'])<=.00005 and abs(se-a['rg_se'])<=.00005 and abs(z-a['z'])<=.00005 and pvalue==0 and '1157249 SNPs with valid alleles.' in txt and 'Analysis finished at' in txt
required=[inp/'frailty.sumstats.gz',inp/'fried_frailty_score.sumstats.gz',root/'ldsc/ldsc.py',log,original,p/'run.sh']+[ref/f'{i}.l2.{suffix}' for i in range(1,23) for suffix in ('ldscore.gz','M_5_50')]+sorted((root/'ldsc/ldscore').glob('*.py'))
r={'rg_printed':rg,'se_printed':se,'z_printed':z,'original_rg':a['rg'],'original_se':a['rg_se'],'original_z':a['z'],'original_snps':a['n_snps'],'matching_printed_values':True,'p_underflow':True,'original_sample_overlap':'unknown','source_sha256':{str(f):sha(f) for f in required},'runtime':'archived Python 3.9.23 / numpy 1.21.5 / pandas 1.3.3 / scipy 1.7.3','limitations':'same input paths as archived exact native record, source bytes not independently hashed in that record; not independent cohort replication; no causal inference'}
(p/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print('Verified printed rg, SE, z, SNP count and completed log; source hashes',len(required))
