import json,math,pathlib,shutil
p=pathlib.Path('work');out=pathlib.Path('outputs');traits=['UKBFI','FI','FFS'];d={t:json.loads((p/f'{t}_R13_ukbcommon.exact.json').read_text())[0] for t in traits}
assert len({x['row_order_SNP_sha256'] for x in d.values()})==1
for key in ['hsq1','hsq2','gencov']:
 assert len({tuple(x[key]['separators']) for x in d.values()})==1
assert len({json.dumps(x['hsq2']['delete_total']) for x in d.values()})==1
assert len({x['hsq2']['total'] for x in d.values()})==1
for x in d.values():
 a=[v[0] for v in x['delete_rg']];n=len(a);m=sum(a)/n;se=math.sqrt((n-1)/n*sum((v-m)**2 for v in a));assert abs(se-x['rg_se'])<1e-10

def compare(a,b):
 x,y=d[a],d[b];v=[u[0]-z[0] for u,z in zip(x['delete_rg'],y['delete_rg'])];n=len(v);m=sum(v)/n;se=math.sqrt((n-1)/n*sum((z-m)**2 for z in v));delta=x['rg']-y['rg'];return {'delta':delta,'se':se,'ci95':[delta-1.96*se,delta+1.96*se],'p_two_sided_normal_unadjusted':math.erfc(abs(delta/se)/math.sqrt(2))}
res={'trait': 'FinnGen R13 clinical insomnia','n_snps':d['UKBFI']['n_snps'],'snp_order_sha256':d['UKBFI']['row_order_SNP_sha256'],'estimates':{t:{'rg':d[t]['rg'],'se':d[t]['rg_se'],'frailty_h2':d[t]['hsq1']['total'],'frailty_h2_se':d[t]['hsq1']['total_se']} for t in traits},'contrast_ukbfi_minus_metafi':compare('UKBFI','FI'),'contrast_ffs_minus_ukbfi':compare('FFS','UKBFI'),'r13_insomnia_h2':d['FI']['hsq2']['total'],'r13_insomnia_h2_se':d['FI']['hsq2']['total_se'],'r13_insomnia_intercept':d['FI']['hsq2']['intercept'],'r13_insomnia_intercept_se':d['FI']['hsq2']['intercept_se'],'same_insomnia_delete_values_and_block_boundaries':True,'block_count':200,'s_sumstats_sha256':'3d7073182abbddfbd83eefbd1187ac999a5d49321ec982d19c851ac1027826c6'}
(out/'ukb_fi_finngen_insomnia_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
for t in traits:shutil.copyfile(p/f'{t}_R13_ukbcommon.exact.json',out/f'{t}_R13_ukbcommon.exact.json')
print(json.dumps(res,indent=2))
