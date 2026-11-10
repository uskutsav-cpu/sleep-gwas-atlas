import json, math, pathlib, shutil
p=pathlib.Path('work'); out=pathlib.Path('outputs'); traits=['UKBFI','FI','FFS']
d={t:json.loads((p/f'{t}_HFRS_ukbcommon.exact.json').read_text())[0] for t in traits}
assert len({x['row_order_SNP_sha256'] for x in d.values()})==1
assert len({tuple(x['hsq2']['separators']) for x in d.values()})==1
assert len({json.dumps(x['hsq2']['delete_total']) for x in d.values()})==1
assert len({x['hsq2']['total'] for x in d.values()})==1
for x in d.values():
 a=[v[0] for v in x['delete_rg']]; n=len(a); m=sum(a)/n; se=math.sqrt((n-1)/n*sum((v-m)**2 for v in a)); assert abs(se-x['rg_se'])<1e-10

def contrast(a,b):
 x,y=d[a],d[b]; deletes=[u[0]-v[0] for u,v in zip(x['delete_rg'],y['delete_rg'])]; n=len(deletes); mean=sum(deletes)/n; se=math.sqrt((n-1)/n*sum((v-mean)**2 for v in deletes)); delta=x['rg']-y['rg']; return dict(delta=delta,se=se,ci95=[delta-1.96*se,delta+1.96*se],p_two_sided=math.erfc(abs(delta/se)/math.sqrt(2)))
receipt=json.loads(pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/support.json').read_text())
res={'traits':{t:{k:d[t][k] for k in ('rg','rg_se','n_snps','row_order_SNP_sha256')} for t in traits},'contrast_ukbfi_minus_metafi':contrast('UKBFI','FI'),'contrast_ffs_minus_ukbfi':contrast('FFS','UKBFI'),'hfrs_h2':d['UKBFI']['hsq2']['total'],'hfrs_intercept':d['UKBFI']['hsq2']['intercept'],'hfrs_intercept_se':d['UKBFI']['hsq2']['intercept_se'],'same_hfrs_delete_values_and_blocks':True,'input_receipt':receipt}
(out/'ukb_fi_hfrs_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
for t in traits: shutil.copyfile(p/f'{t}_HFRS_ukbcommon.exact.json',out/f'{t}_HFRS_ukbcommon.exact.json')
print(json.dumps(res,indent=2))
