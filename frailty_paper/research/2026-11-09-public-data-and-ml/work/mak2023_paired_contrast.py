import json,math,statistics,pathlib,hashlib,shutil
old=json.load(open('work/mak2023_fourway_old.exact.json'))[0];new=json.load(open('work/mak2023_fourway_mak.exact.json'))[0];assert old['n_snps']==new['n_snps']==1129146 and old['row_order_SNP_sha256']==new['row_order_SNP_sha256'];assert old['hsq1']['separators']==old['hsq2']['separators']==new['hsq1']['separators']==new['hsq2']['separators']
def se(d):
 m=statistics.mean(d);return math.sqrt((len(d)-1)/len(d)*sum((x-m)**2 for x in d))
a=[x[0] for x in old['delete_rg']];b=[x[0] for x in new['delete_rg']];assert len(a)==len(b)==200
assert math.isclose(se(a),old['rg_se'],rel_tol=1e-10) and math.isclose(se(b),new['rg_se'],rel_tol=1e-10)
diff=new['rg']-old['rg'];unc=se([y-x for x,y in zip(a,b)]);n=1.959963984540054
ret={'n_snps':old['n_snps'],'order_hash':old['row_order_SNP_sha256'],'n_blocks':200,'old_rg':old['rg'],'old_se':old['rg_se'],'mak_rg':new['rg'],'mak_se':new['rg_se'],'difference_mak_minus_original':diff,'paired_se':unc,'ci_low':diff-n*unc,'ci_high':diff+n*unc,'p_unadjusted':math.erfc(abs(diff/unc)/math.sqrt(2)),'limitations':'Post hoc; correlated studies with overlapping UKB, source definitions/coding and participant subsets differ. Paired jackknife estimates sampling uncertainty under LDSC assumptions and matched blocks, not the reason for any difference. The support is also restricted by earlier HFRS and R13 availability.'}
pathlib.Path('outputs/mak2023_fi_fp_paired_comparison.json').write_text(json.dumps(ret,indent=2)+'\n')
for name in ('mak2023_fourway_old','mak2023_fourway_mak','mak2023_fi_fp_ldsc'):
 shutil.copyfile('work/'+name+'.exact.json','outputs/'+name+'.exact.json');shutil.copyfile('work/'+name+'.log','outputs/'+name+'.log')
print(json.dumps(ret,indent=2))
