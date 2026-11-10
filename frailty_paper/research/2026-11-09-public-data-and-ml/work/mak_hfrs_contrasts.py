import json,math,statistics,pathlib,shutil
names=['old_FI','old_FFS','mak_FI','mak_FP'];a=[json.load(open('work/mak_hfrs_'+n+'.exact.json'))[0] for n in names];assert len(set(d['row_order_SNP_sha256'] for d in a))==1 and len(set(d['n_snps'] for d in a))==1
assert all(d['hsq1']['separators']==a[0]['hsq1']['separators']==d['hsq2']['separators'] for d in a)
assert all(d['hsq2']['total']==a[0]['hsq2']['total'] and d['hsq2']['delete_total']==a[0]['hsq2']['delete_total'] for d in a)
dels=[[x[0] for x in d['delete_rg']] for d in a];assert all(len(x)==200 for x in dels)
def se(v):
 m=statistics.mean(v);return math.sqrt((len(v)-1)/len(v)*sum((x-m)**2 for x in v))
for x,d in zip(dels,a):assert math.isclose(se(x),d['rg_se'],rel_tol=1e-10)
ret={'n_snps':a[0]['n_snps'],'hash':a[0]['row_order_SNP_sha256'],'n_blocks':200,'traits':{n:{'rg':d['rg'],'se':d['rg_se'],'frailty_h2':d['hsq1']['total'],'frailty_h2_se':d['hsq1']['total_se'],'HFRS_h2':d['hsq2']['total'],'HFRS_h2_se':d['hsq2']['total_se'],'frailty_intercept':d['hsq1']['intercept'],'HFRS_intercept':d['hsq2']['intercept'],'cross_intercept':d['gencov']['intercept']} for n,d in zip(names,a)},'contrasts':{}}
for name,co in [('old_FFS_minus_FI',[-1,1,0,0]),('mak_FP_minus_FI',[0,0,-1,1]),('difference_Mak_minus_old',[1,-1,-1,1])]:
 p=sum(c*d['rg'] for c,d in zip(co,a));dp=[sum(c*dels[i][j] for i,c in enumerate(co)) for j in range(200)];s=se(dp);ret['contrasts'][name]={'estimate':p,'se':s,'ci_low':p-1.95996398454*s,'ci_high':p+1.95996398454*s,'p_unadjusted':math.erfc(abs(p/s)/math.sqrt(2))}
pathlib.Path('outputs/mak_hfrs_matched_contrasts.json').write_text(json.dumps(ret,indent=2)+'\n')
for n in names:
 for suff in ('.exact.json','.log'):shutil.copyfile('work/mak_hfrs_'+n+suff,'outputs/mak_hfrs_'+n+suff)
print(json.dumps(ret,indent=2))
