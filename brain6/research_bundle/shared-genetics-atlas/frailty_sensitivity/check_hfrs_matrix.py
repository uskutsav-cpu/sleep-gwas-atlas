import json,pathlib,statistics,math,hashlib
p=pathlib.Path(__file__).resolve().parent;names={'FI_FFS':'FI_FFS_HFRS_common.exact.json','FI_HFRS':'FI_HFRS_common.exact.json','FFS_HFRS':'FFS_HFRS_common.exact.json'};f={}
for k,n in names.items():
 d=json.loads((p/n).read_text());d=d[0] if isinstance(d,list) else d;d['deletes']=[float(x[0] if isinstance(x,list) else x) for x in d['delete_rg']];assert len(d['deletes'])==200 and d['n_snps']==1137927;f[k]=d
assert len({x['row_order_SNP_sha256'] for x in f.values()})==1
assert f['FI_HFRS']['hsq2']['separators']==f['FFS_HFRS']['hsq2']['separators'] and f['FI_HFRS']['hsq2']['delete_total']==f['FFS_HFRS']['hsq2']['delete_total']
def jk(a):m=statistics.mean(a);return math.sqrt(199/200*math.fsum((x-m)**2 for x in a))
for d in f.values():assert abs(jk(d['deletes'])-d['rg_se'])<1e-12
x,y=f['FFS_HFRS'],f['FI_HFRS'];delta=x['rg']-y['rg'];se=jk([a-b for a,b in zip(x['deletes'],y['deletes'])]);res={'rg':{k:{'estimate':d['rg'],'SE':d['rg_se']} for k,d in f.items()},'FFS_HFRS_minus_FI_HFRS':{'delta':delta,'SE':se,'CI95':[delta-1.95996398454*se,delta+1.95996398454*se],'unadjusted_p':math.erfc(abs(delta/se)/math.sqrt(2))},'n_snps':1137927,'snp_order_sha256':f['FI_FFS']['row_order_SNP_sha256'],'hfrs_shared_delete_vector':True,'input_hashes':{k:hashlib.sha256((p/n).read_bytes()).hexdigest() for k,n in names.items()}}
(p/'hfrs_matrix_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
