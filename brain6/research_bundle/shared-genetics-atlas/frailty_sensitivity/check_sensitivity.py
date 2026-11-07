import json, pathlib,math,statistics,hashlib
p=pathlib.Path(__file__).resolve().parent
fits={}
for trait in ('UKBFI','FI','FFS'):
 d=json.loads((p/(trait+'_R13_ukbcommon.exact.json')).read_text());d=d[0] if isinstance(d,list) else d;d['flat_delete']=[float(x[0] if isinstance(x,list) else x) for x in d['delete_rg']];assert len(d['flat_delete'])==200 and d['n_snps']==1134953;fits[trait]=d
assert len({d['row_order_SNP_sha256'] for d in fits.values()})==1
first=fits['FI']['hsq2'];assert all(d['hsq2']['separators']==first['separators'] and d['hsq2']['delete_total']==first['delete_total'] and d['hsq2']['total']==first['total'] for d in fits.values())
def jack(a):
 m=statistics.mean(a);return math.sqrt(199/200*math.fsum((x-m)**2 for x in a))
for d in fits.values():assert abs(jack(d['flat_delete'])-d['rg_se'])<1e-12
results={}
for a,b in [('UKBFI','FI'),('FFS','UKBFI')]:
 x,y=fits[a],fits[b];delta=x['rg']-y['rg'];se=jack([u-v for u,v in zip(x['flat_delete'],y['flat_delete'])]);results[a+'_minus_'+b]={'difference':delta,'SE':se,'CI95':[delta-1.95996398454*se,delta+1.95996398454*se],'p_two_sided_unadjusted':math.erfc(abs(delta/se)/math.sqrt(2))}
receipt={'estimates':{k:{'rg':d['rg'],'se':d['rg_se']} for k,d in fits.items()},'contrasts':results,'matched_snps':1134953,'row_order_sha256':first and fits['FI']['row_order_SNP_sha256'],'aligned_blocks':200,'same_insomnia_h2_delete_vector':True,'inputs_sha256':{k:hashlib.sha256((p/(k+'_R13_ukbcommon.exact.json')).read_bytes()).hexdigest() for k in fits}}
(p/'audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
