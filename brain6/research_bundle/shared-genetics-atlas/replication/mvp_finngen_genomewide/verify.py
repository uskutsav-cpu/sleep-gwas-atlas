import json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent;j=json.loads((p/'stream_receipt.json').read_text());s=(p/'mvp_adhd_clean.log').read_text();c=j['counts'];assert c['rows']==19703815 and c['retained']==1124041 and c['duplicate_ids_excluded']==140 and c['rounded_OR_no_direction']==27416
assert s.count('1124041 SNPs with valid alleles.')==1
m=re.search(r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\s+Z-score: ([\d.\-]+)\s+P: ([\d.eE\-]+)',s);assert m
rg,se,z,pv=map(float,m.groups());papprox=math.erfc(abs(z)/math.sqrt(2));assert abs(papprox-pv)/pv<.001 and abs(rg/se-z)<.01
n=4*78566*329572/(78566+329572);assert abs(j['n_effective_min']-n)<1e-6 and abs(j['n_effective_max']-n)<1e-6
(p/'verification.json').write_text(json.dumps({'retained':c['retained'],'duplicate_ids_excluded':c['duplicate_ids_excluded'],'effective_n':n,'rg':rg,'se':se,'logged_p':pv,'reconstructed_p_from_rounded_z':papprox,'scope':'stream/log count and approximate tail arithmetic; not independent source or LDSC calibration validation'},indent=2)+'\n');print('Verified retained count, sample approximation, and reported p-value arithmetic')
