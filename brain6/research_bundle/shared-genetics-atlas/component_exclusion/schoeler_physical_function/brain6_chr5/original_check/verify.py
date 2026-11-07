import json,math
from pathlib import Path
p=Path(__file__).resolve().parent
source=json.loads((p/'receipt.json').read_text())['rows'];harm=json.loads((p.parent/'receipt.json').read_text())['rows']
expected=[('rs2431108',103947968),('rs77960',103964585)];out=[]
for rs,pos in expected:
 a=[r for r in source if int(r['base_pair_location'])==pos];b=[x for x in harm if x['expected_rsid']==rs]
 assert len(a)==len(b)==1 and b[0]['exact_matching_count']==1
 x=a[0];y=b[0]['all_rows_at_position'][0]
 assert {x['effect_allele'],x['other_allele']}=={y['effect_allele'],y['other_allele']}
 flip=x['effect_allele']!=y['effect_allele'];assert flip and x['effect_allele']==y['other_allele']
 for k in ('standard_error','p_value','n'):assert math.isclose(float(x[k]),float(y[k]),rel_tol=1e-12)
 assert math.isclose(float(x['beta']),-float(y['beta']),rel_tol=1e-12)
 assert math.isclose(float(x['effect_allele_frequency'])+float(y['effect_allele_frequency']),1,abs_tol=1e-12)
 out.append({'rsid_from_harmonized_release':rs,'original_grch37_position':pos,'effect_allele_swapped':flip,'beta_sign_swapped':True,'frequency_complementary':True,'se_p_n_equal':True})
(p/'verification.json').write_text(json.dumps({'scope':'partial selected-row source comparison; does not verify full harmonized file','rows':out},indent=2)+'\n');print(json.dumps(out,indent=2))
