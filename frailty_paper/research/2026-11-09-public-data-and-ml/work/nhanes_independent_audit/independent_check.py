"""Independent integrity, score-denominator, raw-design and exposure audit."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyreadstat

root=Path(__file__).resolve().parent
p=root/'finalized_study'
x=pd.read_csv(p/'finalized_participant_data.csv')
dictionary=pd.read_csv(p/'fi_item_dictionary.csv')
items=dictionary['item'].tolist()
assert len(items)==49 and len(set(items))==49
assert dictionary.loc[dictionary.source.str.contains('DPQ030|DPQ040'), 'item'].tolist()==['fi_sleep','fi_fatigue']
components={
 'fi_all':items,
 'fi_no_sleep':[v for v in items if v!='fi_sleep'],
 'fi_no_fatigue':[v for v in items if v!='fi_fatigue'],
 'fi_neither':[v for v in items if v not in ('fi_sleep','fi_fatigue')],
}
score_checks={}
for name,cols in components.items():
    assert x[cols].stack().dropna().between(0,1).all()
    n=x[cols].notna().sum(axis=1)
    mathceil=int(np.ceil(.8*len(cols)))
    calc=(x[cols].sum(axis=1,min_count=1)/n).where(n>=mathceil)
    stored=x[name]
    assert np.array_equal(n.to_numpy(),x[name+'_n'].to_numpy())
    assert np.array_equal(calc.isna(),stored.isna())
    delta=(calc-stored).abs().max()
    assert delta<1e-12
    score_checks[name]={'components':len(cols),'minimum_observed':mathceil,'max_abs_difference':float(delta)}
raw=[]
for year,suffix in [(2011,'G'),(2013,'H')]:
    dem,_=pyreadstat.read_xport(str(p/'source_files'/f'DEMO_{suffix}.XPT'))
    slq,_=pyreadstat.read_xport(str(p/'source_files'/f'SLQ_{suffix}.XPT'))
    join=dem[['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']].merge(slq[['SEQN','SLQ050']],on='SEQN',how='left',validate='one_to_one')
    join['cycle']=f'{year}-{year+1}'
    raw.append(join)
raw=pd.concat(raw,ignore_index=True)
z=x.merge(raw,on=['cycle','SEQN'],suffixes=('_saved','_raw'),validate='one_to_one')
assert len(z)==len(x)==3632
for name in ['RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR','SLQ050']:
    assert np.allclose(z[name+'_saved'],z[name+'_raw'],equal_nan=True),name
exposure=z.SLQ050_raw.map({1.:1.,2.:0.})
assert np.allclose(exposure,z.sleep_told_doctor,equal_nan=True)
core=(z.WTMEC2YR_raw>0)&z.sleep_told_doctor.notna()
for name in components:core &= z[name].notna()
assert int(core.sum())==3336
design=pd.read_csv(p/'survey_full_design_frame.csv')
assert np.allclose(design.wt4,design.WTMEC2YR/2,equal_nan=True)
assert len(design)==int((raw.WTMEC2YR>0).sum())
assert not design.duplicated(['cycle','SEQN']).any()
psus=design[['stratum','SDMVPSU']].drop_duplicates();n_strata=psus.stratum.nunique();df=len(psus)-n_strata
result={'score_checks':score_checks,'eligible_interview':len(x),'positive_MEC':int((z.WTMEC2YR_raw>0).sum()),'common':int(core.sum()),'exposed':int(z.loc[core,'sleep_told_doctor'].sum()),'full_design_n':len(design),'strata':n_strata,'PSUs':len(psus),'design_df':df,'limitations':'checks recomputation from saved coded items and raw exposure/design, not independent recoding of all 49 source fields'}
(root/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
