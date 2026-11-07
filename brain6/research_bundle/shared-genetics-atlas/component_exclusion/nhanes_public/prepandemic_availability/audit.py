import io,json,hashlib
from pathlib import Path
import pandas as pd
import urllib.request,ssl
p=Path(__file__).resolve().parent
ctx=ssl.create_default_context(cafile='/opt/miniconda3/lib/python3.13/site-packages/certifi/cacert.pem')
cols={'P_DEMO':['SEQN','RIDAGEYR','WTMECPRP','SDMVSTRA','SDMVPSU'],'P_SLQ':['SEQN','SLD012'],'P_WHQ':['SEQN','WHD020','WHD050','WHQ060'],'P_PAQ':['SEQN','PAQ650','PAQ665'],'P_DPQ':['SEQN','DPQ040']}
files={};frames={}
for name,fields in cols.items():
 u=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/{name}.xpt'
 with urllib.request.urlopen(u,context=ctx,timeout=60) as h:data=h.read()
 d=pd.read_sas(io.BytesIO(data),format='xport');assert set(fields).issubset(d.columns) and d.SEQN.is_unique
 for col in fields:
  if pd.api.types.is_numeric_dtype(d[col]):d.loc[d[col].notna() & (d[col].abs()<1e-60),col]=0.0
 frames[name]=d[fields];files[name]={'url':u,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data),'rows':len(d),'columns_used':fields}
x=frames['P_DEMO']
for key in ('P_SLQ','P_WHQ','P_PAQ','P_DPQ'):x=x.merge(frames[key],on='SEQN',how='left',validate='one_to_one')
x=x[x.RIDAGEYR>=60].copy();goodweight=x.WHD020.between(1,7776)&x.WHD050.between(1,7776);loss=goodweight&((x.WHD050-x.WHD020)/x.WHD050>=.05);absolute10=goodweight&(x.WHD050-x.WHD020>=10);intention=x.WHQ060.isin([1,2]);act=x.PAQ650.isin([1,2])&x.PAQ665.isin([1,2]);energy=x.DPQ040.isin([0,1,2,3]);sleep=x.SLD012.between(2,14);counts={'age_60_plus':len(x),'weight_pair_usable':int(goodweight.sum()),'recreational_activity_pair_usable':int(act.sum()),'energy_item_usable':int(energy.sum()),'weekday_sleep_usable':int(sleep.sum()),'at_least_5pct_loss':int(loss.sum()),'5pct_loss_and_10lb_or_more':int((loss&absolute10).sum()),'5pct_loss_and_less_than_10lb':int((loss&~absolute10).sum()),'5pct_loss_with_intentionality_answer':int((loss&intention).sum()),'5pct_loss_without_intentionality_answer':int((loss&~intention).sum()),'5pct_loss_unintentional_answer':int((loss&(x.WHQ060==2)).sum()),'5pct_loss_intentional_answer':int((loss&(x.WHQ060==1)).sum())}
assert counts['5pct_loss_with_intentionality_answer']+counts['5pct_loss_without_intentionality_answer']==counts['at_least_5pct_loss']
(p/'receipt.json').write_text(json.dumps({'files':files,'counts':counts,'pandas_version':pd.__version__,'sas_zero_handling':'Finite near-zero values smaller than 1e-60 were normalized to zero after SAS XPT decoding; CDC DPQ040 code 0 otherwise appeared approximately 5.4e-79.','limitation':'Unweighted descriptive counts; input bytes held in memory; no person-level data saved.'},indent=2)+'\n');print(json.dumps(counts,indent=2))
