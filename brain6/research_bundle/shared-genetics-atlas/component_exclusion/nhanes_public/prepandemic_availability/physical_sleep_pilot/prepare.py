import io,json,hashlib,ssl,urllib.request,sys
from pathlib import Path
import pandas as pd,numpy as np
p=Path(__file__).resolve().parent;receipt=json.loads((p.parent/'receipt.json').read_text());ctx=ssl.create_default_context(cafile='/opt/miniconda3/lib/python3.13/site-packages/certifi/cacert.pem');frames={}
for name,info in receipt['files'].items():
 with urllib.request.urlopen(info['url'],context=ctx,timeout=60) as f:raw=f.read()
 assert hashlib.sha256(raw).hexdigest()==info['sha256']
 d=pd.read_sas(io.BytesIO(raw),format='xport')[info['columns_used']];assert d.SEQN.is_unique
 for col in info['columns_used']:
  if pd.api.types.is_numeric_dtype(d[col]):d.loc[d[col].notna()&(d[col].abs()<1e-60),col]=0.0
 frames[name]=d
x=frames['P_DEMO']
for k in ('P_SLQ','P_WHQ','P_PAQ','P_DPQ'):x=x.merge(frames[k],on='SEQN',how='left',validate='one_to_one')
x=x[x.WTMECPRP>0].copy();valid=x.WHD020.between(1,7776)&x.WHD050.between(1,7776);loss=valid&((x.WHD050-x.WHD020)/x.WHD050>=.05);x['weightloss']=np.nan;x.loc[valid&~loss,'weightloss']=0;x.loc[loss&(x.WHQ060==1),'weightloss']=0;x.loc[loss&(x.WHQ060==2),'weightloss']=1;x['ambiguous']=(loss&~x.WHQ060.isin([1,2])).astype(int);x['inactivity']=np.nan;v=x.PAQ650.isin([1,2])&x.PAQ665.isin([1,2]);x.loc[v,'inactivity']=((x.loc[v,'PAQ650']==2)&(x.loc[v,'PAQ665']==2)).astype(int);x['energy']=np.nan;v=x.DPQ040.isin([0,1,2,3]);x.loc[v,'energy']=(x.loc[v,'DPQ040']>=2).astype(int);x['sleep']=x.SLD012.where(x.SLD012.between(2,14));x[['RIDAGEYR','WTMECPRP','SDMVSTRA','SDMVPSU','sleep','weightloss','ambiguous','inactivity','energy']].to_csv(sys.stdout,index=False,na_rep='')
