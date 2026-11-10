"""Pure NumPy/sklearn inference for frozen measured-data models."""
import pickle
from common import *
from network import forward
def load_member(record):
    path=Path(record['checkpoint']);assert sha(path)==record['checkpoint_sha256']
    if path.suffix=='.pkl':
        with path.open('rb') as f:p=pickle.load(f)
        return {'kind':'tree',**p}
    with np.load(path) as z:p={k:z[k] for k in z.files}
    recipe=json.loads(str(p['recipe_json']));return {'kind':'neural','parameters':p,'recipe':recipe,'columns':p['feature_columns'],'mean':p['feature_mean'],'sd':p['feature_sd'],'cell_names':p['cell_names'].tolist()}
def inputs(member,bank_x,bank_index,cells):
    cols=member['columns'];x=np.asarray(bank_x[np.ix_(bank_index,cols)],dtype='float32');x=np.clip((x-member['mean'])/member['sd'],-20,20)
    cell_names=member['cell_names'];mapping={s:i for i,s in enumerate(cell_names)};one=np.zeros((len(cells),len(cell_names)),dtype='float32')
    for i,c in enumerate(cells):
        if str(c) in mapping:one[i,mapping[str(c)]]=1
    return np.c_[x,one],np.mean(x*x,axis=1)
def predict_member(record,bank_x,bank_index,cells,batch_size=512):
    member=load_member(record);pred=[];distance=[]
    for start in range(0,len(bank_index),batch_size):
        end=min(start+batch_size,len(bank_index));x,dist=inputs(member,bank_x,bank_index[start:end],cells[start:end]);distance.extend(dist.tolist())
        if member['kind']=='tree':p=member['model'].predict_proba(x)[:,1]
        else:
            params=member['parameters'];z=forward(params,x,np.zeros(len(x),dtype='int32'),member['recipe']['architecture'],params['context_columns'],params['variant_columns'],backend='numpy')[0];p=1/(1+np.exp(-np.clip(z,-30,30)))
        pred.extend(p.tolist())
    p=np.array(pred);assert np.isfinite(p).all();return p,np.array(distance)
def predict_registry(registry,bank_x,bank_index,cells):
    families={};members={};distance=None
    for record in registry['records']:
        p,dist=predict_member(record,bank_x,bank_index,cells);families.setdefault(record['name'],[]).append(p);members[record['name']+'_'+str(record.get('seed','tree'))]=p
        if distance is None and record['name'] in registry['selected_family_weights']:distance=dist
    out={name:np.mean(p,axis=0) for name,p in families.items()};out['primary']=sum(weight*out[name] for name,weight in registry['selected_family_weights'].items())
    return out,members,distance
