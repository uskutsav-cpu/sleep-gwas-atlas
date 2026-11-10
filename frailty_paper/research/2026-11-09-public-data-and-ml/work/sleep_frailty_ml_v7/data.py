"""Training-only transformations, explicit effect units and real cell contexts."""
from common import *
from sequence_features import branches
def load():
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    d['n']=len(d['effect']);d['log_effect']=np.log1p(d['effect']).astype('float32');d['cell_names']=d['cell_names'].tolist();d['cell_line_names']=d['cell_line_names'].tolist();d['study_names']=d['study_names'].tolist();d['context_names']=d['context_names'].tolist();d['x']=np.load(DATA/'feature_bank.npy',mmap_mode='r')
    families={'induced_excitatory_neuron':'neural','neuronal_tumor':'neural','neural_progenitor':'neural','glioma':'neural','hepatocyte_tumor':'hepatic','kidney_derived':'renal','hematopoietic':'hematopoietic','colon_tumor':'intestinal'};family_names=sorted({families.get(c,'other') for c in d['cell_names']});lookup={c:i for i,c in enumerate(family_names)};d['family_names']=family_names;d['family_by_cell']=np.array([lookup[families.get(c,'other')] for c in d['cell_names']],dtype='int16');d['n_context_flags']=len(d['cell_names'])+len(d['cell_line_names'])+len(family_names)
    return d
def fit_transform(d,ids):
    unique=np.unique(d['bank_index'][ids]);c,v=branches();total=np.zeros(1145);sq=np.zeros(1145)
    for start in range(0,len(unique),2048):
        x=np.asarray(d['x'][unique[start:start+2048]],dtype='float64');total+=x.sum(0);sq+=(x*x).sum(0)
    mu=total/len(unique);sd=np.sqrt(np.maximum(sq/len(unique)-mu*mu,0));mu[v]=0;sd[v]=np.sqrt(np.maximum(sq[v]/len(unique),0));sd[sd<1e-5]=1
    scales=[];weights=np.zeros(d['n'],dtype='float32');metadata={}
    for unit in [0,1]:
        tr=ids[d['unit'][ids]==unit];scales.append(max(float(np.std(d['log_effect'][tr])),.05));keys=np.unique(d['context'][tr])
        for key in keys:
            these=tr[d['context'][tr]==key];w=min(len(ids)/(len(keys)*len(these)),20.);weights[d['context']==key]=w*(1 if unit==0 else .15);metadata[d['context_names'][key]]={'rows':len(these),'unit':unit,'training_weight':w*(1 if unit==0 else .15)}
    return {'mean':mu.astype('float32'),'sd':sd.astype('float32'),'effect_scales':np.array(scales,dtype='float32'),'context_columns':c,'variant_columns':v,'weights':weights,'fit_cells':np.unique(d['cell'][ids]),'fit_lines':np.unique(d['line'][ids]),'fit_families':np.unique(d['family_by_cell'][d['cell'][ids]]),'normalization_unique_DNA_rows':len(unique),'context_weight_metadata':metadata}
def flags(d,ids,t,no_cell=False):
    result=np.zeros((len(ids),d['n_context_flags']),dtype='float32')
    if no_cell:return result
    cell=d['cell'][ids];line=d['line'][ids];family=d['family_by_cell'][cell];row=np.arange(len(ids));nc=len(d['cell_names']);nl=len(d['cell_line_names'])
    use=np.isin(cell,t['fit_cells']);result[row[use],cell[use]]=1;use=np.isin(line,t['fit_lines']);result[row[use],nc+line[use]]=1;use=np.isin(family,t['fit_families']);result[row[use],nc+nl+family[use]]=1
    return result
def features(d,ids,t,no_cell=False,avi=None):
    raw=np.asarray(d['x'][d['bank_index'][ids]],dtype='float32');raw=np.clip((raw-t['mean'])/t['sd'],-12,12);c=np.c_[raw[:,t['context_columns']],flags(d,ids,t,no_cell)];v=raw[:,t['variant_columns']]
    if avi is not None:
        a=np.array(avi[d['bank_index'][ids]],copy=True);a[:,:2]=np.where(np.isfinite(a[:,:2]),a[:,:2],t['avi_median']);a[:,:2]=(a[:,:2]-t['avi_mean'])/t['avi_sd'];c=np.c_[c,np.clip(a,-12,12)]
    return c.astype('float32'),v.astype('float32')
def add_avi_transform(d,ids,t,avi):
    rows=avi[np.unique(d['bank_index'][ids])];known=rows[:,2]==1;assert known.sum()>100;t=dict(t);t['avi_median']=np.nanmedian(rows[known,:2],axis=0);filled=np.where(np.isfinite(rows[:,:2]),rows[:,:2],t['avi_median']);t['avi_mean']=filled.mean(0);t['avi_sd']=filled.std(0);t['avi_sd'][t['avi_sd']<1e-5]=1;return t
