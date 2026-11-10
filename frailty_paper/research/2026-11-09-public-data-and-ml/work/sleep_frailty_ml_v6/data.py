"""Locus-aware measured outcomes, features and declared difficult pairs."""
import gzip,collections,bisect
from common import *
SPLIT={'train':0,'validation':1,'calibration':2}
def load():
    audit=json.loads((OUT/'measured_data_audit.json').read_text())
    rows=list(csv.DictReader((DATA/'measured_observations.csv').open()))
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    cells=sorted({r['cell_group'] for r in rows});studies=sorted({r['study'] for r in rows});cell_id={s:i for i,s in enumerate(cells)};study_id={s:i for i,s in enumerate(studies)}
    n=len(rows);idx=np.array([int(r['bank_id']) for r in rows],dtype='int32');y=np.array([int(r['label']) for r in rows],dtype='float32');split=np.array([SPLIT[r['split']] for r in rows],dtype='int8');cell=np.array([cell_id[r['cell_group']] for r in rows],dtype='int32');study=np.array([study_id[r['study']] for r in rows],dtype='int32');effect=np.log1p(np.array([float(r['absolute_effect']) for r in rows],dtype='float32'));active=np.array([int(r['active_screen']) if r['active_screen'] else -1 for r in rows],dtype='int8');groups=np.array([r['locus_group'] for r in rows]);context=np.array([r['study']+'|'+r['cell_group'] for r in rows]);activity=np.array([float(r['activity_mean']) if r['activity_mean'] else np.nan for r in rows],dtype='float32')
    return {'rows':rows,'bank':bank,'bank_index':idx,'y':y,'split':split,'cell':cell,'cell_names':cells,'study':study,'study_names':studies,'effect':effect,'active':active,'groups':groups,'context':context,'activity':activity,'n':n,'audit':audit}
def hard_pairs(d):
    rows=d['rows'];bank=d['bank'];near=collections.defaultdict(list);across=collections.defaultdict(list)
    for i,r in enumerate(rows):
        bi=int(r['bank_id']);b=bank[bi];p=int(b['position_hg38']);a=float(r['absolute_effect']);positive=int(r['label'])==1 and a>=.2
        if r['study']=='Kircher2019':
            good_tags=int(r['tags'])>=20;negative=float(r['pvalue'])>=.1 and a<=.05 and good_tags;positive=positive and good_tags
        else:negative=int(r['label'])==0 and int(r['active_screen'])==1 and a<=.05
        if not positive and not negative:continue
        near[(r['study'],r['experiment'],r['cell'],r['locus_group'])].append((p,i,positive))
        across[(r['study'],b['variant_GRCh38'])].append((i,positive))
    pairs=set();kinds=collections.Counter()
    for key,items in near.items():
        neg=sorted((p,i) for p,i,positive in items if not positive);positions=[p for p,i in neg]
        for p,i,positive in items:
            if not positive:continue
            lo=bisect.bisect_left(positions,p-25);hi=bisect.bisect_right(positions,p+25);candidates=sorted(neg[lo:hi],key=lambda t:(abs(t[0]-p),t[1]))[:4]
            for q,j in candidates:pairs.add((i,j));kinds['nearby_same_experiment']+=1
    for items in across.values():
        for i,positive in items:
            if not positive:continue
            for j,p in items:
                if not p and rows[i]['cell']!=rows[j]['cell']:pairs.add((i,j));kinds['same_variant_across_cells']+=1
    p=np.array(sorted(pairs),dtype='int32').reshape(-1,2)
    assert len(p)>0 and np.all(d['split'][p[:,0]]==d['split'][p[:,1]]) and np.all(d['y'][p[:,0]]==1) and np.all(d['y'][p[:,1]]==0)
    return p,dict(kinds)
def columns(receipt,view):
    slices=receipt['slices']
    keys={'local':['local_hyena','local_pattern'],'biological':['local_hyena','full_hyena','local_pattern','local_motif'],'pattern_motif':['local_pattern','local_motif']}[view]
    return np.concatenate([np.arange(*slices[k],dtype='int32') for k in keys])
def statistics(bank_x,ids,cols):
    # Unique sequence rows only: repeated cell measurements cannot change
    # sequence normalization. Population variances, fitted on train only.
    ids=np.unique(ids);total=np.zeros(len(cols),dtype='float64');squares=np.zeros_like(total)
    for start in range(0,len(ids),1024):
        x=np.asarray(bank_x[np.ix_(ids[start:start+1024],cols)],dtype='float64');total+=x.sum(0);squares+=(x*x).sum(0)
    mean=total/len(ids);sd=np.sqrt(np.maximum(squares/len(ids)-mean*mean,0));sd[sd<1e-5]=1
    return mean.astype('float32'),sd.astype('float32')
def feature_rows(bank_x,d,ids,cols,mean,sd):
    x=np.asarray(bank_x[np.ix_(d['bank_index'][ids],cols)],dtype='float32');x=np.clip((x-mean)/sd,-20,20)
    return np.c_[x,np.eye(len(d['cell_names']),dtype='float32')[d['cell'][ids]]]
def weights_and_targets(d,tr):
    keys=np.unique(d['context'][tr]);row_weight=np.ones(d['n'],dtype='float32');pw=np.ones(d['n'],dtype='float32');target=np.zeros(d['n'],dtype='float32');metadata={}
    for key in keys:
        all_ids=np.flatnonzero(d['context']==key);ids=tr[d['context'][tr]==key];n=len(ids);positive=float(d['y'][ids].sum());weight=np.sqrt(len(tr)/(len(keys)*n));positive_weight=np.sqrt((n-positive)/max(positive,1));mu=float(d['effect'][ids].mean());sd=max(float(d['effect'][ids].std()),1e-4);row_weight[all_ids]=weight;pw[all_ids]=positive_weight;target[all_ids]=np.clip((d['effect'][all_ids]-mu)/sd,-5,5);metadata[str(key)]={'n':n,'positive':int(positive),'row_weight':weight,'positive_weight':positive_weight,'effect_mean':mu,'effect_sd':sd}
    row_weight=np.clip(row_weight,.2,5);row_weight/=row_weight[tr].mean()
    return row_weight,pw,target,metadata
