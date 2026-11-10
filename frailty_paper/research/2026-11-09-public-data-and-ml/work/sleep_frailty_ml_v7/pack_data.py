"""Compact numerical training inputs and measured contrast pairs, no held-out outcomes."""
import gzip,collections,bisect
from common import *
ROLE={'train':0,'validation':1,'calibration':2}
def main():
    assert not (OUT/'model_search_started.json').exists();audit=json.loads((OUT/'development_data_audit.json').read_text())
    assert sha(DATA/'measured_magnitudes.csv')==audit['outputs']['measured_magnitudes.csv']
    rows=list(csv.DictReader((DATA/'measured_magnitudes.csv').open()));line_alias={'NGN2_induced_neurons':'WTC11_NGN2'};lines=sorted({line_alias.get(r['cell_line'],r['cell_line']) for r in rows});li={v:i for i,v in enumerate(lines)};cells=sorted({r['cell_group'] for r in rows});studies=sorted({r['study'] for r in rows});experiments=sorted({r['study']+'|'+r['experiment']+'|'+r['cell_line'] for r in rows});loci=sorted({r['locus_group'] for r in rows});contexts=sorted({r['study']+'|'+r['cell_group']+'|'+r['effect_unit'] for r in rows});mapping=lambda names:{s:i for i,s in enumerate(names)};ci,si,ei,gi,xi=map(mapping,[cells,studies,experiments,loci,contexts])
    d={'bank_index':np.array([int(r['bank_id']) for r in rows],dtype='int32'),'effect':np.array([float(r['absolute_effect']) for r in rows],dtype='float32'),'unit':np.array([r['effect_unit']!='absolute_log2_effect' for r in rows],dtype='int8'),'split':np.array([ROLE[r['split']] for r in rows],dtype='int8'),'cell':np.array([ci[r['cell_group']] for r in rows],dtype='int16'),'study':np.array([si[r['study']] for r in rows],dtype='int8'),'experiment':np.array([ei[r['study']+'|'+r['experiment']+'|'+r['cell_line']] for r in rows],dtype='int16'),'group':np.array([gi[r['locus_group']] for r in rows],dtype='int32'),'context':np.array([xi[r['study']+'|'+r['cell_group']+'|'+r['effect_unit']] for r in rows],dtype='int16')}
    d['line']=np.array([li[line_alias.get(r['cell_line'],r['cell_line'])] for r in rows],dtype='int16');d['cell_line_names']=np.array(lines)
    assert np.isfinite(d['effect']).all() and np.all(d['effect']>=0)
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    n=len(bank);position=np.array([int(r['position_hg38']) for r in bank]);p51=np.empty(n,dtype='int16');ref=np.lib.format.open_memmap(DATA/'tokens_ref51.npy',mode='w+',dtype='uint8',shape=(n,51));alt=np.lib.format.open_memmap(DATA/'tokens_alt51.npy',mode='w+',dtype='uint8',shape=(n,51));lut=np.zeros(256,dtype='uint8');lut[[ord(c) for c in 'ACGT']]=[7,8,9,10]
    for i,r in enumerate(bank):
        p=int(r['changed_position']);start=max(0,min(p-25,len(r['reference'])-51));p51[i]=p-start;ref[i]=lut[np.frombuffer(r['reference'][start:start+51].encode(),dtype='uint8')];alt[i]=lut[np.frombuffer(r['alternate'][start:start+51].encode(),dtype='uint8')]
    ref.flush();alt.flush();assert np.all(np.sum(ref!=alt,axis=1)==1)
    near=collections.defaultdict(list);same=collections.defaultdict(list)
    for i in range(len(rows)):
        b=d['bank_index'][i];near[(int(d['experiment'][i]),int(d['group'][i]),int(d['unit'][i]))].append(i)
        if d['unit'][i]==0:same[(int(d['study'][i]),int(b))].append(i)
    local=set();cross=set()
    for ids in near.values():
        ids.sort(key=lambda i:position[d['bank_index'][i]]);pp=[position[d['bank_index'][i]] for i in ids]
        for i,p in zip(ids,pp):
            lo=bisect.bisect_left(pp,p-25);hi=bisect.bisect_right(pp,p+25);candidates=[]
            for j in ids[lo:hi]:
                if i==j or d['effect'][i]-d['effect'][j]<.1 or d['effect'][i]<.2:continue
                candidates.append(j)
            for j in sorted(candidates,key=lambda j:(abs(position[d['bank_index'][j]]-p),j))[:3]:local.add((i,j))
    for ids in same.values():
        for i in ids:
            for j in ids:
                if d['cell'][i]==d['cell'][j] or d['effect'][i]-d['effect'][j]<.05:continue
                cross.add((i,j))
    local=np.array(sorted(local),dtype='int32').reshape(-1,2);cross=np.array(sorted(cross),dtype='int32').reshape(-1,2)
    for pair in [local,cross]:
        assert np.all(d['split'][pair[:,0]]==d['split'][pair[:,1]]) and np.all(d['group'][pair[:,0]]==d['group'][pair[:,1]])
    d.update({'position51':p51,'local_pairs':local,'cross_cell_pairs':cross,'cell_names':np.array(cells),'study_names':np.array(studies),'experiment_names':np.array(experiments),'context_names':np.array(contexts),'locus_names':np.array(loci)})
    np.savez_compressed(DATA/'training_arrays.npz',**d);counts={}
    for name,code in ROLE.items():counts[name]={'observations':int(np.sum(d['split']==code)),'primary_log2_observations':int(np.sum((d['split']==code)&(d['unit']==0))),'local_continuous_contrast_pairs':int(np.sum(d['split'][local[:,0]]==code)),'cross_cell_continuous_contrast_pairs':int(np.sum(d['split'][cross[:,0]]==code)),'locus_groups':int(len(np.unique(d['group'][d['split']==code])))}
    write_json(OUT/'compact_data_and_pair_receipt.json',{'observations':len(rows),'bank_rows':n,'inventory':counts,'cell_groups':cells,'cell_line_names':lines,'line_aliases':line_alias,'studies':studies,'continuous_targets_not_detection_labels':True,'local_pair_rule':'Same experiment/cell/locus/unit;within25bp;larger magnitude>=.2 and magnitude gap>=.1;up to3nearest smaller measured effects. Coefficient-unit pairs stay in the separate auxiliary task.','cross_cell_pair_rule':'Same study and exact physical paired-DNA bank row, different cell group, measured magnitude gap>=.05log2;actual difference is regression target, not a binary call.','outputs':{name:sha(DATA/name) for name in ['training_arrays.npz','tokens_ref51.npy','tokens_alt51.npy']},'observed_development_data_audit_sha256':sha(OUT/'development_data_audit.json'),'reserved_effects_read_or_used':False});print('COMPACT EFFECT DATA AND CONTINUOUS PAIRS READY',json.dumps(counts),flush=True)
if __name__=='__main__':main()
