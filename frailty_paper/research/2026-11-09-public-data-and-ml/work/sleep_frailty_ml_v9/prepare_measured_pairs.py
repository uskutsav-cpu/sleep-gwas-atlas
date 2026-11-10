"""Actual nearby and cross-cell comparisons in the augmented development set."""
from v9_common import *


def main():
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in ['bank_index','cell','line','study','experiment','group','split']}
    with np.load(DATA/'partition_bank_metadata.npz') as z:
        position=np.zeros(int(d['bank_index'].max())+1,dtype='int64');position[z['bank_index']]=z['position']
    order=np.lexsort((d['line'],d['cell'],d['bank_index'],d['study']));cross=[];i=0;n=len(order)
    while i<n:
        j=i+1
        while j<n and d['study'][order[j]]==d['study'][order[i]] and d['bank_index'][order[j]]==d['bank_index'][order[i]]:j+=1
        ids=order[i:j]
        for a in range(len(ids)):
            for b in range(a+1,len(ids)):
                if d['cell'][ids[a]]!=d['cell'][ids[b]]:cross.append((ids[a],ids[b]))
        i=j
    order=np.lexsort((position[d['bank_index']],d['group'],d['experiment']));a,b=order[:-1],order[1:];gap=position[d['bank_index'][b]]-position[d['bank_index'][a]]
    keep=(d['experiment'][a]==d['experiment'][b])&(d['group'][a]==d['group'][b])&(gap>=0)&(gap<=25)&(d['bank_index'][a]!=d['bank_index'][b])
    local=np.c_[a[keep],b[keep]].astype('int32');cross=np.array(cross,dtype='int32').reshape(-1,2)
    for p in [local,cross]:assert np.all(d['split'][p[:,0]]==d['split'][p[:,1]]) and np.all(d['group'][p[:,0]]==d['group'][p[:,1]])
    np.savez_compressed(DATA/'measured_pairs.npz',local=local,cross=cross)
    write_json(OUT/'measured_pair_audit.json',{'local_pairs':len(local),'cross_cell_pairs':len(cross),'genomic_distance_bp':25,'all_pairs_within_same_study_and_original_experiment_for_local_pairs':True,'no_significance_or_effect_gap_selection':True,'all_pairs_in_one_partition_and_genomic_component':True,'data_sha256':sha(DATA/'training_arrays.npz'),'output_sha256':sha(DATA/'measured_pairs.npz')})
    print('AUGMENTED MEASURED PAIRS READY',len(local),len(cross),flush=True)


if __name__=='__main__':main()
