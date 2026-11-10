"""Bound tensor shapes while retaining every real training observation."""
import numpy as np


def epoch_batches(ids,lengths,rng,batch_size):
    batches=[]
    for length in np.unique(lengths[ids]):
        selected=ids[lengths[ids]==length];order=rng.permutation(selected)
        batches.extend((int(length),order[i:i+batch_size]) for i in range(0,len(order),batch_size))
    order=rng.permutation(len(batches))
    return [batches[i] for i in order]


def pair_buckets(pairs,lengths,train):
    result={}
    for kind,values in pairs.items():
        p=values[np.all(train[values],axis=1)];assert np.all(lengths[p[:,0]]==lengths[p[:,1]])
        result[kind]={int(length):p[lengths[p[:,0]]==length] for length in np.unique(lengths[p[:,0]])}
    return result


def assemble(main,length,pools,rng,local_count,cross_count):
    chunks=[main];counts=[]
    for kind,count in [('local',local_count),('cross',cross_count)]:
        p=pools[kind].get(length);n=count if p is not None and len(p) else 0;counts.append(n)
        if n:chunks.append(p[rng.integers(len(p),size=n)].reshape(-1))
    return np.concatenate(chunks),counts
