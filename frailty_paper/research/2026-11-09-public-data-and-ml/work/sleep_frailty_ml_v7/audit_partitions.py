"""Independently recheck saved partitions without opening final effect values."""
import gzip,bisect,collections
from common import *
def main():
    test=list(csv.DictReader((DATA/'reserved_physical_design.csv').open()));points=collections.defaultdict(list);cores=set();rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    for r in test:
        points[r['chromosome']].append(int(r['position_hg38']));p=int(r['changed_position'])
        for seq in [r['reference'],r['alternate']]:
            for j in range(max(0,p-149),min(p,len(seq)-150)+1):cores.add(seq[j:j+150]);cores.add(rc(seq[j:j+150]))
    points={k:sorted(set(v)) for k,v in points.items()};prefix={s[:16] for s in cores};roles={};chrom=collections.defaultdict(list);seen={};badclose=badcore=0;bankgroups=[];bankroles=[]
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:
        for number,r in enumerate(csv.DictReader(f)):
            assert number==int(r['bank_id']);group,role=r['locus_group'],r['split'];assert group not in roles or roles[group]==role;roles[group]=role;bankgroups.append(group);bankroles.append(role);c=r['chromosome'];p=int(r['position_hg38']);chrom[c].append((p,group,role));ps=points.get(c,[]);j=bisect.bisect_left(ps,p);close=any(0<=k<len(ps) and abs(ps[k]-p)<=2000 for k in [j-1,j]);badclose+=int(close);sub=int(r['changed_position']);a,b=r['reference'],r['alternate'];start=max(0,min(sub-75,len(a)-150));core=a[start:start+150];canonical=min(core,rc(core))
            if canonical in seen:assert seen[canonical]==(group,role),'Shared central context crossed a role/group'
            else:seen[canonical]=(group,role)
            match=False
            for s in [a,b]:
                for j in range(max(0,sub-149),min(sub,len(s)-150)+1):
                    if s[j:j+16] in prefix and s[j:j+150] in cores:match=True;break
            badcore+=int(match)
    adjacent=0
    for c,rows in chrom.items():
        rows.sort()
        for a,b in zip(rows[:-1],rows[1:]):
            if b[0]-a[0]<=2000:adjacent+=1;assert a[1:]==b[1:],'Nearby DNA crossed group/role'
    assert badclose==badcore==0
    with np.load(DATA/'training_arrays.npz') as z:
        b=z['bank_index'];g=z['group'];split=z['split'];locus_names=z['locus_names'];role={'train':0,'validation':1,'calibration':2};assert np.array_equal(split,np.array([role[bankroles[i]] for i in b]));assert np.array_equal(locus_names[g],np.array(bankgroups)[b]);paircounts={}
        for key in ['local_pairs','cross_cell_pairs']:
            pairs=z[key];assert np.all(split[pairs[:,0]]==split[pairs[:,1]]) and np.all(g[pairs[:,0]]==g[pairs[:,1]]);paircounts[key]=len(pairs)
    write_json(OUT/'independent_partition_audit.json',{'passed':True,'bank_rows':len(bankroles),'saved_observations':len(b),'nearby_pairs_checked':adjacent,'locus_groups':len(roles),'reserved_coordinate_overlaps_within_2kb':badclose,'reserved_exact_variant_containing_150bp_core_overlaps':badcore,'continuous_contrast_pairs_checked':paircounts,'input_sha256':{n:sha(DATA/n) for n in ['sequence_bank.csv.gz','training_arrays.npz','reserved_physical_design.csv']},'role_and_group_assignments_independently_reconstructed':True,'reserved_effect_values_read':False,'scope':'Physical/context and2kb-neighborhood partition audit. This does not prove ancestry-specific LD, participant or gene-wide independence.'});print('INDEPENDENT PARTITION AUDIT PASSED',len(b),flush=True)
if __name__=='__main__':main()
