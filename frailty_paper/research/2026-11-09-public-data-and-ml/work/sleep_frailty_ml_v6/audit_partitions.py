"""Independent complete split checks and difficult-pair inventory."""
from data import *
def main():
    d=load();rows=d['rows'];bank=d['bank'];a=d['audit']
    for name,h in a['outputs'].items():assert sha(DATA/name)==h
    bybank=collections.defaultdict(set);bylocus=collections.defaultdict(set);byvariant=collections.defaultdict(set)
    for r in rows:
        i=int(r['bank_id']);b=bank[i];assert r['locus_group']==b['locus_group'] and r['split']==b['split']
        bybank[i].add(r['split']);bylocus[r['locus_group']].add(r['split']);byvariant[b['variant_GRCh38']].add(r['split'])
    assert all(len(v)==1 for v in bybank.values()) and all(len(v)==1 for v in bylocus.values()) and all(len(v)==1 for v in byvariant.values())
    chroms=collections.defaultdict(list)
    for b in bank:chroms[b['chromosome']].append((int(b['position_hg38']),b['locus_group']))
    for items in chroms.values():
        items.sort()
        for (p,g),(q,h) in zip(items,items[1:]):
            if q-p<=2000:assert g==h
    pairs,kinds=hard_pairs(d);np.savez_compressed(DATA/'hard_pairs.npz',pairs=pairs)
    inventory={s:{'observations':int(np.sum(d['split']==i)),'called_effects':int(d['y'][d['split']==i].sum()),'hard_pairs':int(np.sum(d['split'][pairs[:,0]]==i))} for s,i in SPLIT.items()}
    train_cells=set(d['cell'][d['split']==0].tolist());out={s:sorted({d['cell_names'][i] for i in d['cell'][d['split']==j] if i not in train_cells}) for s,j in SPLIT.items() if j}
    write_json(OUT/'partition_audit.json',{'passed':True,'all_observations_checked':len(rows),'all_sequence_rows_checked':len(bank),'variant_locus_and_2kb_neighbor_split_checks':True,'inventory':inventory,'hard_pair_kinds':kinds,'hard_pair_definition':'Observed detected effect with abs(source effect)>=0.2 versus active non-detection with abs(effect)<=0.05; for Kircher require >=20tags, negative p>=0.1. Match same experiment within25bp (max4nearest negatives), or same exact variant across cells. Source effect units differ; these are training strata, not a common biological effect threshold.','hard_pairs_sha256':sha(DATA/'hard_pairs.npz'),'cell_groups_unseen_in_training':out,'final_study_outcomes_opened':False})
    print(json.dumps({'inventory':inventory,'unseen_cells':out,'hard_pair_kinds':kinds}),flush=True)
if __name__=='__main__':main()
