"""All matched cells and consecutive measured nearby variants, no hit filter."""
from v8_common import *


def main():
    with np.load(DATA/'measured_training_arrays.npz') as z:
        d = {k: z[k] for k in ['bank_index', 'cell', 'study', 'experiment', 'group', 'split']}
    n = len(d['cell']); position = np.zeros(int(d['bank_index'].max())+1, dtype='int64')
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        if int(r['bank_id']) < len(position):
            position[int(r['bank_id'])] = int(r['position_hg38'])
    order = np.lexsort((d['cell'], d['bank_index'], d['study']))
    cells = []; i = 0
    while i < n:
        j = i+1
        while j < n and d['bank_index'][order[j]] == d['bank_index'][order[i]] and d['study'][order[j]] == d['study'][order[i]]:
            j += 1
        group = order[i:j]
        for a in range(len(group)):
            for b in range(a+1, len(group)):
                if d['cell'][group[a]] != d['cell'][group[b]]:
                    cells.append((group[a], group[b]))
        i = j
    order = np.lexsort((position[d['bank_index']], d['group'], d['experiment']))
    a, b = order[:-1], order[1:]
    gap = position[d['bank_index'][b]]-position[d['bank_index'][a]]
    keep = (d['experiment'][a] == d['experiment'][b])&(d['group'][a] == d['group'][b])&(gap >= 0)&(gap <= 25)&(d['bank_index'][a] != d['bank_index'][b])
    local = np.c_[a[keep], b[keep]].astype('int32')
    cross = np.array(cells, dtype='int32').reshape(-1, 2)
    for pairs in [local, cross]:
        assert np.all(d['split'][pairs[:, 0]] == d['split'][pairs[:, 1]])
        assert np.all(d['group'][pairs[:, 0]] == d['group'][pairs[:, 1]])
    np.savez_compressed(DATA/'measured_pair_arrays.npz', local=local, cross=cross)
    write_json(OUT/'measured_pair_receipt.json', {
        'local_pairs': len(local), 'matched_cell_pairs': len(cross),
        'source_call_or_effect_magnitude_filter': False,
        'local_rule': 'Consecutive measured variants within 25 bp in the same original experiment and genomic component; no significance or minimum effect-gap selection.',
        'cell_rule': 'All complete pairs of different cells for the same physical DNA pair and study.',
        'all_pairs_within_one_genomic_component_and_split': True,
        'output_sha256': sha(DATA/'measured_pair_arrays.npz'),
        'data_sha256': sha(DATA/'measured_training_arrays.npz'),
    })
    print('MEASURED PAIRS READY', len(local), len(cross), flush=True)


if __name__ == '__main__':
    main()
