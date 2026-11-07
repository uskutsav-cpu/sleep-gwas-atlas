#!/usr/bin/env python3
import csv, gzip, hashlib, heapq, json, math
from itertools import zip_longest
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = Path('/Volumes/Brain6Work')
PATHS = {'hfrs': BASE/'ukb_hfrs_public_20261028.sumstats.gz',
         'questionnaire': BASE/'ukb-hfrs-comparison-20261028/questionnaire.sumstats.gz',
         'registry': BASE/'ukb-hfrs-comparison-20261028/registry.sumstats.gz'}
EXPECTED = {'hfrs':'477f245c86aeea045d1226fc827a409b6d875615fdd800597c2611b96d34e78a',
            'questionnaire':'4ff71a0d46ce0357fbc4b9fc6471aaa69b2ec1e5891cc043b6614862adbbddd5',
            'registry':'6fb908da20f8dacf652320028583b61e558897039c18fb01c5fe8f9e33262e99'}
N = 1020621
M = 2*N
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()
for name,path in PATHS.items(): assert sha(path) == EXPECTED[name]

values = np.ones(M, dtype=np.float64)
heaps = {'questionnaire': [], 'registry': []}
hits = []
seen = set()
counts = {'invalid_pair_slots': 0, 'allele_mismatch_rows': 0}
files = [gzip.open(PATHS[k], 'rt', newline='') for k in ('hfrs','questionnaire','registry')]
readers = [csv.DictReader(f, delimiter='\t') for f in files]
try:
    for i, triple in enumerate(zip_longest(*readers)):
        assert i < N and all(x is not None for x in triple)
        h, q, r = triple
        assert h['SNP'] == q['SNP'] == r['SNP'] and h['SNP'] not in seen
        seen.add(h['SNP'])
        ha = {h['A1'],h['A2']}
        assert len(ha)==2 and ha <= set('ACGT') and ha not in ({'A','T'},{'C','G'})
        if any({x['A1'], x['A2']} != ha for x in (q,r)):
            counts['allele_mismatch_rows'] += 1
            counts['invalid_pair_slots'] += 2
            continue
        for j, (label, x) in enumerate((('questionnaire',q),('registry',r))):
            try:
                zh,zs = float(h['Z']),float(x['Z'])
                if not math.isfinite(zh) or not math.isfinite(zs): raise ValueError
            except ValueError:
                counts['invalid_pair_slots'] += 1
                continue
            ph = math.erfc(max(0,abs(zh)-1e-9)/math.sqrt(2))
            ps = math.erfc(max(0,abs(zs)-0.0005)/math.sqrt(2))
            p = max(ph,ps)
            values[2*i+j] = p
            aligned_z = zs if h['A1'] == x['A1'] else -zs
            row = {'source':label,'SNP':h['SNP'],'hfrs_A1':h['A1'],'hfrs_A2':h['A2'],
                   'hfrs_Z':zh,'insomnia_Z_on_hfrs_A1':aligned_z,'hfrs_p_guarded':ph,
                   'insomnia_p_guarded':ps,'max_p':p,'bonferroni_two_pair':min(1,p*M)}
            entry=(-p, -i, row)
            if len(heaps[label]) < 100: heapq.heappush(heaps[label],entry)
            elif entry > heaps[label][0]: heapq.heapreplace(heaps[label],entry)
            if p <= .05/M: hits.append(row)
    assert i+1 == N
finally:
    for f in files: f.close()
ordered=np.sort(values)
ranks=np.arange(1,M+1,dtype=np.float64)
bh_raw=ordered*M/ranks
harmonic=float(np.sum(1/ranks))
bh_count=int(np.max(np.flatnonzero(bh_raw<=.05))+1) if np.any(bh_raw<=.05) else 0
by_count=int(np.max(np.flatnonzero(bh_raw*harmonic<=.05))+1) if np.any(bh_raw*harmonic<=.05) else 0
top=sorted([e[2] for heap in heaps.values() for e in heap],key=lambda x:(x['source'],x['max_p'],x['SNP']))
fields=list(top[0])
for name, records in [('top100_per_source.tsv',top),('bonferroni_hits.tsv',hits)]:
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(records)
summary={'input_sha256':EXPECTED,'rows':N,'family_slots':M,'counts':counts,
         'bonferroni_threshold':.05/M,'bonferroni_hits':len(hits),'bh_rejections':bh_count,
         'by_rejections':by_count,'smallest_bh_adjusted':min(1,float(np.min(bh_raw))),
         'smallest_by_adjusted':min(1,float(np.min(bh_raw))*harmonic),
         'minimum_per_source':{k:min((e[2] for e in heap),key=lambda x:x['max_p']) for k,heap in heaps.items()},
         'output_sha256':{n:sha(ROOT/n) for n in ('top100_per_source.tsv','bonferroni_hits.tsv')}}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
