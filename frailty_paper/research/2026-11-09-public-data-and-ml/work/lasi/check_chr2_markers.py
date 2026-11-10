#!/usr/bin/env python3
"""Scan complete public LASI-DAD LD table for two prespecified candidate pairs."""
import csv, gzip, hashlib, json, os

path = 'work/lasi/chr2_ld.csv.gz'
expected = 1397139462
if os.path.getsize(path) != expected:
    raise SystemExit('Incomplete chromosome 2 download')
targets = {
    '2:66523432': {'rsid': 'rs113851554', 'alleles': {'G', 'T'}, 'counted': 'T'},
    '2:66795102': {'rsid': 'rs78206187', 'alleles': {'A', 'G'}, 'counted': 'G'},
    '2:58695786': {'rsid': 'rs35881094', 'alleles': {'G', 'T'}, 'counted': 'G'},
    '2:58751388': {'rsid': 'rs17550821', 'alleles': {'A', 'G'}, 'counted': 'A'},
}
pairs = [('rs113851554', 'rs78206187'), ('rs35881094', 'rs17550821')]
presence = {v['rsid']: {} for v in targets.values()}
matches = {f'{a}/{b}': [] for a,b in pairs}
digest = hashlib.sha256()
with open(path,'rb') as f:
    for chunk in iter(lambda:f.read(1048576),b''): digest.update(chunk)
with gzip.open(path,'rt',newline='') as f:
    header = next(f)
    for line in f:
        if not any(key in line for key in targets): continue
        row = next(csv.reader([line])); found=[]
        for idx,afidx in ((5,7),(6,8)):
            identifier=row[idx]; parts=identifier.split(':'); key=':'.join(parts[:2])
            if key not in targets or set(parts[2:4]) != targets[key]['alleles']: continue
            rsid=targets[key]['rsid']; found.append(rsid)
            item=presence[rsid].setdefault(identifier,{'rows':0,'alternate_frequency':row[afidx]})
            item['rows']+=1
        for a,b in pairs:
            if a in found and b in found:
                flips = []
                for idx in (5, 6):
                    parts = row[idx].split(':'); target = targets[':'.join(parts[:2])]
                    flips.append(1 if target['counted'] == parts[3] else -1)
                matches[f'{a}/{b}'].append({
                    'source_row': row,
                    'counted_allele_haplotype_r_approx': float(row[3])**0.5 * int(row[11]) * flips[0] * flips[1]
                })
receipt={'bytes':expected,'sha256':digest.hexdigest(),'header':header.strip(),
         'provisional_grch38_targets':{k:{'rsid':v['rsid'],'alleles':sorted(v['alleles']),'counted':v['counted']} for k,v in targets.items()},
         'candidate_rows':presence,'candidate_pair_rows':matches,
         'limitations':'Public thresholded lookup has r2 >= 0.2 and 1 Mb distance filters; no score or outcome data. Coordinates and alleles require verification; absent records need not mean r2=0. Signed r is derived from rounded r2 and source haplotype Sign under biallelic allele flipping, not an individual genotype correlation.'}
with open('outputs/stage6_lasi_chr2_ld_receipt.json','w') as f:json.dump(receipt,f,indent=2)
print(json.dumps(receipt,indent=2))
