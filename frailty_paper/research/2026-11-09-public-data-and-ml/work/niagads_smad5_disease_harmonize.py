#!/usr/bin/env python3
"""Align open brain SMAD5 eQTL rows to existing local disease GWAS caches."""
import csv, gzip, hashlib, json, math, sys
from pathlib import Path
from pyliftover import LiftOver

root = Path(__file__).resolve().parent.parent
out = root / 'outputs'
assert len(sys.argv) <= 2 and (len(sys.argv) == 1 or sys.argv[1] == '--tad')
suffix = '_tad' if len(sys.argv) == 2 else ''
qfile = out / f'stage6_niagads_smad5_exc_region{suffix}.tsv.gz'
lo = LiftOver(str(root / 'work/hg38ToHg19.over.chain.gz'))
with gzip.open(qfile, 'rt') as f: qtl = list(csv.DictReader(f, delimiter='\t'))
variants = {}
for r in qtl:
    c, pos, ref, alt = r['variant_id'].split(':')
    if len(ref) != 1 or len(alt) != 1 or {ref, alt} in ({'A', 'T'}, {'C', 'G'}): continue
    af = float(r['non_ref_af'])
    if not .01 <= af <= .99: continue
    m = lo.convert_coordinate('chr5', int(pos)-1)
    if len(m) != 1 or m[0][0] != 'chr5' or m[0][2] != '+': continue
    key = (int(m[0][1])+1, frozenset((ref, alt)))
    if key in variants: raise ValueError(f'duplicate lifted variant {key}')
    variants[key] = r
report = {}
for trait in ['Mak_FI', 'Mak_FP', 'FinnGen_insomnia']:
    source = root / f'work/stage3_region_cache/{trait}.regional.tsv.gz'
    matched, duplicates = {}, set()
    with gzip.open(source, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if trait == 'FinnGen_insomnia':
                if r['#chrom'] != '5': continue
                pos = int(r['pos']); ref, alt = r['ref'], r['alt']
                m = lo.convert_coordinate('chr5', pos-1)
                if len(m) != 1 or m[0][0] != 'chr5' or m[0][2] != '+': continue
                key = (int(m[0][1])+1, frozenset((ref,alt)))
                effect, af, beta, se = alt, r['af_alt'], r['beta'], r['sebeta']
            else:
                if r['CHR'] != '5': continue
                key = (int(r['POS']), frozenset((r['A1'],r['A2'])))
                effect, af, beta, se = r['A1'], r['AF1'], r['BETA'], r['SE']
            if key not in variants: continue
            q = variants[key]; chrom, pos38, qref, qalt = q['variant_id'].split(':')
            if trait == 'FinnGen_insomnia' and (pos != int(pos38) or ref != qref or alt != qalt): continue
            try: af,beta,se = float(af),float(beta),float(se)
            except ValueError: continue
            if not all(math.isfinite(v) for v in (af,beta,se)) or se <= 0: continue
            if effect != qalt: af, beta = 1-af, -beta
            if not .01 <= af <= .99: continue
            row = {'variant_GRCh38': q['variant_id'], 'gwas_pos37': key[0],
                   'qtl_alt_af': float(q['non_ref_af']), 'qtl_beta_alt': float(q['beta_non_ref']),
                   'qtl_se': float(q['beta_se_non_ref']), 'gwas_af_alt': af,
                   'gwas_beta_alt': beta, 'gwas_se': se}
            if key in matched: duplicates.add(key)
            else: matched[key] = row
    for key in duplicates: matched.pop(key, None)
    rows = list(matched.values())
    path = out / f'stage6_niagads_smad5_exc_{trait}_aligned{suffix}.csv'
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    deltas = sorted(abs(r['qtl_alt_af']-r['gwas_af_alt']) for r in rows)
    report[trait] = {'matched':len(rows), 'frequency_difference_gt_0_1':sum(d>.1 for d in deltas),
                     'max_abs_af_difference':max(deltas), 'candidate_present':any(r['variant_GRCh38']=='chr5:136202653:A:C' for r in rows),
                     'duplicates_excluded':len(duplicates),'output_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
receipt = {'qtl_rows':len(qtl), 'nonpalindromic_qtl_liftover_eligible':len(variants), 'traits':report,
           'input_qtl_sha256':hashlib.sha256(qfile.read_bytes()).hexdigest(),
           'limitations':'Postselected window; source GWAS caches may be narrower; no clinical score validation.'}
(out / f'stage6_niagads_smad5_harmonization{suffix}.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
