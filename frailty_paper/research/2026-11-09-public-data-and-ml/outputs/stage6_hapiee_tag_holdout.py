"""Select array-position tags in non-Finnish 1000G EUR; evaluate in FIN.

This assesses reference-panel sensitivity, not cohort imputation accuracy.
"""
import csv, io, json, zipfile
from pathlib import Path
import numpy as np
import pysam

root = Path(__file__).resolve().parents[1]
outdir, work = root / 'outputs', root / 'work'
targets = {}
for file, c, p in [('unvalidated_15_variant_score_prototype.csv', 'chr', 'position_GRCh37_inferred'),
                    ('unvalidated_22_shared_direction_score.csv', 'chr_GRCh37_assumed', 'position_GRCh37_assumed')]:
    with (outdir / file).open() as fh:
        for row in csv.DictReader(fh):
            x = (str(row[c]), int(row[p]), frozenset((row['counted_allele'], row['other_allele'])))
            if row['variant'] in targets and targets[row['variant']] != x:
                raise ValueError('Inconsistent variant definitions')
            targets[row['variant']] = x
array = {str(i): set() for i in range(1, 23)}
with zipfile.ZipFile(work / 'gsa_v3_a1_hg19.zip') as archive:
    fh = io.TextIOWrapper(archive.open(archive.namelist()[0]), encoding='utf-8-sig', newline='')
    for line in fh:
        if line.strip() == '[Assay]':
            break
    for row in csv.DictReader(fh):
        if row['IlmnID'] == '[Controls]':
            break
        try:
            position = int(row['MapInfo'])
        except (ValueError, TypeError):
            continue
        snp = row['SNP']
        if row['Chr'] in array and len(snp) == 5 and snp[0] == '[' and snp[2] == '/' and snp[4] == ']' and set(snp[1:4:2]) <= set('ACGT') and row['Intensity_Only'] == '0':
            array[row['Chr']].add(position)
with (work / '1kg_phase3_panel.tsv').open() as fh:
    panel = list(csv.DictReader(fh, delimiter='\t'))
train_ids = {r['sample'] for r in panel if r['super_pop'] == 'EUR' and r['pop'] != 'FIN'}
test_ids = {r['sample'] for r in panel if r['pop'] == 'FIN'}
url = 'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

def dosage(record, samples):
    values = []
    for sample in samples:
        gt = record.samples[sample].get('GT')
        if gt is None or len(gt) != 2 or None in gt:
            return None
        values.append(sum(a == 1 for a in gt))
    values = np.asarray(values, dtype=float)
    return values if values.std() > 0 else None

rows = []
for chrom in sorted({v[0] for v in targets.values()}, key=int):
    vcf = pysam.VariantFile(url.format(chrom=chrom))
    train = sorted(train_ids & set(vcf.header.samples))
    test = sorted(test_ids & set(vcf.header.samples))
    for name, (ch, pos, alleles) in targets.items():
        if ch != chrom:
            continue
        candidates = [v for v in vcf.fetch(ch, pos - 1, pos) if v.pos == pos and v.id and name in v.id.split(';') and len(v.alleles) == 2 and set(v.alleles) == set(alleles)]
        result = dict(variant=name, chr=ch, pos37=pos, train_n=len(train), holdout_n=len(test), status='candidate_missing_or_nonunique', tags_evaluated=0, tag_id='', tag_pos37='', train_r2='', holdout_r2='')
        if len(candidates) == 1:
            target_train = dosage(candidates[0], train)
            target_test = dosage(candidates[0], test)
            result['status'] = 'candidate_incomplete_or_invariant_train'
            if target_train is not None:
                result['status'] = 'no_eligible_tag'
                best = (-1., None)
                for v in vcf.fetch(ch, max(0, pos - 100001), pos + 100000):
                    if v.pos not in array[ch] or len(v.alleles) != 2 or v.pos == pos:
                        continue
                    x = dosage(v, train)
                    if x is None:
                        continue
                    r2 = float(np.corrcoef(target_train, x)[0, 1] ** 2)
                    result['tags_evaluated'] += 1
                    if r2 > best[0]:
                        best = (r2, v)
                if best[1] is not None:
                    tag = best[1]
                    tag_test = dosage(tag, test)
                    result.update(status='holdout_invariant_or_incomplete', tag_id=tag.id or '.', tag_pos37=tag.pos, train_r2=best[0])
                    if target_test is not None and tag_test is not None:
                        result.update(status='evaluated', holdout_r2=float(np.corrcoef(target_test, tag_test)[0, 1] ** 2))
        rows.append(result)
        print(name, result['status'], result['train_r2'], result['holdout_r2'], flush=True)
with (outdir / 'stage6_hapiee_tag_holdout.csv').open('w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
(outdir / 'stage6_hapiee_tag_holdout_receipt.json').write_text(json.dumps({'vcf': url, 'selection': 'non-Finnish EUR', 'holdout': 'FIN', 'window_bp': 100000, 'limitation': 'Reference-panel check only; FIN is not a Czech cohort and array sites are position-only, without assay allele or strand verification. Not imputation quality or score validation.'}, indent=2) + '\n')
