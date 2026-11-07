"""Descriptive systematic sample; association tests need not be Wald tests."""
import csv, gzip, hashlib, json, math, sys
from pathlib import Path

name = sys.argv[1]
base = Path(__file__).resolve().parent
receipt = json.loads((base / name / 'receipt.json').read_text())
source = Path('/Volumes/Extreme SSD/brain6-public-component-gwas-20261029') / ('GWAS_' + name + '_maf0.01_info0.9.txt.gz')
h = hashlib.sha256()
with source.open('rb') as f:
    while chunk := f.read(1048576):
        h.update(chunk)
if h.hexdigest() != receipt['source_sha256']:
    raise ValueError('source changed since checksum audit')

counts = {'rows': 0, 'sampled': 0, 'valid': 0, 'p_underflow': 0,
          'wald_chisq_rel_diff_gt_1e-3': 0, 'wald_chisq_rel_diff_gt_0_1': 0,
          'logp_abs_diff_gt_0_01': 0, 'logp_abs_diff_gt_0_1': 0}
max_rel = max_logp = 0.0
with gzip.open(source, 'rt') as f:
    for i, row in enumerate(csv.DictReader(f, delimiter='\t')):
        counts['rows'] += 1
        if i % 100:
            continue
        counts['sampled'] += 1
        beta, se, chi, logp = (float(row[x]) for x in ('BETA', 'SE', 'CHISQ', 'LOG10P'))
        if not all(map(math.isfinite, (beta, se, chi, logp))) or se <= 0 or chi < 0:
            continue
        counts['valid'] += 1
        wald = (beta / se) ** 2
        rel = abs(wald - chi) / max(chi, 1e-8)
        max_rel = max(max_rel, rel)
        counts['wald_chisq_rel_diff_gt_1e-3'] += rel > 0.001
        counts['wald_chisq_rel_diff_gt_0_1'] += rel > 0.1
        p = math.erfc(math.sqrt(chi / 2))
        if p == 0:
            counts['p_underflow'] += 1
            continue
        difference = abs(-math.log10(p) - logp)
        max_logp = max(max_logp, difference)
        counts['logp_abs_diff_gt_0_01'] += difference > 0.01
        counts['logp_abs_diff_gt_0_1'] += difference > 0.1
out = {'source_sha256': h.hexdigest(), 'stride': 100, 'counts': counts,
       'max_wald_chisq_relative_difference': max_rel, 'max_logp_difference': max_logp,
       'interpretation': 'Systematic sample, not independent tests; a non-Wald association test may legitimately differ from beta/SE.'}
(base / name / 'statistic_check.json').write_text(json.dumps(out, indent=2) + '\n')
print(name, json.dumps(out), flush=True)
