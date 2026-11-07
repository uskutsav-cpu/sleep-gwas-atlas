#!/usr/bin/env python3
"""Build a descriptive chromosome 5 case study from the prepared evidence."""
import csv, html
from pathlib import Path
P=Path(__file__).resolve().parent; O=P/'results'; CID='placo_insomnia__adhd_chr5_103447968_104447968'
with (O/'variant_evidence.tsv').open() as f: rows=[r for r in csv.DictReader(f,delimiter='\t') if r['candidate_locus_id']==CID]
with (O/'regions.tsv').open() as f: regions=list(csv.DictReader(f,delimiter='\t'))
r=next(x for x in regions if x['candidate_locus_id']==CID)
assert len(rows)==int(r['n_union_credible_variants'])
ranked=sorted(rows,key=lambda z:-max(float(z['insomnia_conditional_abf_pip'] or 0),float(z['adhd_conditional_abf_pip'] or 0)))[:8]
lines=['# Brain6 × AlphaGenome: chromosome 5 case study','', '**This report presents the original GWAS-based baseline.** See [frontal-cortex predictions](../jbrowse_predictions/FRONTAL_CORTEX_REPORT.md) and [AVI results](avi_mirror/AVI_PILOT_REPORT.md) for subsequent AlphaGenome analyses.', '', 'The frozen candidate is insomnia–ADHD, chr5:103,447,968–104,447,968 (GRCh37, one-based). Its 71 variant rows are the union of the separately reported exploratory 95% ABF credible sets; 71 passed the reference and coordinate checks. The eight ADHD credible-set members are included within that union. These are not independently confirmed causal variants.', '',f"The existing single-signal trait–trait coloc posterior PP.H4 is {float(r['PP_H4_p12_1e-5']):.4f} at p12=10⁻⁵ and {float(r['PP_H4_p12_1e-6']):.4f} at p12=10⁻⁶. These are model- and prior-dependent values, not evidence of a demonstrated mechanism. Brain6's LAVA family failed its QC gate; the alternative local analysis did not support this candidate at the frozen FWER threshold; no compatible independent two-trait locus replication was established.",'', '![Exploratory ABF PIPs](chr5_pip.svg)', '', 'The plot shows model-conditional PIPs only for the union of published credible sets. Points outside this selected union are not shown. Larger values indicate more posterior mass under the single-causal-variant ABF model; they are not measured functional effects.', '', '| SNP | GRCh38 variant | Insomnia PIP | ADHD PIP | Credible-set membership |','|---|---|---:|---:|---|']
for x in ranked:
    fmt=lambda v:f'{float(v):.4f}' if v else 'not available'
    membership=', '.join(t for t in ['insomnia','adhd'] if x[t+'_credible_set_member']=='YES')
    lines.append(f"| {x['snp']} | [{x['atlas_variant']}]({x['atlas_url']}) | {fmt(x['insomnia_conditional_abf_pip'])} | {fmt(x['adhd_conditional_abf_pip'])} | {membership} |")
lines += ['', 'The displayed rows are ordered by the larger of the two conditional PIPs solely to make this case study readable. It is not a shared-causality score. Alleles in the GWAS beta columns refer to the source A1; the mapped GRCh38 REF/ALT is a separate orientation.', '', 'This baseline table does not use AlphaGenome output or assert a gene, tissue or mechanism. The later pilot scored all selected variants with overall AVI and obtained tissue predictions for one variant per region. All remaining five regions appear in `regions.tsv` and `variant_evidence.tsv`.','']
(O/'CHR5_CASE_STUDY.md').write_text('\n'.join(lines))
print('Wrote case study')
