from pathlib import Path
import csv,hashlib
out=Path('outputs')
def read(name):
 with (out/name).open(newline='') as f:return list(csv.DictReader(f))
def p(v):
 a=float(v)
 return f'{a:.3g}' if a>=.001 else f'{a:.2e}'.replace('e-0','e−').replace('e-','e−')
r=read('prototype_variant_five_trait_evidence.csv');tri={x['variant']:x for x in read('prototype_five_source_descriptive_triage.csv')};d=read('direct_fi_ffs_genetic_correlation.csv')
assert len(r)==len(tri)==15 and len(d)==2
rows=[]
for x in r:
 name=x['variant']; q=tri[name]; values=[p(x[k+'_p']) for k in ['FI','FFS','Jansen_insomnia','R13_insomnia','MVP_insomnia']];rows.append('| '+ ' | '.join([name,x['counted_allele'],*values,'yes' if q['all_five_nominal_p_lt_0_05']=='True' else 'no'])+' |')
text='''# Insomnia and frailty: risk comparison and candidate evidence

Research briefing for review • 27 October 2026

## Key findings

1. **FI and FFS have substantial genome-wide genetic correlation in the project's existing data.** Auxiliary LDSC estimates were 0.757 (SE 0.021; 95% CI 0.715–0.799) and 0.755 (SE 0.020; 95% CI 0.716–0.794) on two matched SNP sets. They use the same two frailty GWAS, so they are sensitivity checks, not independent replications. These are genetic correlations, not correlations between measured scores in the same people or percentages of shared individual risk.
2. **Insomnia correlates more strongly with this FI source than with this FFS source** in each of three source-specific comparisons. Why remains unknown. The UK Biobank component of the FI contains an insomnia item and a fatigue item; FFS includes exhaustion. Score content, participant overlap, age, ascertainment and other differences are plausible explanations that this analysis cannot separate.
3. **There is no validated 15-marker shared insomnia–frailty risk score.** A research-only prototype and 15-marker evidence table have been prepared. The markers were chosen from a published post-selected cross-trait list and their weights have not been tested in independent genotyped participants. They do not yield an individual's calibrated risk.

## FI versus FFS

| Feature | Frailty Index (FI) | Fried Frailty Score (FFS) |
|---|---|---|
| Construct | Accumulation of multiple health deficits | Five physical domains: weakness, slow walking, exhaustion, low activity and weight loss |
| Source used here | Atkins et al. FI meta-GWAS: 164,610 UK Biobank + 10,616 TwinGene participants, European descent | Ye et al. quantitative 0–5 FFS GWAS: 386,565 UK Biobank participants, European descent |
| Content concern | UKB FI includes an insomnia and a tiredness item; TwinGene's FI is not identical | Exhaustion overlaps conceptually with fatigue; detailed original field/cutoff crosswalk remains unresolved |
| Study overlap | FI and FFS both include UK Biobank; exact participant intersection unknown | Same limitation |

Sources: [Atkins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/), [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/). The genetic-correlation evidence is reproduced from the project's saved native estimates; the FI–FFS point estimate and 200-block jackknife SE were recalculated from those outputs, not refit from raw GWAS. See `direct_fi_ffs_genetic_correlation.csv` and `fi_ffs_risk_followup.md` for checksums and input details.

## Insomnia correlations with each measure

| Insomnia source | Insomnia–FI r_g | Insomnia–FFS r_g | FFS minus FI | 95% CI for difference |
|---|---:|---:|---:|---:|
| UK Biobank symptom measure | 0.640 | 0.483 | −0.156 | −0.213 to −0.100 |
| FinnGen clinical insomnia | 0.426 | 0.245 | −0.181 | −0.260 to −0.102 |
| MVP European clinical insomnia | 0.549 | 0.398 | −0.151 | −0.220 to −0.082 |

The difference uses matched SNP support and paired jackknife covariance. All three pass the existing conservative 396-test Bonferroni overlay. The same FI and FFS GWAS are reused, and the comparison has not tested cross-source differences. The original UKB symptom contrast had paired-uncertainty validation; the two clinical contrasts additionally passed the saved three-trait matrix checks. These are descriptive comparisons of specified data, not independent replication of frailty, evidence of causality, or proof that removing an FI item would eliminate the difference.

## Fifteen-marker risk evidence

The 15 rows come from Song et al.'s published cross-trait list, taking the 30 lowest reported CPASSOC p-values and retaining those available in the three project inputs after existing QC. This is a **post-result, availability-driven research set**. It is not LD-pruned or an independently validated score. Nine of the 15 rsIDs have more than two alleles listed in the Ensembl reference. The FI and FFS inputs overlap, exact overlap with insomnia sources remains unknown, and the raw non-HapMap3 FinnGen/MVP markers have not undergone full project QC. Source p-values below are **unadjusted**. A p-value describes a statistical association, not a person's chance of disease or a causal effect.

| Variant | Counted allele | FI p | FFS p | UKB insomnia p | FinnGen insomnia p | MVP insomnia p | All five p<0.05? |
|---|---|---:|---:|---:|---:|---:|---|
'''+ '\n'.join(rows)+'''

Two variants—rs12536335 and rs7734045—show p<0.05 in all five sources with the same aligned sign. This is a **descriptive observation**, not a multiple-testing-adjusted discovery: the markers were selected after earlier results, the data are not fully independent, and LD and allele representations need validation. None of the 15 satisfies p<0.05/15 in all five sources; that numerical benchmark is also not a valid post-selection corrected claim. Non-significance is not evidence of no effect.

A cautionary example is rs113851554: despite strong insomnia p-values in all three sleep sources, its FI p is 0.1987 and its FFS p is 0.19 in the project inputs. Its sleep evidence must not be presented as a confirmed shared frailty association.

The 124 original Song candidates and the full 15-marker estimates, SEs and approximate intervals are preserved separately. Song's cross-trait approach does not prove that each listed variant independently affects both traits or that a shared causal variant exists. [Song et al.](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full). The project's own FI×sleep LAVA run completed all 29,940 slots, but all 12 prespecified failure gates failed; those results cannot be used to establish or exclude local sharing.

## Score status and what the published scores show

The `unvalidated_15_variant_score_prototype.csv` gives allele-aligned FI, FFS and insomnia coefficients for three separate raw weighted sums. It has not been applied to a person; it has no validated scaling, calibration, discrimination, missing-genotype handling, LD pruning, or independent evaluation. The three phenotype-specific weights are on different scales and must not be added into a single score. The prototype is not for clinical or individual risk use.

Published genome-wide scores exist: [PGS005228](https://www.pgscatalog.org/score/PGS005228/) was derived from FI and [PGS005229](https://www.pgscatalog.org/score/PGS005229/) from FFS. In the published study these were evaluated against FI in older-adult cohorts; reported performance varies by cohort. Their whole-score results cannot be assigned to the 15-marker subset and do not establish a joint insomnia–frailty score. [Foote et al.](https://www.nature.com/articles/s41588-025-02269-0)

A file audit found that all 2,379,429 rows in the original FI PGS file repeat the effect allele in the `other_allele` field. This does not prove incorrect weights or invalidate the published study, but the two columns cannot be treated as a distinct allele pair. In a deterministic sample, 2,322 of 2,323 score rows matched a project FI variant at the same position, and every matched counted allele agreed with one allele in that FI input. This supports the counted allele within that filtered sample only; it is not a validation of all source rows or predictive performance. The original score was not changed.

## What would resolve the main uncertainty

* **Why FI and FFS differ:** obtain appropriately authorized, same-participant FI and FFS data; reconstruct source-authentic or clearly named alternative scores; compare full FI with FI excluding the insomnia item and then insomnia plus fatigue, using covariance-aware estimates and a fixed analysis plan. Exact source scoring and adequate power must be verified first.
* **A genuine shared score:** define the outcome and time horizon, select and estimate weights with an LD-aware, prespecified procedure, then evaluate calibration and incremental performance against existing FI, FFS and insomnia scores in independent genotyped cohorts. No such evaluation is available here.
* **Candidate loci:** retain the 15 and 124 tables as hypotheses. A separately justified, QC-passing analysis with dense inputs, ancestry-matched LD, overlap review and fine-mapping/colocalization is required before claiming a shared causal variant. The failed LAVA gate remains in force.

## Supporting files

* `prototype_variant_five_trait_evidence.csv` — 15 variants, five source effects and intervals where SEs exist.
* `prototype_five_source_descriptive_triage.csv` — transparent numerical screen; not adjusted significance.
* `published_insomnia_fi_124_cross_trait_candidates.csv` — complete transcribed published candidate list.
* `unvalidated_15_variant_score_prototype.csv` — exploratory coefficients and allele warnings, not a clinical score.
* `published_fi_ffs_pgs_source_audit.json` and `fi_pgs_allele_sample_summary.json` — published-score provenance and allele audit.
* `fi_ffs_risk_followup.md` — full methods, source limitations, example calculations and links.
'''
(out/'insomnia_frailty_risk_briefing.md').write_text(text)
print('Wrote consolidated briefing:',len(r),'variant rows')
