# Confirmatory LAVA rescue planning, before additional source selection

**RESCUE_PLANNING_ONLY. Frozen canonical v3 remains FAILED_QC_NOT_PROMOTED.**

The immutable seven-trait × 2,495-locus aggregate has 13,745 TESTED, 3,720
NOT_RUN, and zero FAILED trait-locus cells. The frozen limit is 873 NOT_RUN;
a valid rescue must reduce the count by at least 2,847 without altering the
family, thresholds, or criteria. No bivariate v3 tests were executed.

## Trait failure matrix

| Trait | Planned | TESTED | NOT_RUN | FAILED | Numeric testable | Low local h² | Shared-reference K | Component K |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| longsleep | 2,495 | 1204 | 1291 | 0 | 48.26% | 1271 | 20 | 0 |
| insomnia | 2,495 | 1824 | 671 | 0 | 73.11% | 651 | 20 | 0 |
| parkinson | 2,495 | 1859 | 636 | 0 | 74.51% | 613 | 23 | 0 |
| mdd | 2,495 | 1906 | 589 | 0 | 76.39% | 565 | 24 | 0 |
| adhd | 2,495 | 2272 | 223 | 0 | 91.06% | 200 | 23 | 0 |
| bipolar | 2,495 | 2327 | 168 | 0 | 93.27% | 145 | 22 | 1 |
| scz | 2,495 | 2353 | 142 | 0 | 94.31% | 119 | 22 | 1 |

Observed reasons total 3,564 low-local-h², 154 fewer shared-reference variants,
and two insufficient components. These LAVA outcomes do not identify whether
underlying GWAS power, variant coverage, harmonization, allele ambiguity,
LD-reference compatibility, or binary N handling caused a given low-h² cell.
The sparse-reference group is a direct overlap/representation problem; the
source-side gap at loci 950–969 is already documented separately. No canonical
numerical failure was observed. Sample overlap cannot explain this trait-only
univariate missingness, but any replacement changes the pairwise overlap
contract and requires new review.

## Pair eligibility and unexecuted tests

| Frozen pair | Planned | Bivariate TESTED | Bivariate NOT_RUN by family gate | Bivariate FAILED observed | Both numeric | At least one NOT_RUN | Both strict gates passed |
|---|---:|---:|---:|---:|---:|---:|---:|
| insomnia__adhd | 2,495 | 0 | 2,495 | 0 | 1687 | 808 | 0 |
| insomnia__mdd | 2,495 | 0 | 2,495 | 0 | 1423 | 1072 | 0 |
| longsleep__scz | 2,495 | 0 | 2,495 | 0 | 1153 | 1342 | 1 |
| longsleep__bipolar | 2,495 | 0 | 2,495 | 0 | 1135 | 1360 | 0 |
| longsleep__parkinson | 2,495 | 0 | 2,495 | 0 | 913 | 1582 | 1 |

The 12,475 pair slots were not executed because the family failed first.
Zero observed bivariate failures therefore does not establish bivariate
numerical stability. The detailed pair-reason combination TSV accounts for
every locus and both observed trait reasons without double-counting pair
eligibility. Numeric testability and Bonferroni strict-gate eligibility are
reported separately.

## Source-independent rescue bound

Ranking by NOT_RUN contribution alone gives longsleep, insomnia, parkinson,
mdd, adhd, bipolar, then scz. Perfectly repairing the three largest leaves
1,122 NOT_RUN, still above 873. The unique smallest perfect-repair set has
four traits: insomnia, longsleep, mdd, and parkinson; it would leave 533.
This is an arithmetic lower bound, not evidence that any such replacement
exists or will pass. MDD alone leaves 3,131. Engineering the 156 sparse-K
cells alone leaves 3,564. Do not run a full rescue family until source
semantics, multi-trait pilots, and a pre-run receipt justify it.
