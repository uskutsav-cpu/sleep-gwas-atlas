# Phase 0/1 results for the six public substitutions

The six replacement sources were streamed against the exact HapMap3 allowlist,
lifted from GRCh38 to GRCh37 with the registered UCSC chain, filtered under the
fixed Phase 0 rules, and munged with the pinned LDSC runtime. Complete raw
sources and full-source audit ledgers remain retained; the streaming prefilter
only bounds the LDSC working set.

| Trait | Source rows | HapMap3 prefilter | Harmonized | Allele-matched effects | h² scale | h² (SE) | Z | Intercept | Gate |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| `ms` | 20,170,207 | 1,196,088 | 1,159,339 | 1,159,114 | liability | 5.9823 (1.1958) | 5.00 | 1.0288 | PASS under the predefined Z/intercept gate; liability magnitude >1 is a major interpretation warning |
| `asthma` | 20,165,155 | 1,196,085 | 1,159,258 | 1,159,033 | liability | 0.2342 (0.0161) | 14.55 | 1.1568 | PASS |
| `t2d` | 20,170,006 | 1,196,088 | 1,159,319 | 1,159,094 | liability | 0.3329 (0.0176) | 18.91 | 1.3153 | DROP: intercept above 1.20 |
| `cad` | 20,170,236 | 1,196,088 | 1,159,334 | 1,159,109 | liability | 0.1346 (0.0100) | 13.46 | 1.1618 | PASS; endpoint is the broader IHD proxy |
| `telomere_length` | 15,022,702 | 1,204,235 | 1,174,303 | 1,174,105 | observed | 0.0634 (0.0038) | 16.68 | 1.0893 | PASS |
| `melanoma` | 20,167,616 | 1,196,087 | 1,159,216 | 1,158,991 | liability | 3.7310 (1.0456) | 3.57 | 1.0343 | DROP: h² Z below 4; liability magnitude >1 is also nonphysical |

The locked sleep×non-sleep family contains all 396 pairs. The 372 pairs in
which both traits pass h² QC are labelled `PRIMARY_PHASE1`; the 24 T2D or
melanoma pairs are labelled `QC_FAILED_SENSITIVITY` and
`EXCLUDED_FROM_PRIMARY_INFERENCE`. Benjamini–Hochberg correction over the full
396 tests yields 153 primary FDR<0.05 pairs; the primary-only 372-test family
yields 155. The complete-family table preserves both corrections and the exact
trait-level QC reason for every sensitivity row.

The very high MS and melanoma liability-scale estimates arise in extremely
imbalanced, low-case-count Finnish endpoints and must not be interpreted as
plausible variance proportions. No prevalence, sample-size convention, or gate
was changed after observing the estimates. Their raw LDSC logs, observed-scale
quantities used internally by cross-trait LDSC, and sensitivity labels remain
available for the robustness audit.
