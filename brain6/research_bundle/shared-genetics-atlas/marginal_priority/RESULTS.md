# Marginal-association patterns in common HapMap3 markers

This post-hoc screen identifies markers to review, not independently replicated or independent genetic loci. Across 1,138,262 shared markers and three traits, the Bonferroni cutoff for the full 3,414,786 SNP–trait family is 1.4642×10⁻⁸. The 368 unique SNPs with at least one crossing have these patterns:

| Pattern of marginal threshold crossings | SNPs |
|---|---:|
| FinnGen R13 insomnia only | 30 |
| FI only | 21 |
| FFS only | 313 |
| FI and FFS, not insomnia | 4 |
| Insomnia and either frailty measure | 0 |

The table shows **significance patterns**, not trait-specific effects: a non-crossing can reflect limited power or measurement differences, and the FI includes sleep items. SNP counts cannot be interpreted as locus counts because neighboring markers may be in linkage disequilibrium. No sleep–frailty association can be ruled out by this table; the separately described max-p analysis tests their conjunction using its own two-pair family.

Examples for follow-up include insomnia-ranked rs9611542 (GRCh37 chr22:41677317; insomnia p upper bound 4.12×10⁻¹⁰, FI 0.0223, FFS 0.00210); FI-ranked rs362273 (chr4:3227419; FI p bound 3.35×10⁻¹⁰, FFS 3.02×10⁻⁸, insomnia 0.556); and FFS-ranked rs13107325 (chr4:103188709; FFS p bound 7.62×10⁻¹⁹, FI 0.00797, insomnia 0.669). These examples were selected after the scan, are not LD-clumped or independently replicated, and their signs and p-values depend on the listed effect alleles and source methods. The four FI+FFS SNPs map to a short interval on chr4 and may tag a common signal; their number does not establish four loci. The complete list with alleles, all three Z scores and p-value bounds, and positions is `marginal_threshold_crossings.tsv`.

The normal-tail p bounds conservatively offset the three-decimal input Z by 0.0005; they do not fix model miscalibration or source rounding in beta and SE. Positions are joined to the FI source by unique rsID and matching unordered allele pair, not independently validated with a reference genome; all 368 matched. Input hashes, source paths, row counts, code and the pre-result plan are in `summary.json`, `screen.py` and `PROTOCOL.md`. This is a common HapMap3 subset, not dense whole-genome coverage, and excludes HFRS. Source population ancestry, sample overlap, case ascertainment and power differ. Fine mapping with study-appropriate LD and independent cohorts is needed before a locus is called shared or specific.
