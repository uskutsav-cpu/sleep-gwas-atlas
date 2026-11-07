# Reference LD for selected sleep and depression molecular leads

Using 503 European participants in the 1000 Genomes Phase 3 reference, the two FinnGen depression–RNA top shared variants showed **moderate**, not near-perfect, linkage with the Brain6 AlphaGenome candidate rs62264764:

| Published depression–RNA top variant | r² with rs62264764 | r² with insomnia lead rs1456193 |
|---|---:|---:|
| rs1585084, hippocampus feature | 0.397 | 0.291 |
| rs9822630, anterior-cingulate feature | 0.349 | 0.425 |

The two depression markers have r²=0.880 with each other in this reference. The two reported molecular lead variants, rs13073833 and rs4597702, have identical observed haplotypes to each other and to rs9822630 in these 503 people. The two insomnia-region markers rs1456193 and rs62264764 have r²=0.821, consistent with the earlier project check. Frequencies, counted alleles, signed correlations and all 15 pairwise comparisons are in [ld.tsv](ld.tsv).

These values describe linkage in a small external European reference, not the FinnGen, UKB, ADHD or GTEx cohorts. A top shared variant is not necessarily the causal variant; different causal variants, multiple signals and study-specific LD can produce similar or different regional patterns. Moderate LD neither confirms nor rules out a shared signal, gene mediation or a three-trait sleep–ADHD–depression mechanism. The selected depression features were chosen after inspecting prior findings, and the two molecular–depression comparisons reuse GTEx data.

[Code](calculate.py) records unique forward-strand Ensembl mappings from GRCh38 to GRCh37, exact VCF allele matches and phased calls. All six sites had 503 eligible reference participants. [Receipt](receipt.json) saves the raw mapping responses, rsIDs and source details; [verification](verification.json) checks the saved arithmetic. Remote VCF/index access was regional and not whole-file hashed. Rounding beyond ±1 is clipped after a numerical-tolerance check. The prior [FinnGen supplement review](../RESULTS.md) records why CLPP and insomnia coloc PP4 cannot be equated.
