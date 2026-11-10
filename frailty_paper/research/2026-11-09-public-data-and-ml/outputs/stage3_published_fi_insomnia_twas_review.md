# Published FI and insomnia TWAS: shared gene-level leads

I compared the **four** genes appearing in both cross-tissue UTMOST and whole-blood FUSION in [Lin et al.'s FI study](https://doi.org/10.3389/fgene.2024.1404456), Table 1, with the **15** genes in the corresponding UTMOST/FUSION overlap table of [Li et al.'s insomnia study](https://doi.org/10.1152/jn.00490.2024). Two gene symbols occur in both published tables: **LRPPRC and CCDC134**. This is an overlap of two author-selected lists, not a new TWAS, genome-wide cross-phenotype test, or proof of a shared causal variant. The reported association values are transcribed in [the comparison](stage3_published_fi_insomnia_twas_intersection.csv).

| Gene | FI: UTMOST p / whole-blood FUSION p | Insomnia: UTMOST p / whole-blood FUSION p | FUSION Z: FI / insomnia |
|---|---|---|---|
| LRPPRC | 3.32 × 10⁻⁵ / 1.96 × 10⁻⁷ | 1.46 × 10⁻⁴ / 2.14 × 10⁻⁴ | −5.20 / −3.70 |
| CCDC134 | 1.94 × 10⁻⁴ / 7.25 × 10⁻⁵ | 4.62 × 10⁻⁵ / 0.00206 | 3.97 / 3.08 |

The FI paper reports its whole-blood and cross-tissue analyses at FDR < 0.05; the insomnia paper describes FDR < 0.10. For CCDC134, the insomnia FUSION FDR is 0.07 and therefore does not meet an FDR < 0.05 threshold. Similar FUSION Z signs are descriptive gene-model associations and do not establish the same causal effect, gene mediation, direction of intervention, or a joint effect across traits. Both studies rely on the same or overlapping UK Biobank-related GWAS and GTEx prediction resources used elsewhere in this project; their methods are not independent cohorts or four independent validations.

Lin et al. report FI–whole-blood eQTL coloc H4 = 0.98 and FOCUS PIP = 0.997 for LRPPRC under their own model. This is **FI–molecular evidence**, not insomnia–molecular colocalization or a three-trait analysis, and was not reanalyzed here. Their published descriptions of the colocalization threshold and MR interpretation warrant careful source-method review; neither a TWAS association nor a reported MR association alone proves the gene causes a trait. LRPPRC is an additional gene-level lead **outside the 22-marker shortlist**; do not silently add it as a candidate-variant hit.

## A direct link to the existing MMRN1-labeled marker

The insomnia TWAS paper also reports **MMRN1** in its 15-gene table: cross-tissue UTMOST p = 1.38 × 10⁻⁶ (reported FDR 4.29 × 10⁻⁴) and whole-blood FUSION Z = −4.03, p = 5.48 × 10⁻⁵ (reported FDR 0.01). The project’s rs7691487 is labeled near MMRN1 and has an exact-marker GTEx v10 spleen eQTL for MMRN1, providing a reason to examine this region. The TWAS is a gene-model signal, not a result for rs7691487. A gene-name match and a nearby variant do not establish they represent one causal signal, or that MMRN1 is shared with FI/FFS. [Existing variant–gene table](stage3_initial_variant_gene_table.csv).

The insomnia paper’s FOCUS reporting is internally inconsistent: its text attributes PIP 0.998 to MMRN1, while Table 2 assigns **0.849 to MMRN1** and **0.998 to VRK2**. Record both and do not choose one as verified without the authors' underlying results. FOCUS gene-model PIP is distinct from a variant PIP or probability of illness. The paper also describes these analyses as using UK Biobank insomnia `ukb-b-3957`, so its association with the project's UKB-based findings is not independent.

## Scope and provenance

The intersection was computed from the published **Table 1** lists and recorded with the published rounded statistics. It does not cover all 36 UTMOST FI genes, 107 FUSION FI genes, 195 insomnia UTMOST genes, or 332 insomnia FUSION genes reported elsewhere in the papers; absence from these short tables is not evidence of no TWAS association. The source tables differ in test families, thresholds, study definitions and tissue models. The comparison does not discover new variants, infer independent loci or validate a risk score. Script: `work/stage3_published_twas_intersection.py`; output SHA-256 `64e8a20a2b23286c32151dc8e08670312fdc9eb9c2ff42727ebab7bea7951f98`.

A full analysis would use complete downloadable TWAS models and results, harmonized GWAS and LD references, and multi-signal GWAS–QTL colocalization in relevant tissues, with explicit phenotype and sample-overlap accounting. These published results are an additional evidence layer, not a substitute for that work.
