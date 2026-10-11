# Predefined source-discrepancy adjudication

Date: 11 October 2026. Independent automated scientific review; not human peer review. This report records a material source-direction interpretation defect and the limits of its resolution. It does not report a new association, independent replication, effect-size contrast, mechanism, or causal result.

## Scope and selection

The parent assignment fixed five comparisons before extracting HGA values: insomnia–IBD, insomnia–Crohn's disease, long sleep–Crohn's disease, long sleep–ulcerative colitis, and a concordant insomnia–depression control. They were chosen to investigate a previously observed discrepancy, not to estimate its prevalence across resources. All 25 loci in Table 1 of de Lange 2017 were then used as source-side positive controls, without selecting on their sleep associations. One indel and two strand-ambiguous SNVs were excluded before source-direction comparison, leaving 22 eligible controls. This is a documented case study, not an externally registered benchmark.

The review read actual frozen V4 core full-precision results, the locally verified HGA v20191115 pair/source tables, preserved de Lange raw/harmonized/munged files, and the primary de Lange paper. It performed no LDSC/GWAS fit, no GWAS source download, and no change to historical files. Third-party numerical extracts and the source article remain private on the SSD.

## Source identity is stronger than PMID agreement, but equivalence is unproven

HGA IDs 2029/2030/2031 specify the same de Lange 20161107 Crohn/UC/IBD source basenames and disease-specific sample counts used by the atlas. ID 3786 explicitly specifies the Jansen UK Biobank-only insomnia source, excluding 23andMe, with the atlas's case/control counts. ID 4273 specifies the Dashti long-sleep source and the same tail-stratum counts. ID 4293 specifies Howard's UKB+PGC depression release excluding 23andMe, with the same case/control counts.

These metadata establish candidate release and stratum identity more directly than matching PMID alone. They do not prove that HGA consumed the same file bytes, interpreted the same effect allele/phenotype direction, used identical filtering or an identical estimator/reference contract. No HGA processed summary statistics or source-specific munging log was available; HGA documentation says processed bodies are not distributed. Exact V7 class A source/estimand equivalence remains uncertified. [HGA documentation](https://atlas.ctglab.nl/documentation)

The named sources also clarify phenotype meaning: Jansen UKB is frequent initiation/maintenance complaints, not diagnosed chronic insomnia; Dashti long sleep is >=9 reported hours versus 7–8-hour controls, not clinical hypersomnia. The sleep sources share UKB participants, and Crohn/UC/combined IBD analyses share source cohorts. These estimates cannot supply independent two-trait replication. The depression source combines broad UKB help-seeking with primarily clinical PGC ascertainment, so it should not be described as an entirely clinically homogeneous MDD phenotype. [Dashti primary paper](https://www.nature.com/articles/s41467-019-08917-4)

## Concrete direction defect

The four bowel comparisons have opposite point-estimate signs in HGA and V4, generally with similar magnitudes. The insomnia–depression control is concordant. This pattern alone was insufficient to assign fault or infer a biological difference.

Primary-paper validation changes that assessment. De Lange Table 1 identifies a risk allele, odds ratio above one, phenotype, GRCh37 coordinate, and association P value for each selected locus. For every one of the 22 eligible nonpalindromic SNV controls, interpreting raw `Effect` as a disease-increasing effect for `Allele1` gives the opposite direction to the published risk allele. This holds in 13/13 IBD, 5/5 Crohn, and 4/4 UC controls. Raw P values agree with the table to its rounding (maximum relative discrepancy 0.283%). Exponentiating the negative of the raw effect oriented to the published risk allele recovers the published risk OR within the table's two-decimal rounding (maximum absolute difference 0.00418). Thus this is not explained by an arbitrary sign in an uncertain cross-trait point estimate. [De Lange primary paper, Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC5289481/); [institutional accepted manuscript](https://www.repository.cam.ac.uk/bitstream/1810/262817/1/deLange_et_al-2017-Nature_Genetics-AM.pdf)

Four of those unambiguous paper controls were retained in the frozen HapMap3 products (three IBD, one Crohn). All four retain the direction conflict in both harmonized beta and munged Z after accounting for allele order. None of the UC paper anchors was retained, so no retained UC paper-anchor claim is made. A separate, value-independent sample of the first 20 frozen UC harmonized records shows exact preservation of raw Allele1/Allele2/Effect and compatible munged-Z signs in 20/20 records. That check establishes pipeline preservation on sampled records; it does not independently validate UC disease direction.

The current schema assertion that raw `Effect` per `Allele1` points toward increased disease liability is therefore contradicted by external primary-paper controls. Copying those same fields into a GWAS Catalog formatted file does not resolve the contradiction. The historical runner accurately reproduces its input convention; reproduction accuracy and correct phenotype direction are separate properties.

## Mechanism remains partly unresolved

The de Lange Methods report SNPTEST v2.5 additive association followed by standard-error weighted METAL (2011-03-05) meta-analysis. Official SNPTEST documentation defines additive beta for allele_B, whereas METAL relies on explicitly supplied effect-allele labels and effect fields. This makes an allele-label convention mismatch plausible. The original SNPTEST-to-METAL conversion script, exact METAL configuration, and definitive source README connecting final `Allele1/Allele2` to tested allele were not recovered here. The controls support reversed disease orientation under the atlas's current assumption, but do not uniquely distinguish an effect-per-Allele2 convention from reversed phenotype coding or an upstream label translation. Do not assert a proved source-author mistake or an HGA software defect. [SNPTEST v2.5 documentation](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest_v2.5.html); [METAL documentation](https://genome.sph.umich.edu/wiki/METAL_Documentation)

A mere allele-order change is not a correction unless its beta/Z is re-expressed consistently. Reversing both traits' signed statistics leaves rg unchanged; reversing only one phenotype/allele effect convention changes rg sign. Ordinary positive liability-scale or sample-size scaling cannot explain this global direction reversal.

## Narrow scientific consequence and remedy

Direction claims for all 36 core rg rows using these three de Lange outcomes (12 sleep constructs × 3 bowel traits) should be quarantined pending a separately versioned source-semantic remedy. Their frozen historical values should remain available as reproduced historical results. No rows in the current 62-sensitivity or 41 outcome-side validation tables use these three outcomes by source trait ID. This does not certify every downstream artifact: figures, literature classifications, evidence tables, candidate release, and biological-use-case text that depend on these 36 directions need targeted reconciliation.

The remedy should document the published-risk-anchor evidence, the final effect-allele/disease-direction contract, and whether it swaps alleles or negates signed effects. It should verify the corrected contract against the positive controls and unchanged N, variant identity and filtering, then validate one bounded corrected native fit before regenerating dependent signed presentations. Do not simply edit the public rg column to match HGA. Under an exact uniform re-expression of one trait, rg/Z/genetic covariance/cross-intercept and their intervals change sign, while SE/two-sided P/BH significance and single-trait h2 remain unchanged; that invariance must be checked on actual consumed inputs before it supports a remedy.

The long-sleep–UC V4 interval includes zero, so its opposite point sign is not evidence of a resolved statistical effect or a between-resource difference. No between-resource difference P value was calculated for any case: source compatibility and covariance of the estimates are unestablished. Any renewed direction statement should use corrected direction with its original uncertainty and testing family.

The demonstrated utility is specific: following source definitions and primary risk-allele controls prevents a misleading disease-direction interpretation that numeric reproduction alone would not detect. It is not a representative usability benchmark or proof of superiority over HGA. The technical correction, once independently confirmed, improves the scientific validity of the atlas; it is not itself new biological discovery.

## Outputs

- `public_safe_case_adjudications.tsv`: five-case decisions and atlas-owned uncertainty metadata; excludes HGA numerical estimates and source-paper numerical rows.
- `private_predefined_case_values.tsv`: actual paired numeric estimates and source references, retained privately on SSD.
- `private_hga_selected_source_metadata.tsv`: six exact HGA source metadata rows, private.
- `private_delange_paper_direction_controls.tsv`: all 25 Table 1 source controls, eligibility and raw-to-munged trace, private.
- `private_uc_pipeline_preservation20.tsv`: bounded UC preservation trace, private.
- `direction_trace_summary.json` and `uc_pipeline_preservation_summary.json`: derived count summaries.
- `trace_predefined_cases.py` and `trace_uc_pipeline_preservation.py`: read-only extraction scripts.

The original V4 analyses, HGA snapshot, historical source bodies, and prior reports were not changed.
