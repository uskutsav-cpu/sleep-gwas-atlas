# Published walking-pace comparison: insomnia and ADHD

[Timmins et al. 2020](https://www.nature.com/articles/s42003-020-01357-7) analyzed self-reported walking pace in 450,967 European-ancestry UK Biobank participants, coding slow, steady/average and brisk as 0, 1 and 2. This individual physical-function trait is relevant to the slowness component of physical frailty and does not itself include an insomnia or exhaustion question. It is **not** a reconstructed FFS without exhaustion or an FI with sleep items removed.

The following are the authors' published LDSC estimates, extracted from their Supplementary Data 8, not new calculations by Brain6. Negative values mean the genetic direction associated with faster self-reported pace is correlated with the genetic direction associated with the other trait; they do not establish causality.

| Other trait | Walking pace r_g (SE), p | Walking pace adjusted for BMI r_g (SE), p |
| --- | --- | --- |
| Insomnia | −0.2998 (0.0227), 1.129×10⁻³⁹ | −0.2264 (0.0247), 4.4213×10⁻²⁰ |
| ADHD | −0.3170 (0.0283), 3.336×10⁻²⁹ | −0.1307 (0.0313), 2.9363×10⁻⁵ |
| Hand grip strength | 0.2606 (0.0197), 4.1156×10⁻⁴⁰ | 0.2962 (0.0228), 1.3493×10⁻³⁸ |

All six displayed p-values are below the paper's stated Bonferroni threshold of 0.05/53 ≈ 9.43×10⁻⁴. The paper reports 53 traits; this threshold does not establish a separately corrected family across every comparison or every analysis in this atlas. BMI-adjusted and unadjusted estimates use related data: the apparent differences are descriptive, not tests of a difference or proof of BMI mediation. The article reported continued insomnia correlation after BMI adjustment. These prior findings mean a broad claim that genetic overlap between insomnia and a physical-frailty-related component is new would be incorrect.

The table identifies insomnia as the Jansen et al. 2019 meta-analysis (N=1,331,010), ADHD as Demontis et al. 2018 (N=53,293), and hand grip as Tikkanen et al. 2018 (N=334,825). Those are the source table's sample sizes and references, not independently reconstructed participant counts. The insomnia meta-analysis includes UK Biobank, as does the walking-pace study; the exact overlap and handling of cross-trait dependence were not independently established here. The ADHD and grip references are distinct input studies, not proof of independent sleep–ADHD–frailty locus replication. Slow walking can be affected by many health and behavioral factors, and walking pace is a subjective measure, so these results do not isolate purely physical biology or a sleep-independent frailty construct.

The [published supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-020-01357-7/MediaObjects/42003_2020_1357_MOESM4_ESM.xlsx) was saved on the SSD at `/Volumes/Brain6Work/walking_pace_published_supplement.xlsx`. Its MD5 matched the value in the Europe PMC article XML; its SHA-256 and extracted row positions are in the [receipt](receipt.json). [Extraction code](extract.py) reads SD8 rows 30, 33 and 49 without editing the workbook; [extracted data](published_correlations.tsv) retain the source labels and numbers. This single published comparison does not resolve the [matched component-exclusion analysis](../CONSTRUCTION_PROTOCOL.md), local-correlation failures, or the need for independent validation.
