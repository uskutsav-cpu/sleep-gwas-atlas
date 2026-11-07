---
license: other
license_name: alphagenome-services-additional-terms
license_link: https://deepmind.google.com/science/alphagenome/terms
tags:
  - genomics
  - variant-effect-prediction
  - alphagenome
size_categories:
  - n>1T
---

# AlphaGenome AVI scores, re-encoded

AlphaGenome's **Variant Impact (AVI)** scores for **8,812,917,339 SNVs** on GRCh38, re-encoded from
the 88.5 GB published tabix TSV into ~34 GB of parquet by
[`just-dna-enricher`](https://github.com/just-dna-seq/just-dna-format).

**This is a re-encoding, not a re-analysis.** No score is changed, recomputed or filtered.

## What is in it

| | |
| --- | --- |
| `data/alphagenome_avi-<contig>.parquet` | `chrom`, `pos` (1-based VCF), `ref`, `alt`, `raw_score_e5` |
| `avi_knots.parquet` | the `PHRED` reconstruction curve — **not optional**, see below |
| `release.json` | row counts, the source digest, and the artifact's own timestamp |
| `LICENSE.txt` | the Output Terms' "Use restrictions" section, verbatim |

**`raw_score` is stored as `Int32` at a scale of 10^5.** The published file prints at most five
decimals, so the integer is exactly lossless where a `Float32` would be both larger and lossy.
**Compare in the integer domain** — `score >= 0.1` is `raw_score_e5 >= 10_000`. Dividing back
(`raw_score_e5 / 1e5`) disagrees with the printed value on ~53% of rows, because the division rounds
a second time.

**`PHRED` is not stored.** Measured over every row it is an exact within-corpus rank — `PHRED >= p`
keeps `10^(-p/10)` of the corpus — so it is a function of `raw_score`, and storing it costs 24.7 GB.
`avi_knots.parquet` carries the curve in 466 KB instead, as an **interval** per printed score rather
than a point, because the source prints `raw_score` to four significant digits and `PHRED` to six.

That interval makes threshold safety **decidable in advance**: a threshold is unsafe iff it falls
inside some knot's span. Genome-wide, exactly **one** knot straddles any integer threshold from 1 to
50 — `0.00076`, 676,356 rows, spanning `PHRED` 2.99961 to 3.00027. Every other integer threshold is
decided.

**Absence is row-absence.** AVI covers ~95% of the assembly and includes 672,931 genuine zeros, so a
position with no row is *unscored* and a row holding zero is *scored zero*. They are not the same.

## Terms — read these

Governed by the **AlphaGenome Services Additional Terms of Service** and the **AlphaGenome Output
Terms of Use**. The AVI SNV scores are classified by the AlphaGenome download page as a *Permissive
Use Downloadable Artifact*, **for commercial and non-commercial use**; the merged-splicing and
feature-importance artifacts are **not** and are absent here deliberately.

`LICENSE.txt` in this repository carries the Output Terms' "Use restrictions" section as an
enforceable provision, which restriction 3b requires of anyone attaching their own terms. **By using
this data you agree to the AlphaGenome Output Terms of Use** at
<http://deepmind.google.com/science/alphagenome/output-terms>.

Two bars worth stating plainly because they are easy to miss:

- **No training** of machine-learning models for predicting genetic variant effects (prohibition 4).
- Google may **request deletion** of Output and Derivatives in your possession on breach — not only
  on termination.

The applicable version of the terms is the one **effective on the date the Output was generated**,
which is why `release.json` records the source artifact's own timestamp.

## Citation

> Avsec et al., Advancing regulatory variant effect prediction with AlphaGenome, Nature 649(8099):1206-1218, 2026. doi:10.1038/s41586-025-10014-0

## Provenance

Built by `just-dna-enricher alphagenome build --input <the artifact>` from the publisher's own
download; the source digest is in `release.json`. Nothing in this repository was fetched by a tool —
the source is behind a sign-in whose eligibility clause bars classes of holder, so acquisition is the
operator's act under their own acceptance of the terms.
