# Brain6 cohort-overlap follow-up: long sleep–schizophrenia

**Audit date:** 2026-09-25. **Scope:** source-level cohort documentation for the exact locked pair; no participant-level data or identity checks were performed.

## Locked sources

- Long sleep: Dashti et al. 2019, UK Biobank self-reported sleep duration tail (PMID 30846698; [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6405943/)).
- Schizophrenia: Trubetskoy et al. 2022 PGC3 Wave 3 European autosomal stratum, 53,386 cases and 77,258 controls (PMID 35396580; [Nature article](https://www.nature.com/articles/s41586-022-04434-5)). The locked source card identifies the European autosomal file `PGC3 SCZ wave3. european. autosome.public. v3.tsv.gz`.

## Evidence review

The primary article’s data-availability statement identifies the public PGC source for its core, extended, ancestry-specific, and sex-stratified summary statistics. Its linked [Supplementary Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-022-04434-5/MediaObjects/41586_2022_4434_MOESM1_ESM.pdf) contains a “Case-control sample descriptions” section (PDF pages 35–54) describing the contributing schizophrenia samples. That source description does not name UK Biobank as a contributing cohort. The article’s separate references to UK Biobank concern comparator GWAS, not the listed PGC3 schizophrenia sample.

This establishes that the checked publication does not document UK Biobank cohort membership for the selected schizophrenia GWAS. It does not prove zero participant overlap: the public summary-statistic file does not provide individual-level identifiers or an exact participant intersection, and the cohort descriptions cannot exclude undocumented sample reuse or control overlap.

## Decision

Keep `longsleep__scz` classified `UNKNOWN` in `manifests/cohort_overlap.tsv`. Record the narrower finding as **no UK Biobank cohort named in the checked PGC3 source descriptions; participant-level overlap remains unverified**. Do not call this pair independent on this evidence alone. No GWAS statistics were downloaded or analyzed, and no global rg, replication, LAVA, PLACO, or downstream result was changed.
