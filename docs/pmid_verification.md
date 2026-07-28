# PMID verification — Phase 0 trait manifest

Independent check of `phase0_trait_manifest - Phase0 manifest.csv`.

**Method.** Every PMID queried against the NCBI E-utilities `esummary`
endpoint (`db=pubmed`), comparing returned journal, year and title against
the source the manifest claims. Rows marked `RESOLVE` were searched with
`esearch` by DOI, author+title, or exact title. This is reproducible: the
queries are plain HTTP GETs against a public API, no key required.

**Date run:** 2026-07-28.

**Not verified here.** Sample sizes, case/control counts, genome build, and
ancestry were *not* checked — confirming that a PMID points at the right
paper says nothing about whether the numbers next to it are correct. Those
still require the papers themselves.

## Result: all 39 `VERIFIED` PMIDs are correct

Every PMID the manifest marks `VERIFIED` resolves to the claimed paper, in
the claimed journal, in the claimed year. No exceptions, no near-misses. The
manifest's verification discipline holds up.

Includes confirmation that **CDG3 = PMID 41372416**, *Nature* 2026, "Mapping
the genetic landscape across 14 psychiatric disorders" — matching the
citation in `CLAUDE.md` and `methods_map.md`.

## Newly resolved (11 of 15 `RESOLVE` rows)

| Row | Trait | PMID | Journal / year | Note |
|---|---|---|---|---|
| D2d | Autism | 30804558 | Nat Genet 2019 | Grove 2019, as expected |
| D2h | OCD | 40360802 | Nat Genet 2025 | Strom 2025, 30 loci |
| D2i | Substance use | 36477530 | Nature 2022 | GSCAN/Saunders, matches manifest DOI |
| D3g | Asthma | 32296059 | Nat Commun 2020 | Han 2020 |
| D3h | Atopic dermatitis | 37794016 | Nat Commun 2023 | Budu-Aggrey 2023 |
| D4d/D4e | Fasting glucose / insulin | 34059833 | Nat Genet 2021 | Chen 2021, "trans-ancestral genomic architecture of glycemic traits" |
| D6d | Frailty index | 34431594 | Aging Cell 2021 | Atkins 2021, matches manifest DOI |
| D6e | Epigenetic age | 34187551 | Genome Biol 2021 | McCartney 2021 |
| D7a | Breast cancer | 32424353 | Nat Genet 2020 | Zhang 2020 BCAC |
| D7d | Lung cancer | 28604730 | Nat Genet 2017 | McKay 2017 |
| D7e | Melanoma | 32341527 | Nat Genet 2020 | Landi 2020 |

## Still unresolved (4)

- **SL08 Snoring (Campos 2020).** The title given in the manifest returns
  nothing. The paper is likely *"Insights into the aetiology of snoring from
  observational and genetic investigations in the UK Biobank"* (Nat Commun
  2020) — the manifest's working title looks paraphrased. Confirm before
  citing.
- **D7b Prostate (Conti 2021).** Search reaches only PMID 33473200, a
  *Publisher Correction* to the trans-ancestry prostate meta-analysis. The
  primary article PMID still needs retrieving. Do not cite the correction.
- **D7c Colorectal (Fernandez-Rozadilla 2023).** Same situation: PMID
  36782065 is an *Author Correction* to "Deciphering colorectal cancer
  genetics through multi-omic analysis". Primary article PMID outstanding.
- **D2f Anxiety.** No candidate paper named in the manifest, so nothing to
  verify. Needs a decision on which release to use before it can be looked up.

## What this does not fix

`pop_prev` appears nowhere in the manifest, and `h2_SNP` is empty for every
row. Neither is a lookup — see the mentor questions in `README.md`. Phase 0
cannot close on the strength of this manifest alone.
