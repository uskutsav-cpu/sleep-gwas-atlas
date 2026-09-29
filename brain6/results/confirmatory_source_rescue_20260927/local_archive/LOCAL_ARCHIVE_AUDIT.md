# Exact Dashti long-sleep local archive audit, 2026-09-27

## Question and frozen contract

Can a local historical copy supply the missing native per-variant analyzed sample size for the exact Dashti et al. 2019 UK Biobank contrast, self-reported sleep **≥9 h versus 7–8 h** (34,184 cases and 305,742 controls), without replacing the locked phenotype? The original source and every archived raw candidate examined here have the ten-column header `SNP CHR BP ALLELE1 ALLELE0 A1FREQ INFO BETA_LONGSLEEP SE_LONGSLEEP P_LONGSLEEP`. None contains native `N`, `NMISS`, `N_analyzed`, `OBS_CT`, case/control counts, or per-SNP missingness. `A1FREQ` and `INFO` are present but do not determine analyzed N.

## Audited candidate files

The [machine-readable inventory](LOCAL_LONG_SLEEP_CANDIDATES.tsv) hashes seven separate loose source containers on the Extreme SSD and two historical tar members. The three ZIP archives are byte identical: SHA-256 `0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885`. Four loose gzip copies have differing compressed hashes from repackaging. Streaming SHA-256 of the **uncompressed** content in all four gives the same `09e8dea3fb192f3021894e470a3130370ef152c9e5174f18e50abb2aafc36db8`, 1,094,455,329 bytes and 14,661,601 data rows, identical to the ZIP member. The July tar and 56-GB August tar.gz members hash as the same gzip container as the August loose and September Utsav copies. In the July tar, `longsleep_dashti.txt.gz` is a **hard link** to `longsleep.txt.gz`, not a second export.

Directories examined included the current Brain6 repository and `/Users/swethasunilkumar/Downloads`, the Extreme SSD `brain6-work/atlas-raw`, `sleep-gwas-atlas-frailty-v1`, `Codex Archive/2026-08-26-sleep-gwas-atlas`, `Utsav-Research-Archive/Sleep-GWAS`, `Mac Offload 2026-09-04/Research Downloads` and `Codex Archives`, and `Codex-Archive/2026-08-26.tar.gz`. Seven dated/offload tar files were listed without extraction. The `Research Downloads/sleep-gwas-atlas.tar` contains a derived LDSC `longsleep.sumstats.gz` but no original long-sleep raw source file. The 2026-07-28 tar contains the raw gzip and hard-link alias; the August tar.gz contains the same raw gzip bytes and derived files. The `Downloads/sleep-gwas-atlas-main` copies expose paper/supplement and analysis logs but no raw long-sleep file.

## Derived-file trap

The archived September `data/harmonized/longsleep.harmonized.tsv.gz` and LDSC `data/munged/longsleep.sumstats.gz` *do* have `N` columns. They are downstream transforms, not new source information. The harmonized first row has `N=122985.40891841165`; its archived QC receipt says `sample_size_schema_mapping=CONSTANT_N_EFF_FROM_MANIFEST`, with all 6,549,769 retained rows in that mode. The LDSC first row rounds this to `122985.409`. Thus these files cannot supply a source-native N or validate the BOLT-LMM export's sample-size semantics.

## Authorized cohort access

The project [external blocker record](../../lava_longsleep_source_rescue_v1/LONG_SLEEP_EXTERNAL_BLOCKER.md) states a prior documentation search found no Brain6-approved UK Biobank application or technically available authorized individual-level environment. The local [replication feasibility audit](../../../qc/replication_source_feasibility_20260923.md) records that MVP-only results are tied to `phs001672` and the dbGaP/MVP application route; no application was submitted. This audit found no additional approval receipt in the reviewed Brain6 manifests, QC notes, or source documents. This is a statement about project evidence, not a claim that the investigators have no institutional access. No credentials or restricted participant data were inspected.

## Local admission decision

**No richer exact source was found locally.** The only raw exact-phenotype content identified is the byte-identical original release in multiple containers. The archived `N` values were calculated by a later harmonizer from study-level case/control counts. The local search therefore does not resolve the missing per-variant analyzed N or exact released-model transform. It admits no new long-sleep LAVA input and leaves the protected and canonical outputs untouched. The next admissible evidence remains a custodian export/clarification or an authorized exact-phenotype rerun with verified per-variant N and model provenance.

## Read-only verification method

Files were inventoried by path, ZIP central directory, tar headers, gzip/ZIP first line, SHA-256 of compressed containers, and streaming SHA-256 of decompressed content. No archive was extracted to the internal disk and no source bytes were modified. The prior [original release inventory](../../lava_longsleep_source_rescue_v1/original_release_inventory_v1.json) supplies the independent full-field audit of the same uncompressed hash.
