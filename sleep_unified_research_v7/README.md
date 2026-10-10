# Sleep GWAS Atlas V7: empirical resource benchmark

V7 is a separate, bounded evaluation of source-aware resource coverage and evidence interpretation. It leaves V4/V5/V6 estimates, source records, and testing families unchanged. V7 runs no GWAS or LDSC calculations and cannot establish that a historical association is novel.

The protocol was frozen before the first V7 cross-resource match or task answer. See [`protocol/empirical_benchmark_preregistration_v1.md`](protocol/empirical_benchmark_preregistration_v1.md) and the dated, result-aware access decisions in [`protocol/protocol_deviations.md`](protocol/protocol_deviations.md).

## Build source-aware cross-resource matches

The builder reads frozen V4 results, the already-extracted Morrison/Goodman tables, the official Human GWAS ATLAS release 3 files, and the three previously verified Fan 2026 publisher supplements from an attached SSD. The Fan workbooks are not copied into this repository. The output stores atlas estimates, source metadata, comparator row references, and conservative match classes; it never copies the Fan rg/SE/P/FDR values. SleepChart is carried forward from the completed V4 audit as partial/related comparator context, without claiming a V7 row-level join or re-downloading the workbook.

Requires Python 3.10+ and `openpyxl` for reading Fan workbooks. From the repository root, supply the verified external input directories/files:

```sh
python3 sleep_unified_research_v7/scripts/build_cross_resource_matches.py \
  --gwas-atlas-dir "/path/to/verified/human_gwas_atlas_release3" \
  --fan-data2 "/path/to/verified/fan_sd2.xlsx" \
  --fan-data4 "/path/to/verified/fan_sd4.xlsx" \
  --fan-data11 "/path/to/verified/fan_sd11.xlsx"
```

The builder checks the frozen 1,637-record universe, comparator table dimensions, Human GWAS ATLAS pair-key integrity, and Fan SD11 pair-key uniqueness. It records every candidate row reference. Independent reviewers found that PMID, approximate sample-size, ancestry, and label agreement did not support 30 initial Human GWAS ATLAS exact-match labels; these are classified E pending source-release, ascertainment/coding, effect-direction, estimator/reference, and estimand evidence. Goodman SD24 contains six composite scores. The build manifest records source, code, protocol, and generated-table hashes.

Publisher workbooks carry their own rights restrictions. Fan's article states CC BY-NC-ND 4.0. Although V7 does not export Fan statistics, derived row references and classifications still require human rights review before public redistribution. Source tables remain with their original publisher/download terms.

Run `python3 sleep_unified_research_v7/scripts/render_resource_coverage.py` after the build to regenerate the vector figure. See `reports/resource_benchmark_interpretation.md` and `reports/authoring_handoff.md` for final coverage counts, the scientific interpretation, the seven-priority status, review findings, and human release gates. V7 establishes no exact HGA source/estimand matches and claims no new association, method, independent replication, or human usability benefit.
