# Derived-workspace capacity recheck — 2026-09-26 12:24 UTC

The configured analysis workspace is `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace`. `make -C frailty_paper storage-check FRAILTY_ANALYSIS_WORKSPACE='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace'` passed, reporting 1,629 GiB free against the 20 GiB minimum. The volume reports 1.6 TiB available. The internal repository filesystem has 7.5 GiB free; keep large temporary and derived outputs on the configured external workspace.

This is an operational capacity check only; it does not resolve data, review, source-access, replication or scientific QC gates.
