# Brain6 replication source access recheck

Audit date: 2026-09-25. Scope: whether official MVP summary-statistic access has changed enough to unblock the independent insomnia–MDD replication candidate.

## Evidence

- The current dbGaP public analysis list for `phs001672.v12.p1` identifies `pha005122.1` as `MDD.EUR.MVP.NatNeuro2021`, an ICD-code-defined European-ancestry MDD GWAS: [dbGaP analysis list](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetListOfAllObjects.cgi?object_type=analysis&study_id=phs001672.v12.p1).
- The current official [MVP Discover MVP Data page](https://www.mvp.va.gov/pwa/discover-mvp-data) says summary results for many projects are available through `phs001672.v11.p1`, while full summary-statistic downloads require a short dbGaP application. A public catalog entry therefore does not establish that the complete `pha005122.1` file is available without the access workflow.
- The separately listed `pha012639.1` analysis is maternal-grandparent depression in three ancestry groups (8,448 cases among 400,487 veterans), not participants' own ICD-coded MDD; it is not a phenotype-compatible replacement: [dbGaP analysis record](https://ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?pha=12639&study_id=phs002453.v1.p1).

## Decision

Keep the MVP European ICD-coded MDD replication candidate **BLOCKED** pending the official full-statistics access workflow and a source receipt that binds released bytes, ancestry, genome build, sample counts, and processing. No application was submitted and no access control was bypassed. No file was downloaded or used, and no replication estimate was computed. Existing cohort-overlap exclusions and the Brain6 global result remain unchanged.
