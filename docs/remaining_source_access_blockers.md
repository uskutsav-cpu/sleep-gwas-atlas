# Resolved source-access and ancestry blockers

No selected atlas-v1.0 source remains access-pending. The formerly blocked or
ancestry-ineligible records were replaced by public, European-ancestry sources
without accepting a download form, requesting controlled data, or changing the
locked 45-trait panel.

| Trait | Promoted source | Resolution and retained caveat |
| --- | --- | --- |
| `ms` | FinnGen R9 `G6_MS` | Public Finnish register endpoint, 2,182 cases and 373,987 controls. It is smaller and differently ascertained than the unavailable IMSGC discovery meta-analysis. |
| `asthma` | FinnGen R9 `J10_ASTHMA_EXMORE` | Public Finnish register endpoint, 42,163 cases and 202,399 controls. Its additional control exclusions differ from the replaced Han main-UKB phenotype. |
| `t2d` | FinnGen R9 `T2D` | Public Finnish register endpoint, 57,698 cases and 308,252 controls. It is smaller than the terms-gated DIAMANTE European meta-analysis. |
| `cad` | FinnGen R9 `I9_IHD` | Public Finnish wide ischaemic-heart-disease endpoint, 63,744 cases and 313,533 controls. It is broader than CAD and remains explicitly labelled as a proxy in all interpretations. |
| `melanoma` | FinnGen R9 `C3_MELANOMA_SKIN_EXALLC` | Public Finnish register endpoint, 2,993 cases and 287,137 cancer-free controls. It is smaller and uses a different control definition than the controlled Landi analysis. |
| `telomere_length` | Burren 2024 GWAS Catalog `GCST90435144` | Public non-Finnish-European UK Biobank stratum, N=438,351. It replaces the mixed-ancestry full-UKB Codd release, which is retained only as archived provenance. |

The six complete source files total 115,865,922 traversed data rows. Their
audits found zero malformed rows and confirmed autosomal coordinates, explicit
rsIDs, single-base SNP alleles, and valid effect, standard-error, P-value, and
frequency fields. FinnGen's published object MD5 values were reproduced; every
file also has a locally registered SHA-256 and passed gzip CRC validation. The
Burren file had no located upstream checksum, so its first acquisition used the
published HTTPS byte count before registering local MD5 and SHA-256 values.

All six sources are GRCh38 and use checksum-pinned, predeclared hg38-to-hg19
liftover plans for the EUR LDSC workflow. Source and schema readiness are now
45/45. Harmonization and LDSC h² are also complete for all six replacements;
the remaining work is downstream analysis, not source access.
