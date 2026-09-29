# Public pQTL resource access recheck

Audit date: 2026-09-24. Scope: current official/publisher resource pages for the two indexed 2023 proteomics pQTL resources. No statistics or protein-assay files were downloaded.

## UK Biobank Pharma Proteomics Project

The AWS Registry of Open Data identifies the UKB-PPP, 54,219 UK Biobank participants and 2,923 proteins, and gives a CC BY license label. On the same page, however, the described population-specific GWAS summary-statistic S3 resource is marked “Controlled Access.” The Synapse project page is JavaScript-rendered to this reader. A project discussion documents an entity-view route for querying file metadata (`syn53038826`) and states that users can search/list the files and download selected files by Synapse ID. It does not establish current permissions or reuse terms for any specific pQTL file. The registry license label therefore cannot be used alone to conclude that a given object is anonymously downloadable or carries that license.

Operational disposition: after a qualifying locus exists, locate only its relevant protein assay and ancestry-specific file, then confirm the object's access state, license, schema, coordinate build, effect scale and sample size before retrieval or analysis. Do not download the bundle or treat this UK Biobank sample as independent replication.

## deCODE 2023 plasma proteomics

The publisher paper’s data-availability statement says GWAS summary statistics for all 2,931 Olink and 4,907 SomaScan assays are available from deCODE. The official deCODE summary-data listing identifies the combined 2023 plasma-proteomics release as 19.7 TB. The linked download form is JavaScript-only in the current text view; the source pages do not establish a locus-targeted transfer method, object-level license, per-file genome build, or file schema.

Operational disposition: retain it as a potential non-UKB molecular-context resource, but check overlap with Icelandic aging/frailty outcomes and confirm exact file metadata before use. Do not initiate a 19.7 TB transfer. No UK Biobank overlap is not equivalent to independence from all Icelandic cohorts.

## Sources

- AWS Registry of Open Data, UKB-PPP: <https://registry.opendata.aws/ukbppp/>
- Synapse UKB-PPP project: <https://www.synapse.org/Synapse:syn51364943>
- Synapse discussion on querying specific files: <https://www.synapse.org/Synapse:syn51364943/discussion/threadId=10184>
- UKB-PPP primary article: <https://doi.org/10.1038/s41586-023-06592-6>
- deCODE/SomaScan/Olink comparison article and data availability: <https://doi.org/10.1038/s41586-023-06563-x>
- deCODE official summary-data listing: <https://www.decode.com/summarydata/>
- deCODE 2023 release form: <https://download.decode.com/form/largescaleplasma-2023>
