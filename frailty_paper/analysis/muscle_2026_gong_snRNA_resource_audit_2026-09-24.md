# 2026 human skeletal-muscle snRNA-seq resource audit

## Decision

Register this study as a candidate human muscle-aging resource, but do not treat it as an acquired processed single-cell dataset or use it for downstream functional interpretation at this stage. The project currently has no eligible frailty-specific loci.

## Verified study facts

The publisher article (Gong et al., *Aging Cell*, first published 2026-04-14, DOI [10.1111/acel.70485](https://doi.org/10.1111/acel.70485)) reports eight male Chinese participants: five adults aged 22–60 and three adults aged 99–101. It reports SPPB, DXA skeletal muscle index, hand-grip strength and walking pace; these functional measurements do not constitute a validated frailty instrument. The methods report 80,277 captured nuclei/cells, with 71,373 retained after quality control. The assay is 10x Chromium 3-prime snRNA-seq; alignment is to GRCh38, Ensembl 105. The paper's data statement cites NGDC GSA-Human HRA009010 and says raw data were made public.

## Access and reuse limits

The article does not identify a processed expression-matrix location. The supplement is described as a 36.4-MB Word document of figures. General [GSA-Human documentation](https://ngdc.cncb.ac.cn/education/tutorials/gsa-human/) describes that archive as a secure human-data service with DAC-reviewed access; the accession-level record and access state for HRA009010 were not verified in this audit. Therefore the article's statement is retained as author-reported public availability, not independently confirmed open access. The article's CC BY license does not establish separate reuse terms for archived sequence data. No sequence files were downloaded.

The eight-donor, male-only, extreme-age design is descriptive and has limited generalizability. It is neither an independent genetic replication nor a frailty GWAS, and it contains no sleep phenotype or QTL results. Any cell-context use would remain exploratory and would require eligible upstream loci plus accession/file metadata.

## Repository changes

Registered as `GONG_2026_CENTENARIAN_MUSCLE_SNRNA` in `manifests/muscle_single_cell_resource_index.tsv` and `MUSCLE_AGING_2026_GONG` in `manifests/targeted_resource_audit.tsv`. Access and processed-matrix status are explicitly unverified; no data were acquired.
