# Human skeletal muscle ageing atlas: processed-object access and terms audit

Audit date: 2026-09-24

## Scope

This audit checked the official Human Skeletal Muscle Aging Atlas portal, the linked publisher article, the study code repository, and the raw-data accession. It did not download or inspect the H5AD object bytes.

## Findings

- The author-hosted portal labels one full human-muscle H5AD and six subset H5AD links as “Processed (h5ad)”. A direct content request for the full human object and the MuSC subset returned HTTP content type `application/x-hdf` to the research browser. This verifies that the linked endpoints serve HDF content, but does not verify successful complete transfer, archive integrity, schema, object contents, or reuse terms.
- The portal summary describes the atlas as 17 donors (8 young, approximately 20–40 years; 9 old, approximately 60–80 years). The publisher article reports 90,902 single cells and 92,259 single nuclei from 17 donors and explains that most donors received CPAP or mechanical ventilation during their hospital stay. The project index preserves the article/supplement age details and records the hospital-context limitation.
- The article’s open-access statement specifies CC BY 4.0 for the article and material included in that article, subject to third-party credit-line exclusions. The data availability statement says that processed objects are hosted at the atlas portal. Neither the article nor the portal page reviewed here explicitly states that the separately hosted H5AD objects themselves are CC BY 4.0. Therefore the article license cannot be assumed to license those downloadable objects.
- The linked GitHub repository says it contains code to reproduce analyses and points to the portal for processed AnnData objects. Its Apache-2.0 license applies to repository contents; it does not establish an object-specific license for the H5AD downloads.
- ArrayExpress E-MTAB-13874 is identified as the repository for raw sequencing data from newly generated libraries. This is separate from the processed H5AD routes and does not resolve their terms.

## Disposition

Keep the atlas indexed as a useful public discovery route for later locus-to-cell-context work, but mark processed-object license as unverified and do not redistribute or build an analysis package around downloaded H5AD files until the portal terms or rights-holder confirms applicable reuse terms. Do not download raw FASTQs for this project stage. No genetic association, frailty phenotype, or independent replication is provided by this atlas.

## Primary sources

- Atlas portal and its H5AD links: <https://www.muscleageingcellatlas.org/>
- Kedlian et al., *Nature Aging* 4, 727–744 (2024), DOI: <https://doi.org/10.1038/s43587-024-00613-3>
- Author code repository: <https://github.com/Teichlab/SKM_ageing_atlas>
- Raw library accession: <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13874>
