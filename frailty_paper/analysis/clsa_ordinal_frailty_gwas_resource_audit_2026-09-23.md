# CLSA ordinal frailty GWAS resource audit — 2026-09-23

## Source and phenotype

Borhan et al. reported a 2026 ordinal GWAS of the Fried frailty phenotype using 23,105 European-ancestry participants from the Canadian Longitudinal Study on Aging (CLSA) comprehensive cohort. Participants were classified as non-frail, prefrail, or frail using the five Fried criteria. The study is a related frailty-phenotype resource, but it does not analyze a sleep exposure and therefore is not itself eligible for the registered sleep–frailty systematic-review question. Its operational definitions are study-specific: weight loss, CES-D exhaustion, PASE activity, 4-m walk time, and sex/BMI-stratified grip strength.

## What is reported

The paper reports one genome-wide significant lead SNP, rs147311617, with *P* = 4.98 × 10⁻⁸, and reports no validation sample. The authors state that they could not provide effect estimates and 95% confidence intervals because of computational cost, and did not estimate SNP heritability with their OrdinalGWAS workflow. These limitations prevent standard genome-wide rg/heritability QC and materially restrict effect-direction or quantitative replication checks. The paper's data-availability statement directs researchers to apply for access to de-identified CLSA data; the article does not identify a public full GWAS summary-statistics file or accession. The article links a supplementary PDF, but the ordinary publisher-content retrieval attempted for this audit did not return it; no file was archived or inspected.

## Interpretation for this project

CLSA is distinct from the UK Biobank/TwinGene FI discovery by source description, but exact participant overlap with any paired sleep dataset remains unknown. The phenotype is relevant as a Fried-style construct, but differences in operationalization, the small frail subgroup (1,547 participants), absence of effect estimates, no reported SNP heritability, and lack of a validation cohort mean it cannot currently support independent genome-wide sleep–frailty replication or a new LDSC analysis. Do not treat rs147311617 alone as replication of any project-wide genetic-correlation result.

There is a source-text inconsistency in the lead-signal cytoband: the title and Results use 12q22, while the abstract says 12p22. Until the underlying full summary data or an independently verified variant coordinate/build is available, do not propagate a cytoband-level locus claim from this article. The paper also has inconsistent PLXNC1/PLXCN1 spelling in its narrative; this audit does not endorse a gene assignment.

**Decision:** retain as contextual, phenotype-side evidence only; classify as `NOT_YET_ELIGIBLE` for pairwise replication. Reassess only if a complete, accessible summary-statistics source with effect estimates, standard errors, build, alleles, and sample/provenance metadata is verified, and if a compatible independent sleep dataset and overlap assessment are available. No CLSA application was submitted and no controlled data were accessed.

## Sources checked

- [Borhan et al., npj Aging (2026), version of record](https://www.nature.com/articles/s41514-026-00363-z)
- [Borhan et al., PubMed record (PMID 42020437)](https://pubmed.ncbi.nlm.nih.gov/42020437/)
- The article's linked supplementary-information PDF was attempted through the publisher's ordinary link; retrieval returned a cache-miss error in the web reader. This is an access limitation of this audit, not evidence that no supplement or data exist elsewhere.
