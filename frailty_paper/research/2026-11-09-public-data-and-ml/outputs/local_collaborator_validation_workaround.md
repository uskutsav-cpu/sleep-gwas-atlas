# Local collaborator route for sleep–frailty score evaluation

Research draft, 1 November 2026. This is a route around transferring individual-level data, not around consent, access rules or study approval. No collaborator has agreed and no validation has been performed.

## The practical proposal

Ask an authorized cohort investigator to run the project’s frozen scoring and analysis plan inside their approved environment. They would retain genotypes and phenotypes and share only disclosure-reviewed results. A collaboration or amendment may still be required, and compute or staff costs are unknown. Encrypted federated computation is possible but unnecessary for a single cooperating site; local analysis is simpler.

**Prioritize LBC1936** for a scientific feasibility conversation: its published work reports a frailty index, and separate research describes physical frailty and sleep assessments. The candidate Wave 3 items, published FI item audit and access caveats are documented in [the LBC1936 draft](lbc1936_validation_research_request_draft.md). Its modest size and potentially small number of frailty events limit precision; the linked sample and score-variant coverage have not been established. It may already appear in frailty follow-up literature, so specifically investigate participant and relative overlap with *every* discovery and selection dataset and all previous uses of this project's candidates before calling the test independent.

**Framingham is an alternative** with sleep and physical-function data in public dictionaries. Ask a researcher who already has an appropriate project whether they can obtain approval for this exact analysis; see [exam-specific proposal](framingham_mentor_access_proposal.md). The FHS says external dbGaP and BioLINCC downloads have no data-access service fees, but institutional and computing costs and approvals remain possible ([FHS](https://www.framinghamheartstudy.org/fhs-for-researchers/data-available-overview/)). Availability of linked eligible participants and exact variants is not demonstrated. Check overlap with previous CHARGE gait work and all other discovery and selection sources.

**NHATS is another alternative** if a collaborator already has approval for the genetic enclave. NHATS says there is currently no data-access charge but genetics, restricted geography and CMS restricted data cannot be combined ([NHATS](https://nhatsannex.org/restricted-data)). Project-specific permissions and the genotype–survey linkage need confirmation; see [feasibility](stage6_nhats_genetic_validation_feasibility.md). Follow all restrictions including any prohibition on using AI with NHATS data.

## Request to a potential collaborator — unsent draft

> Dear [researcher],
>
> We are studying genetic associations of sleep symptoms and frailty. We have exploratory 15- and 22-variant research score definitions and are looking for a genuinely independent evaluation. Would your team be interested in assessing feasibility and, if permitted, running a prespecified analysis locally in your approved environment? We would not request individual-level genetic or clinical data.
>
> Before any outcome analyses, could you advise whether the cohort has jointly usable genotypes, sleep measures and either a published frailty index or physical-frailty components; whether the data can be used for this question; whether there is overlap with the relevant discovery and candidate-selection studies; and whether staff time or other costs would apply? We can provide a variant manifest, provenance, score definitions and a draft protocol for review. We would welcome your input on phenotype definitions and appropriate authorship and approvals.
>
> Thank you, [name]

## Minimum analysis and release gates

1. Agree the primary score and outcome, cohort/visit, date, exclusions and models before accessing test outcome associations. Document prior explorations; a new registration does not erase earlier target lookups. Separate 15-marker FI-, FFS- and insomnia-weighted scores from the 22-marker equal-count and weighted variants; the 22-marker conventions differ. Neither is an established clinical test.
2. Resolve original effect-allele conventions, especially the FFS weight direction, and each target's genome build, strand and exact biallelic or multiallelic contrast. Do not use rsID alone. If an essential allele remains unresolved, do not compute the affected score. Determine in advance what to do for absent/low-quality variants; a partial score must be identified as such.
3. Establish approved linked genotype–phenotype counts, missingness, events, ancestry, family structure, usable marker coverage, and the time ordering of any longitudinal endpoint. Verify overlap with all training, selection and previous target sources. Establish feasibility before testing; no assumption that all available individuals have complete data.
4. Prespecify one frailty endpoint and report FI and Fried-type measures separately. Where possible request item-level derivations of the published FI plus sleep-excluded and sleep-and-fatigue-excluded sensitivities, on the same participants, with denominators and missingness rules. A modified FI requires its own description and does not by itself provide a sleep-excluded GWAS. A physical-frailty phenotype usually includes exhaustion, so it is not automatically fatigue-free.
5. Evaluate association and uncertainty beyond appropriate age, sex, ancestry, genotyping and kinship covariates. If the sample and events permit, evaluate discrimination, calibration and incremental performance against the same baseline and participants, with uncertainty. Avoid overfitting, retuning weights or asserting a negative result from an underpowered test. Report all prespecified scores and outcomes and relevant deviations.
6. Release only cohort-approved aggregate results, including effect scales, confidence intervals, sample flow and marker QC, consistent with disclosure rules. No individual records, small-cell disclosures, private files or external AI processing unless the cohort explicitly permits them.

## What the public evidence already accomplishes

Public GWAS summary-based analyses have already been performed in the [score evidence inventory](stage6_score_validation_evidence_summary.md). Their LD-adjusted weighted-Z checks are postselected approximations and do not measure score calibration, individual prediction, independent replication, or a joint sleep–frailty outcome. NHANES provides observational sensitivity evidence only. The public resources do not remove the need for an independent phenotype-linked test to establish prediction.

The prepared [15-variant score](unvalidated_15_variant_score_prototype.csv), [22-variant equal-count score](unvalidated_22_shared_direction_score.csv), [22-variant weighted score](unvalidated_22_variant_score_prototype.csv), and [allele input gate](stage6_score_allele_input_gate.md) are working research materials, not confirmed analysis-ready clinical instruments.
