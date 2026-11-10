from pathlib import Path
import re
import shutil

source=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-29/investigate-whether-we-can-directly-reconstruct-2/work/sleep-gwas-atlas/frailty_paper/paper/manuscript.md')
target=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/outputs/nhanes_manuscript_audit/manuscript_with_nhanes.md')
text=source.read_text()
def before(marker, addition):
    global text
    if text.count(marker)!=1: raise ValueError(f'marker not unique: {marker}')
    text=text.replace(marker,addition+'\n\n'+marker,1)

old='''84-pair family. A separate FI × 12-sleep local-sharing'''
new='''84-pair family. A separate retrospective NHANES observational sensitivity analysis examined sleep-trouble disclosure and an adapted FI with direct sleep and fatigue items omitted. A separate FI × 12-sleep local-sharing'''
assert text.count(old)==1
text=text.replace(old,new,1)
old='''No independent
pairwise sleep–frailty replication, eligible physical-component comparison,
complete family-level local-sharing result, shared-locus result, or molecular
follow-up is available.'''
new='''No independent
pairwise sleep–frailty genetic replication or participant-level score validation is available. Exploratory follow-ups do not establish a complete family-level local-sharing result or a validated shared causal locus.'''
assert text.count(old)==1
text=text.replace(old,new,1)
old='''37 also met the
same-family Bonferroni sensitivity. The refreshed PubMed retrieval contains'''
new='''37 also met the
same-family Bonferroni sensitivity. In NHANES, sleep-trouble disclosure was associated with the sleep/fatigue-excluded adapted FI (adjusted difference 0.0698, 95% CI 0.0528–0.0869; n=3,336). This observational result does not evaluate genetic overlap. The refreshed PubMed retrieval contains'''
assert text.count(old)==1
text=text.replace(old,new,1)
old='''dominance, unavailable alternative frailty GWAS and incomplete literature
screening limit interpretation.'''
new='''dominance, unavailable matched sleep/fatigue-excluded GWAS statistics for this analysis and incomplete literature
screening limit interpretation.'''
assert text.count(old)==1
text=text.replace(old,new,1)

before('### GWAS sources, phenotype hierarchy and cohort overlap', '''### Exploratory NHANES phenotype-overlap sensitivity

Separately from all genetic analyses, we used public NHANES 2011–2014 data to examine whether the observational association of reported sleep trouble with a health-deficit index persists after removing directly overlapping sleep and fatigue items. The analysis was retrospective and not prospectively registered. The exposure was whether a participant had ever told a doctor or other health professional about trouble sleeping (SLQ050); it is not a clinical insomnia diagnosis. The outcome was an adapted 49-item NHANES frailty index based on Liu et al. (2025; https://doi.org/10.1186/s12877-025-05752-9), alongside a 47-item version omitting the two depression-questionnaire items on sleep or excessive sleep and tiredness/little energy. Each version was the mean of observed deficit scores and required at least 80% of its own items; the 47-item version required at least 38. We retained participants meeting all four FI versions' requirements for the paired comparisons. This adapted FI is not the Atkins FI or the original Fried physical phenotype.

We created the full survey design using positive examination weights (WTMEC2YR divided by two), cycle-specific masked strata and primary sampling units, then restricted to adults aged 60 years or older. Survey-weighted linear models adjusted for age, sex and cycle; the stated confidence intervals and p-values use 32 design degrees of freedom. A secondary complete-covariate analysis additionally adjusted for reported race/ethnicity, education, income-to-poverty ratio and smoking; this is not a causal adjustment set. Direct component omission cannot remove indirect overlap through health conditions, medicines, disability or care-seeking. Full item definitions, missingness, source provenance and sensitivity analyses are supplied in the NHANES supplementary audit.''')

before('### Frozen Frailty Index estimates and workflow verification', '''### Exploratory observational phenotype-overlap sensitivity

Among 3,472 older adults with positive examination weights, 3,336 had the sleep disclosure item and met the completeness rule for all compared index versions; 969 reported ever disclosing sleep trouble. The survey-weighted age-, sex- and cycle-adjusted difference in the adapted 47-item FI was 0.0698 (95% CI 0.0528–0.0869; p=1.54×10⁻⁹) for participants reporting sleep-trouble disclosure versus those not reporting it. The analogous difference for the full 49-item index was 0.0736 (95% CI 0.0566–0.0907). Index composition and denominators differed, so the difference between these coefficients is not a causal decomposition. In the common complete-covariate subsample (n=3,035), the more extensively adjusted difference was 0.0677 (95% CI 0.0512–0.0843). The independently rerun computation reproduced these figures; using the survey-regression residual-df convention instead of the reported design-df convention gave 0.0698 (95% CI 0.0527–0.0869; p=4.42×10⁻⁹).

These results reduce concern that direct sleep and tiredness items fully account for this **observational** relationship under this adapted NHANES measure. They do not quantify or correct phenotype overlap in the Atkins FI GWAS, establish genetic correlation or causality, or validate the genetic scores. Other health deficits can still reflect sleep-related conditions, fatigue, medication use and healthcare contact.''')

before('### Conclusion', '''### Observational sensitivity and genetic inference

The NHANES sensitivity study is a useful separate test of item overlap in observed data, with an independently reproducible survey analysis and a transparent 47-item adaptation. The result is cross-sectional, retrospectively defined and vulnerable to health-selection and healthcare-disclosure effects. Its 80%-observed rule permits variable denominators. Mak et al. reported a GWAS of a modified FI excluding eight mental-wellbeing items including insomnia and fatigue (https://doi.org/10.1111/acel.13868); the full summary statistics for that modified phenotype have not been identified in their listed release. That broader exclusion cannot isolate the contribution of the two sleep and fatigue items. Obtaining those statistics and establishing sample comparability could enable a broader genetic sensitivity analysis; a same-participant GWAS excluding only sleep and fatigue would address the more specific question. NHANES does not provide either GWAS or independent score validation in the present analysis.''')

before('### Conclusion', '''### Exploratory regional functional follow-up

Separately from the frozen FI analyses, an exploratory postselected FOXP2-region comparison used MVP insomnia and FinnGen hospital frailty summary statistics with external ancestry-oriented LD. In a one-megabase window, one modelled signal pair had PP.H4=0.979 under p1=p2=10⁻⁴ and p12=10⁻⁵; PP.H4 changed to 0.322 at p12=10⁻⁷. The narrower window omits a distal frailty signal, and wider analyses were sensitive to the LD reference. MVP sample size was approximated, and ancestry-matched in-sample LD and complete cohort independence were not established. This model-dependent statistical follow-up uses a distinct hospital-record phenotype and does not establish a shared causal variant, an effect on FOXP2, a relevant cell type, or score validation. It includes no molecular QTL in this signal-pair comparison. Detailed diagnostics and limitations are in the [regional audit](functional_followup/foxp2_gnomad_coloc_review.md).''')

before('### Conclusion', '''### Published observational context in an independent setting

In a previously published study of 4,744 community-dwelling adults in Taiwan's Healthy Aging Longitudinal Study, Wu et al. used repeated-measures models and reported an adjusted association between insomnia and Fried-criteria frailty (odds ratio 1.41, 95% CI 1.16–1.72; https://doi.org/10.1016/j.jagp.2023.01.028). This is the authors' published observational result, not an analysis of our genetic scores. It concerns a physical frailty phenotype rather than the adapted NHANES or Atkins FI; its sleep and fatigue definitions and possible overlap require further comparison. The association cannot establish causation or genetic replication, and the published analytic sample does not establish how many participants have usable genotypes for this project.''')

before('### Supplementary-table legends', '''### NHANES supplementary materials

The accompanying `supplement/` directory contains the 49-item source and coding audit, survey model results, participant flow, item missingness by cycle and exposure, figure 1 and its source data, and the independent raw-GWAS score-source audit. The full original reconstruction package remains unchanged. See `audit_and_readiness.md` for source limitations, denominator and survey checks, and the genetic-score provenance assessment.''')

before('### Frozen Frailty Index estimates and workflow verification', '''### Exploratory CLSA physical-component follow-up

After the CLSA public PheWeb2 service became available, on 6 November 2026 we screened the 33 distinct previously selected prototype variants against European-ancestry, combined-sex marginal associations for average grip strength, maximum grip strength and four-metre walking time. Existing reciprocal GRCh37-to-GRCh38 position mappings were matched to the target by exact position, allele pair and rsID. All 99 variant–phenotype pairs were returned. Reported variant-level sample sizes were 22,809 for each grip measure and 24,231 for walk time, rather than the slightly larger catalog totals. The public methods describe rank inverse normal transformation, so effects are not in kilograms or seconds. The frontend labels ALT as the effect allele; full source genotype coding was not independently reconstructed. Five displayed p values were below 0.05, but none survived Holm correction across the complete 99-test family (all adjusted p values were 1). Rounded source p values were used for correction rather than ratios of rounded effects and standard errors.

These physical measurements are not FI, FFS or a joint sleep–frailty endpoint; the two grip outcomes overlap. The screen was exploratory, and the rs113851554 response was inspected during endpoint discovery before the full screen. Participant and family independence from every discovery/selection source remains unverified. The results do not test score calibration, discrimination or aggregate-score association, and absence of corrected marker-level significance does not establish absence of smaller effects. The public catalog had no code or description containing “frail” or “insomn”; this does not establish that the separate ordinal-frailty GWAS is available or that all potentially related labels have been exhausted. Complete results and source receipts are supplied in the [CLSA supplement](../clsa_public_followup_20261106/README.md).''')

# Keep the review copy self-contained for the explicit original-project links.
# Fail rather than quietly leave a broken link or omit a bibliography entry.
for relative in sorted(set(re.findall(r'\.\./(?:analysis|config)/[A-Za-z0-9_.-]+', text))):
    original = source.parent / relative
    destination = target.parent / ('original_' + relative[3:])
    if not original.is_file():
        raise FileNotFoundError(original)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(original, destination)
    text = text.replace(relative, destination.relative_to(target.parent).as_posix())
bibliography = source.parent / 'references.bib'
bib_text = bibliography.read_text()
cited = set(re.findall(r'@([a-zA-Z0-9_:-]+)', text))
entries = set(re.findall(r'@\w+\s*\{\s*([^,]+)', bib_text))
if missing := cited - entries:
    raise ValueError(f'Unresolved bibliography keys: {sorted(missing)}')
shutil.copy2(bibliography, target.parent / 'references.bib')
functional = target.parent / 'functional_followup'
functional.mkdir(exist_ok=True)
shutil.copy2(target.parent.parent / 'stage6_foxp2_gnomad_matched_coloc_review.md', functional / 'foxp2_gnomad_coloc_review.md')
target.write_text(text)
print('Created',target,'from',source)
