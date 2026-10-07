#!/usr/bin/env python3
"""Deterministic reports and claim/source ledgers from audited frozen outputs."""
import csv,json,hashlib,shutil
from pathlib import Path
from audit_evidence import ROOT,AREA,write_tsv
P=ROOT/'brain6/paper/final_package_v1';M=AREA/'MANUSCRIPT';M.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def put(p,s):p.write_text(s.strip()+'\n')
base=json.loads((AREA/'qc/baseline_verification.json').read_text());diag=json.loads((AREA/'qc/secondary_covariance_diagnostics.json').read_text())
put(AREA/'SCIENTIFIC_STATUS.md',r'''# Brain6 scientific status — 7 October 2026

**No validated new discovery emerged. Track A is NO-GO. Track B has a complete bounded manuscript for human scientific review; submission approval is on hold.** This continuation rechecks numerical tables, identifies prior regional reports, repairs source-free software validation, and preserves the failed historical endpoint. It does not rerun native GWAS/LD analyses.

| Component | Status | Quantitative evidence and limit |
|---|---|---|
| Global table | COMPLETE_VALIDATED at table level | 72 comparisons; 35 significant under inherited 396-test FDR; ADHD 9, MDD 10, SCZ 7, bipolar 8, PD 1, AD 0. Raw LDSC replay incomplete. |
| PLACO candidates | COMPLETE_EXPLORATORY | 25 pair candidates, 20 merged display regions; no candidate met higher validation classes. |
| Canonical LAVA | FAILED_QC | 17,465 cells; 13,745 TESTED, 3,720 NOT_RUN versus maximum 873. Immutable FAILED_QC_NOT_PROMOTED. |
| Source rescue | NOT_RUN full replacement | At least 2,847 recoveries required; insomnia 88-locus native-N pilot 47/41 versus 44/44, insufficient advancement. Long-sleep 1,291 failures alone exceed ceiling. |
| SUPERGNOVA | COMPLETE_EXPLORATORY secondary same-source analysis | 8,465 planned slots; 3,304 estimated; 61 low-SNP; 21 low-rank; 5,079 inapplicable. Three FWER-significant blocks; zero same-pair protected candidate overlaps. |
| Local correlation diagnostic | COMPLETE_VALIDATED arithmetic only | 3,304 covariance P values replayed with max error 3.33e-16. Of 2,199 numeric derived correlations, 723 exceed [-1,1]; 1,105 missing. One significant ADHD block has corr=1.1037; covariance wording only. |
| LD / SuSiE | FAILED_QC / BLOCKED_DATA | 0/32 matrices passed; genotype provenance/allele-order/valid reference reconstruction unavailable. No multi-signal PIPs, credible sets or coloc-SuSiE produced. |
| Single-signal ABF | COMPLETE_EXPLORATORY | 28/50 trait units estimated, 24 conditional credible sets. Six pair coloc estimates; three H4≥0.8 at default prior, one at lower prior. |
| Molecular QTL | COMPLETE_EXPLORATORY | Four ADHD-component tests; max H4=0.0824. No pair-level effector or molecular mechanism. |
| Positional enrichment | FAILED_QC full set; COMPLETE_EXPLORATORY subset | Two of 20 regions lack ten matched controls; separate 18-region/21-candidate sensitivity: zero significant of 54 tissues and 1,680 pathways. |
| Independent replication | BLOCKED_DATA / BLOCKED_ACCESS | No source-verified independent two-trait locus replication. Reused sleep, shared PGC cohorts and changed long-sleep definitions remain limited sensitivities. |
| Novelty | NEEDS_REVIEW for unsurveyed loci | chr5 shared region already published by Zu et al.; five of six ADHD intervals have prior-paper lead coordinates consistent with the interval. This is not five proven identical causal variants. |
| Native software | BLOCKED_SOFTWARE locally | R exists; LAVA, susieR, coloc, data.table and PLINK prerequisites absent. Installing software alone cannot repair absent genotypes/input contracts. |

The machine inventory uses COMPLETE_VALIDATED narrowly for table/receipt checks and conservatively NEEDS_REVIEW for uninspected implementation components. Its file enumeration is not a claim of line-by-line review of all 3,274 files. Review covered the named reports, final package, source-rescue and local-method protocols, numerical tables, frozen configurations, provenance, workflows, relevant historical commits and test failure paths. The broad atlas acceptance report is 3/23 and is separate from Brain6 readiness; no atlas completion is claimed.

Local retained archives were inspected before declaring gaps. The earlier v03 checkout contains large PLACO output files, continuous sleep data, and some raw sources; these are downstream/phenotype-specific and do not supply exact missing insomnia/ADHD inputs or genotype LD. Its raw/harmonized/reference links resolve to an unmounted Extreme SSD. None of the historical checkouts or unrelated frailty work was altered.
''')
put(AREA/'NOVELTY_AUDIT.md',r'''# Novelty audit — cutoff 7 October 2026

**There is no defensible first-discovery claim for Candidate C. No candidate currently supports Track A.** This is a date-bounded, documented search and targeted primary-source review, not an exhaustive systematic review. Five broad Europe PMC queries returned more than 1,000 results in four categories and were truncated at 1,000; 3,020 unique retrieved metadata records are retained. Focused title/genetics and sleep–psychiatric cross-trait queries retrieved all 269 and 20 results, respectively. Metadata/title screening is separate from full-text review. Searches omit unindexed records and cannot establish novelty from absence. Search and access receipts preserve query, date, limits, status and hashes.

| Prior work | Actual review | Implication |
|---|---|---|
| Jia et al. SLEEP 2025, DOI 10.1093/sleep/zsae209, PMID 39243390 | Primary indexed full text plus existing Brain6 comparison; seven sleep traits and three disorders; LDSC/HDL, GPA, PLACO, MAGMA/POPS/FUMA/deTS. Direct HTML download returned 403; supplemental locus table not fully acquired. | Psychiatric sleep sharing and cross-trait/gene prioritization precede Brain6. Existing 21 matched pair directions concordant with this work; neither new biological principle nor independent replication. |
| Xue et al. SLEEP 2026, DOI 10.1093/sleep/zsaf317, PMID 41065713 | PubMed/core abstract reviewed; insomnia Neff 314,149 and 12 disorders; MiXeR, conjFDR, ASSET, MAGMA, SEISMIC; 70 shared loci/97 candidate SNPs. Exact supplement loci unavailable in this session. | Large direct competitor with cellular claims. Unresolved exact loci preclude novelty declarations for unmatched Brain6 intervals. Similar Neff does not prove identical release or independent cohorts. |
| Zu et al. Transl Psychiatry 2026, DOI 10.1038/s41398-026-04166-4, PMID 42297780 | Primary main article, actual Table 3 and downloaded DOCX supplement inspected. Table 3: 14 insomnia–ADHD lead rows in 12 reported locus groups (one group has two leads at chr3 and one has two at chr5). | Five of six Brain6 ADHD candidate intervals contain a published lead at coordinate-consistent GRCh37 positions; two exact Brain6 lead rsIDs recur (rs6452785 and rs1476535). See machine crosswalk. No shared causal identity is inferred. |
| Schipper et al. preprint v1, DOI 10.1101/2025.10.18.25338281 | Abstract and indexed primary-source text; direct medRxiv HTML/PDF returned 403. Primary search-index table excerpt mentions INS–DEP chr11:112755447–113889019, local rg 0.87, P=1.66e-7, H4=0.01 and NCAM1. Full table context/build and subsequent versions not fully verified. | Candidate A already has a closely overlapping insomnia–depression regional report; this is a prior-art warning, not an allele-aligned equivalence or novel effector assignment. Preprint is not peer reviewed. |
| Lin et al. Mamm Genome 2026, DOI 10.1007/s00335-026-10232-5, PMID 42115503 | Primary subscription preview and PubMed abstract; 718 pleiotropic loci, 226 co-located loci; POPS/SMR and neurodevelopment interpretation. Full methods/supplement restricted. | Six-disorder breadth or neurodevelopment annotation alone is insufficient novelty. Exact locus comparison remains unresolved. |
| Fan et al. SleepChart 2026, DOI 10.1038/s43856-026-01656-w | Existing source-grounded Brain6 comparison and current indexed metadata. | Same UK Biobank with >8 h versus 6–8 h long sleep, distinct from ≥9 h versus 7–8 h; phenotype sensitivity, not replication. |
| Kozhemiako et al. 2025, DOI 10.1016/j.biopsych.2025.06.024 | Indexed review metadata and abstract. | Broad shared liability is established; demographic/source limitations remain relevant. |

Candidate C: Zu Table 3 reports rs77960 at chr5:103,964,585, PLACO P=6.54e-15, PP3=0.01, PP4=0.99. Brain6 tests chr5:103,447,968–104,447,968 with lead rs2431108, H4=0.9911/default and 0.9174/lower prior. Ensembl GRCh37 maps rs77960 to the published coordinate; Catalog 2026-10-03 reports GRCh38.p14 chr5:104,628,884 and ADHD risk allele G from Demontis 2023. Builds must not be mixed. The paper table does not list effect alleles or explicitly attest its build; GRCh37 consistency was independently checked. Different lead rsIDs are not independent discoveries without LD/allele evidence. The Zu insomnia source is a different UKB sample size (440,422 versus Brain6 386,533); this does not establish independent participants.

Other coordinate-consistent prior leads: chr3 rs62264764:117,639,575; chr5 rs6452785:87,685,500 and rs681446:88,054,292; chr7 rs1476535:114,071,035; chr11 rs12577364:28,658,030. The chr14 interval has no lead match in this particular table; it is not declared novel. Pair intervals are used, never merged geographic display spans.

Candidate A: insomnia–ADHD and insomnia–MDD share the SUPERGNOVA block chr11:112,459,489–114,257,728. Shared boundaries do not imply one causal variant, an ADHD–MDD independent signal, or a common effector. The long-sleep–SCZ protected candidate geographically lies in this area but belongs to a different pair; it cannot be counted as secondary validation. Correlated psychiatric liability, overlap/intercept uncertainty, LD and the out-of-range ADHD correlation are competing explanations. Prior Schipper INS–DEP regional evidence lowers the novelty expectation.

Candidate B: insomnia–MDD chr6:100,630,147–102,636,772 is family-significant covariance but not independently replicated or fine-mapped. Available searches do not settle novelty. No regional gene is claimed causal. Exact comparisons against Xue/Lin/Schipper supplements remain blocked or incomplete.

Ranking for future work: C has the strongest existing ABF stability and tractable interval but demonstrated low novelty; A has two secondary pair effects but a correlation diagnostic and prior regional warning; B has unresolved novelty and no multi-signal support. D, a psychiatric-versus-neurodegenerative difference, is currently confounded by unequal power, phenotype/source composition and selected tracks. No new contrast test was run; established psychiatric genetic sharing is not a new result. A bounded validation audit is the only completed defensible story.

GWAS Catalog API v2 was used because v1 was retired in August 2026. The live release is 2026-10-03, gene build GRCh38.p14, dbSNP156. Three rsID queries have 10/7/7 records and one page each; four competitor PMID queries have zero indexed associations. Those zeros do not indicate absent loci, absent papers or novelty. Derived fields and response checksums are in review/catalog_association_context.tsv and acquisition receipts.

Primary sources: [Jia](https://academic.oup.com/sleep/article/48/1/zsae209/7750695), [Xue PubMed](https://pubmed.ncbi.nlm.nih.gov/41065713/), [Zu Table 3](https://www.nature.com/articles/s41398-026-04166-4/tables/3), [Schipper v1](https://www.medrxiv.org/content/10.1101/2025.10.18.25338281v1.full), [Lin](https://link.springer.com/article/10.1007/s00335-026-10232-5), [Catalog API](https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/).
''')
put(AREA/'RESULTS_AND_VALIDATION.md',f'''# Results and validation

Table audit identity: {base['run_identity']}. Frozen protocol SHA256: {base['protocol_sha256']}. Starting public main commit: {base['starting_commit']}; current remote was checked and remained this commit.

The original package has 77 available SHA256 matches, one missing original canonical source and zero mismatches. The original package auditor was run and correctly fails on that missing source; its failure is retained. Our audit reports PASS_AVAILABLE_TABLES_ORIGINAL_INPUT_REPLAY_INCOMPLETE. Neither duplicate exported tables nor a new validator is a replacement for the missing provenance source.

Actual independent numerical checks:

- Recounted 72 global rows/35 inherited significant, 25 candidates/20 regions and all local-method statuses.
- Replayed all 3,304 SUPERGNOVA covariance P values from rho/sqrt(var); max absolute difference {diag['maximum_absolute_p_error']:.3g}. Recomputed Bonferroni family values, max error {base['max_fwer_error']:.3g}; recomputed same-pair protected interval overlaps: zero. Missing correlation values stay missing.
- Analytically reweighted all six H0–H4 vectors at lower p12; posterior agreement within 1e-10. This checks posterior arithmetic, not Bayes-factor derivation from raw inputs or the one-signal assumption.
- Recomputed BH in the original 54-tissue and 1,680-pathway families and empirical P=(extreme+1)/10001. Max BH difference <1.4e-12; no significant result. This checks archived permutation counts, not rerunning 10,000 null sets.
- Recovered source-first failure arithmetic: 2,847 recoveries needed; long-sleep alone leaves 1,291 NOT_RUN; unique smallest hypothetical perfect repair set is long sleep + insomnia + PD + MDD, leaving 533. No hypothetical recovery was admitted.
- Re-examined recorded LD factor/normalization audits. A diagonal reset without corresponding off-diagonal covariance normalization explains the mathematical scale failure in 19 diagnostic blocks. A normalized Gram matrix is PSD by construction; this does not validate genotype identity, alleles, cohort reference, or the originally exported LD. No empirical corrected LD or native SuSiE was run.

New diagnostic: 723 of 2,199 numeric local correlation estimates fall outside [-1,1]; 1,105 are absent. One significant chr11 insomnia–ADHD block has derived corr 1.103737. Its covariance test is retained with explicit inferential limits; no clipping or biological amplification claim.

Historical secondary results (GRCh37; same-source covariance, not replication):

| Pair | Interval | Covariance | Variance | P | Family-adjusted P |
|---|---|---:|---:|---:|---:|
| Insomnia–ADHD | chr11:112459489–114257728 | 0.0002691266 | 2.306485e-9 | 2.097239e-8 | 0.0001775313 |
| Insomnia–MDD | chr6:100630147–102636772 | 0.0001160579 | 5.911901e-10 | 1.813085e-6 | 0.0153477611 |
| Insomnia–MDD | chr11:112459489–114257728 | 0.0002255223 | 1.262670e-9 | 2.200532e-10 | 1.862750e-6 |

Software validation: 410 source-free contract tests pass with zero skips; 49 loaded test objects explicitly require archived sources/native runtimes (wildcard modules may contain more individual cases). The initial unpartitioned CI run failed 12 tests and errored on 39; this is an engineering availability problem, not evidence of scientific failure or success. Initial Brain6 pytest run: 10 failed, 534 passed, 3 skipped. Its ten real-data tests remain integration requirements; native tests remain separately labeled. The synthetic zero-family fixture now mocks available disk space; production disk checks and the insufficient-space test remain unchanged. No original assertion was removed. Python compilation, shell syntax, panel lock validation, phase0 strict declaration audit and synthetic end-to-end smoke pass. Phase0 reports 0/45 real harmonization/LDSC-ready traits in the new checkout; synthetic smoke has no biological meaning. Broad atlas report: 3/23 gates.

See qc logs and machine test partitions for final numerical test counts, CI execution status and unresolved native checks. No failed native analysis was recast as a pass. The original full command remains `python -m unittest discover -s tests -v`; it must be run after exact archives/runtimes are restored. Source-free CI is a separate named gate.
''')
put(AREA/'CRITICAL_REVIEW.md',r'''# Automated adversarial scientific review

This is one agent's structured hostile review, not independent human peer review or a panel of actual specialists.

| Objection | Attempted falsification | Disposition |
|---|---|---|
| Global findings may be inherited/known | Recounted original FDR labels and compared psychiatric papers | Supported only as a descriptive 72-row map; no new 72-test correction or first-report claim. |
| Candidate C is already known | Retrieved actual Zu Table 3; matched rs77960 GRCh37 and Catalog build/allele | First-region-discovery claim rejected. Different lead rsID is insufficient novelty. |
| Candidate A is one causal signal across disorders | Same-source blocks overlap; prior Schipper region and unstable ADHD derived correlation | Shared causal identity rejected; independent, allele-aligned multivariate work required. |
| Failed LAVA is a null | Recounted all statuses; original ceiling 873 | Rejected. 3,720 cells not estimated; no bivariate promotion. |
| Pilot rescues whole family | Paired 88-locus improvement and seven-trait arithmetic | Rejected. Pilot below fixed advance alternatives; long-sleep alone exceeds ceiling. |
| Secondary family cherry-picked | Recomputed 0.05/8,465 including inapplicable tracks | Correction verified; no denominator shrinkage. |
| Merged regions validate protected candidates | Computed same-pair overlap rather than geographic merge | Zero supported candidates; long-sleep–SCZ chr11 geographic overlap does not validate insomnia pairs. |
| Correlation may be numerically invalid | Replayed covariance tests; counted missing/out-of-range correlations | Covariance P arithmetic passes; correlation claim restricted. Human statistical interpretation remains required. |
| LD can be repaired by clipping/PSD projection | Read scale/Gram diagnostic, original failure gate, available reference inventories | No valid reference/allele reconstruction. Multi-signal claims remain blocked. |
| High ABF H4 proves causal variant | Reweighted six posterior vectors and retained full H0–H4 | Three default-support intervals become one at lower prior; all conditional on single-signal model. |
| QTL/nearest gene proves mechanism | Reviewed four actual component tests and tested universes | Max H4 0.0824; no pair effector; annotations are coordinate context. |
| Removing hard enrichment regions improves outcomes | Preserved full-set gate and separately labeled 18-region sensitivity | Zero FDR-significant tissues/pathways; no full-set inference. |
| FinnGen/MDD noUKB are independent replication | Read phenotype and cohort ledgers | FinnGen reuses sleep; MDD retains nine named PGC cohort overlaps. At most outcome-only/partial sensitivity. |
| Alternative long sleep is same exposure | Compared ≥9h/7–8h with ≥10h/7–8h and >8h/6–8h | Phenotype sensitivity only; original source model not repaired. |
| CI makes scientific archive complete | Ran original auditor and real source tests before partition | Missing receipts remain visible. Source-free pass cannot certify empirical results. |
| Full text gaps imply novel loci | Logged denied/missing supplements, truncation and index zeros | Rejected. Unmatched A/B/chr14 novelty unresolved. |

No unsupported claim remains in Abstract or Conclusions. Regression tests cover BH ties versus SciPy, posterior odds, pair-preserving overlap and numerical/claim validators; full native pipeline validation is blocked, not simulated.

Questions for a real collaborator: Is the secondary covariance estimator adequately calibrated under these meta-analysis N/intercept approximations? What explains 723 out-of-range correlations and the significant ADHD example? Are retrospective negative-validation lessons sufficiently generalizable for publication? Can allele-aligned, participant-independent releases be obtained for both traits? Can archived genotypes and RDS/factor provenance establish a valid replacement LD? Which prior-paper signals are actually equivalent after build/allele/LD alignment? These questions affect inference and editorial merit and are not resolved by software tests.
''')
put(AREA/'PUBLICATION_READINESS.md',r'''# Publication readiness — two separate tracks

Track A, stronger original discovery: **NO-GO** for Genome Medicine, Molecular Psychiatry and a specialist discovery article. No new independent locus validation, valid multi-signal result or confirmed effector was produced. The strongest ABF interval has prior regional evidence. Novelty for other loci is unresolved rather than established. A completed manuscript does not close these scientific gaps.

Track B, bounded secondary-validation paper: **GO for human scientific review and submission preparation; NO-GO for submission approval now.** The editable manuscript, numerical sources, vector figures, supplementary provenance, references, reporting checklist and unsigned cover letter are complete. The evidence is a case-study audit of validation limits, not a general benchmark proving that methods fail. Authors must judge whether this contribution warrants an original research paper, restore/verify the missing original canonical source, reconcile interpretation with an independent collaborator, and complete ethics/authorship/funding/COI/licensing declarations.

| Journal | Current editorial assessment | Required next decision |
|---|---|---|
| Genome Medicine | NO-GO as a strong original discovery. The Mental health and neuropsychiatric disorders collection is topically relevant but does not reduce scientific novelty/validation standards. Collection open, deadline 12 May 2027 as checked 7 Oct 2026. | An investigator must substantiate a meaningful new insight or methodological contribution beyond this single archive audit. |
| Molecular Psychiatry | NO-GO as presently supported. Broad psychiatric sharing, positional genes and exploratory single-signal coloc have close competitors. | Stronger independent evidence and a clear biological or methodological advance required; a case-study audit alone has uncertain editorial fit. |
| BMC Medical Genomics | Conditional specialist fallback for Track B, subject to human scientific review and contribution assessment. Its genomics/medical scope is relevant. | Determine whether the secondary-analysis contribution and reproducibility meet the journal's standard; no guarantee or acceptance probability. |

Current official guidelines: [Genome Medicine research](https://link.springer.com/journal/13073/submission-guidelines/research) requires structured abstract ≤350 words, 3–10 keywords, Background/Methods/Results/Discussion/Conclusions and full Declarations. [Collection](https://link.springer.com/collections/dbjjedaefg). [Molecular Psychiatry preparation](https://www.nature.com/mp/authors-and-referees/preparation-of-articles) currently lists Articles up to 5,000 words, unstructured abstract 150–250 words, six display items and 100 references; linked article-type PDF lists 3,500 words/five items/75 references. We use conservative 3,500/five/75 pending author confirmation of this inconsistent guidance. [BMC Medical Genomics scope](https://link.springer.com/journal/12920/aims-and-scope). Actual retrieved HTML hashes are retained where acquisition succeeded.

Manuscript authorship and institutional statements remain explicitly unverified. AI is not an author. Neither collection suitability nor a compiled PDF is submission authorization. No manuscript has been submitted and no third party has been contacted.
''')
put(AREA/'BLOCKERS_AND_ACTIONS.md',r'''# Terminal blockers and exact restart actions

1. **Original archive/provenance**: mount the authorized Extreme SSD or restore its exact archival contents. Required canonical source: work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv; compare against final_package_v1 provenance, then rerun original package audit. Restore exact insomnia/ADHD dense raw files, harmonized sources, 2,495-locus file, reference genotypes and production receipts. Retained downstream PLACO files are not replacements. Do not repoint frozen paths to unrelated files or alter historical checksums.
2. **Long-sleep eligibility**: request source authors' exact Dashti ≥9 h versus 7–8 h export contract: model, case/control encoding, effect units/A1 orientation, variant-specific analyzed N or validated constant-N justification, effective-N mapping, missingness and release identity. Current absence of released N/model attestation blocks a full new LAVA experiment. Insomnia native N is documented but its pilot cannot rescue the family; MDD effective-N and PD cases/proxies still need exact source semantics. ADHD/bipolar/SCZ cohort conventions remain source-specific. Do not substitute an adjacent sleep definition.
3. **Independent data**: FinnGen current access requires its official form/terms; candidate object HTTP403 was excluded. Request authorized release-specific ADHD/MDD and clinically comparable insomnia sources, population/ancestry strata, per-variant beta/logOR, SE, P, A1/A2, frequency, INFO, N and case/control fields, build, complete coverage, licenses and cohort membership. Exact sleep equivalence and both-trait independence must be documented before any claim. No form or message was submitted.
4. **LD repair**: restore allele-ordered genotypes and reference variant/sample manifests. Recompute signed Pearson LD with verified sample count/ancestry, build, frequency, reference A1/A2, SNP order and summary orientation. Diagnose covariance scaling separately. Require symmetry/unit diagonal/PSD/conditioning and RSS consistency. No clamp/projection promotion. Install pinned LAVA/susieR/coloc/PLINK only after viable data/resources are confirmed, then freeze new run identity before outcomes.
5. **Novelty full texts**: acquire authorized Xue and Lin supplements and full Schipper table/build context; inspect all candidate leads/alleles/regions and later preprint versions. Verify competing locus identities using source build and LD, not current Catalog GRCh38 coordinates directly against Brain6 GRCh37. Access restrictions remain genuine.
6. **Human review**: statistical geneticist review of covariance calibration/correlation anomalies and cohort overlap; investigator assessment of negative-validation contribution. Complete author names/affiliations, CRediT roles, ethics/consent basis, funding, COI, data licenses and final wording. These are prerequisites for submission approval.

Prepared investigator request (unsigned, not sent):

> We are evaluating sleep–psychiatric genetic sharing in an explicitly secondary-validation analysis. Please provide or identify the exact release of [phenotype/source/release selected by investigator], its phenotype thresholds and exclusion rules, cohort memberships and overlaps with UKB/PGC/iPSYCH/deCODE/FinnGen, ancestry strata, association model and effect units, genome build, A1/A2 orientation, imputation/reference metadata, and variant-level rsID/CHR/BP/beta or logOR/SE/P/EAF/INFO/N plus case/control counts where applicable. For Dashti long sleep we require ≥9-hour cases and 7–8-hour controls and the mapping of exported effect/sample-size fields to the fitted model, rather than a changed threshold. Please provide genome-wide tabular compressed files, version/date, SHA256 and permitted research/redistribution terms. We will not represent outcome-only or cohort-overlapping analysis as independent two-trait replication.

No large GWAS acquisition is authorized by availability alone. A restart must record access, exact bytes, hashes, storage headroom for intermediates, and source prerequisites. Current analyses are terminally blocked pending these inputs; bounded manuscript preparation is complete.
''')
# Copy all numerical final tables and provenance for a standalone, navigable handoff.
S=M/'source_tables';S.mkdir(exist_ok=True)
for p in P.glob('*.tsv'):shutil.copy2(p,S/p.name)
for n in ['BRAIN6_FINAL_PROVENANCE.json','BRAIN6_FINAL_AUDIT.json']:shutil.copy2(P/n,S/n)
for n in ['significant_secondary_blocks.tsv','coloc_prior_replay.tsv']:shutil.copy2(AREA/'qc'/n,S/n)
# Carry original exhaustive source fields, with explicit present availability and replication semantics.
ledger=[]
for r in rows(ROOT/'brain6/manifests/gwas_external_availability.tsv'):
 ledger.append({**r,'continuation_access':'ORIGINAL_ARCHIVE_UNAVAILABLE_IN_CURRENT_CHECKOUT','replication_class':'INELIGIBLE_AS_ITS_OWN_REPLICATION','current_hash_verification':'HISTORICAL_RECEIPT_ONLY_NOT_RAW_REHASHED'})
fields=list(ledger[0])
extras=[('FinnGen ADHD R13','INDEPENDENT_ONE_TRAIT_ONLY','Official form/terms; reused original insomnia; full participant independence unverified','https://www.finngen.fi/en/access_results'),('PGC MDD2025 noUKB','PARTIAL_COHORT_OVERLAP','Nine named PGC cohort overlaps with Howard2019; reused sleep','https://pgc.unc.edu/for-researchers/download-results/'),('Fan2026 SleepChart long sleep','PHENOTYPE_SENSITIVITY','UKB >8 h vs6–8 h; Ncase25049/Ncontrol300420; not original phenotype','https://doi.org/10.1038/s43856-026-01656-w'),('Austin-Zimmerman2023 long sleep','PHENOTYPE_SENSITIVITY','UKB+MVP ≥10 h vs7–8 h; UKB overlap','https://sleep.hugeamp.org/downloads.html'),('HUNT/23andMe independent sleep candidate','UNVERIFIED','No authorized complete compatible summary source acquired; exact cohorts/license/power unresolved','https://cncr.nl/research/summary_statistics/')]
for name,cl,notes,url in extras:
 x={k:'UNVERIFIED_OR_NOT_APPLICABLE' for k in fields};x.update(trait=name,study=name,source_URL=url,source_page_URL=url,replication_class=cl,continuation_access=notes,current_hash_verification='NO_NEW_GENOME_WIDE_OBJECT_ACQUIRED',current_readiness='BLOCKED_ACCESS_OR_ELIGIBILITY');ledger.append(x)
write_tsv(AREA/'SOURCE_AND_ACCESS_LEDGER.tsv',ledger)
claims=[]
for ident,file,method,result,allowed,failure in [
 ('global','BRAIN6_FINAL_GLOBAL.tsv','Inherited LDSC table FDR','72 rows;35 inherited significant','The selected display contains 35 inherited FDR-significant entries.','Raw/log replay incomplete; no new multiplicity family'),
 ('candidates','BRAIN6_FINAL_CANDIDATES.tsv','Frozen PLACO candidate rule','25 candidates;20 display regions;0 promoted','These are screen candidates without independent causal-locus confirmation.','Candidate status not a biological null'),
 ('lava','BRAIN6_FINAL_LAVA.tsv','Canonical frozen eligibility QC','13745 tested;3720 notrun;max873','The canonical family failed QC and was not promoted.','Never claim tested null or successful rescue'),
 ('secondary','BRAIN6_FINAL_ALT_LOCAL_BLOCKS.tsv','SUPERGNOVA;5*1693 family','3 significant;0 same-pair candidate overlaps','Three secondary covariance blocks passed the original family threshold.','Same discovery sources; N/intercept approximation; derived correlation instability'),
 ('coloc','BRAIN6_FINAL_TRAIT_COLOC.tsv','Single-signal ABF with lower prior sensitivity','3 default;1 lower prior H4>=0.8','One interval retained descriptive single-signal support under the lower prior.','Unverified one-signal assumption; invalid LD; prior regional publication'),
 ('qtl','BRAIN6_FINAL_EQTL_COLOC.tsv','Component GWAS QTL coloc','4 combined e/sQTL tests;maxH4 .0824','No pair-level effector mechanism was established.','Component-only;single-signal;not proof of molecular absence'),
 ('enrichment','BRAIN6_FINAL_PATHWAYS.tsv','18-region matched subset;BH','0/1680 pathways;0/54 tissues','No tissue or pathway survived its original subset FDR family.','Full20 gate failed;positional not causal-weighted'),
 ('replication','BRAIN6_FINAL_REPLICATION.tsv','Source/cohort audit','No independent two-trait locus replication','Independent two-trait locus replication was not established.','Reused sleep and cohort overlap; phenotype mismatch')]:
 claims.append(dict(claim_id=ident,exact_input=str((P/file).relative_to(ROOT)),input_sha256=sha(P/file),method=method,result=result,independent_evidence='NONE_AT_TWO_TRAIT_LOCUS_LEVEL',failure_conditions=failure,allowed_wording=allowed))
write_tsv(AREA/'CLAIM_TO_EVIDENCE.tsv',claims)
put(M/'SUPPLEMENTARY_METHODS_AND_RESULTS.md',r'''# Supplementary methods and results

The main manuscript is a secondary table/provenance audit. Historical native analyses are described from frozen protocols/receipts; this continuation independently replays numerical arithmetic only. Full native input replay remains incomplete.

S1. source_tables/BRAIN6_FINAL_GLOBAL.tsv has all 72 global estimates, standard errors, original q-values, intercepts, phenotype/sample/source fields. Original family size is 396, not 72. Six displayed disorders span psychiatric and neurodegenerative outcomes; detection counts cannot be compared as a powered architecture contrast.

S2. BRAIN6_FINAL_CANDIDATES.tsv and REGIONS.tsv preserve 25 pair identities and 20 display merges. Five screen tracks: insomnia–ADHD, insomnia–MDD, long sleep–bipolar, long sleep–PD, long sleep–SCZ. PLACO_PLUS with protected legacy insomnia–ADHD lineage; P<1e-8 is the predeclared five-track lead threshold (5e-8/5). Greedy clumping: LAVA UKB EUR v1.1 reference, 100,000 samples, r²=0.1, 500kb; only exact rsID/CHR/BP reference matches. Lead flanks ±500kb merge within pair; across-pair merges are display only. This locus rule was frozen after initial PLACO outputs and before annotation; it must not be described as fully pre-discovery.

S3. BRAIN6_FINAL_LAVA.tsv contains 17,465 trait–locus cells across seven traits and 2,495 loci. Canonical local-h² gate and complete family QC precede bivariate inference; historical decision remains failed. NOT_RUN reasons: 3,564 low-local-h², 154 fewer-than-minimum shared reference variants, two fewer-than-minimum K. By trait: long sleep1,291; insomnia671; PD636; MDD589; ADHD223; bipolar168; SCZ142. No new threshold or source substitution was used. Native-N insomnia pilot 88 loci: 47/41 tested/notrun versus44/44;3.41percentage-point gain and6.98%relative low-h² reduction; advancement required10pp or25%. No material improvement.

S4. ALT_LOCAL_BLOCKS.tsv contains all8,465slots, REGIONS/VALIDATION tables all candidate states. Historical SUPERGNOVA reference: 503-person1000Genomes EUR GRCh37, common unique nonpalindromic rsIDs, genetic map, per-chromosome reconstruction following recorded failed cM-reset attempt. Minimum shared SNPs120. Signed Z/study-total N were admitted only for two insomnia pairs as limited meta-analysis approximations; three long-sleep pairs inapplicable before results. MHC excluded;1,693blocks/pair;fullfive-pairBonferroni. Fourworkers/eligiblepair sequentially. The patched implementation and source hashes are in historical run_provenance.json. Continuation Wald replay usesP=2Φ(-|rho|/sqrt(var));correlation missingness/out-of-range diagnostic is retrospective and does not alter this family.

S5. FINE_MAPPING.tsv:50candidate–traitunits,28ABFestimates,24single-signal95%credible sets,4unsupportedposteriorassociations,9input-QCholds,13long-sleepmodelholds. Credible-setsize3–1,445,median53. These are model-conditional outputs from historicalcoloc5.2.3, notnewSuSiEPIPs. All32assessedLDmatricesfailed;nosuccessfulnativeSuSiE. TRAIT_COLOC.tsv andPRIOR_SENSITIVITY.tsv includeH0–H4,sixestimatedintervals,notzerosfor19unestimatedpairs. Fixedp1=p2=1e-4,p12=1e-5 with1e-6 and5e-5sensitivities. Posteriorreweight: multiplyH4bypriorratio,renormalizeallfivecomponents. ValidateoriginalserializedposteriorsratherthanreconstructingrawBayesfactors.

S6. EQTL/SQTL/REGULATORY tables retaincomponenttestandannotationstatus. Fourrelease7GTExbraincontexts: cortex/frontalcortex eQTL andLeafCutter sQTL. GRCh37windowsliftedtoGRCh38;position/unorderedallelesmatched. Atleast500sharedvariants,QTLandcomponentGWASP≤5e-8fornumericaladmission;QTLN=AN/2andestimatedsdY. FouractualADHDcomponenttests;maximumH4.0824. No multi-trait molecular claim, cell-specific effector or experiment. Published QTL credible-set membership and six exact-lead regulatory overlap rows are coordinate evidence only.

S7. TISSUE_CELLTYPE/PATHWAYS tables retain54GTExv8bulk-tissuetests/1,680Reactomepathwaysandallfailures. Full20-regionmatchedcontrolsgatefailedtworegions(4and7anchors;minimum10). Separatelyfrozen18-region/21-candidatesensitivityused10,000acceptednonoverlappingmatchednullsets. ReplayedP=(extreme+1)/10001andBHwithinseparatecompletefamilies. Minimumq.113388661134tissue/.2426424024267pathway. No significantsubsetfinding,nofullsetP,noregulatory/cellresolvedenrichment. Annotating nearest genes is not a causal weighted universe.

S8. REPLICATION.tsv and source/access ledger distinguish FinnGen outcome-onlyglobalcontext(rg.3817,SE.0521,P2.327e-13; original insomnia reused), PGC2025noUKBMDD(rg.4771,SE.0228,P4.65e-97;ninenamedPGCoverlapentries), andUKB-overlappingadjacentlongsleep. These are historical table results, not newly rerun independent replications. No locus-leveltwo-traitreplication exists. No MR was run in this continuation.

S9. NOVELTY_AUDIT.md, search manifests/receipts and coordinate crosswalk preserve actual access and uncertainty. CatalogAPIv2release2026-10-03GRCh38.p14cannotbedirectlycomparedtoGRCh37. Originaloutput/nativeinputhashesremainseparatefromnewtablehashverification. The missingcanonicalsourceisexplicitin77MATCH/1MISSINGchecks. Copyrightedrawliterature/fullmetadataresponsesareignoredscratch; derivedfacts and primary links are delivered.

S10. Reproduction commands are in REPRODUCE.md. Source-free tests exercise artificial fixtures/numerical arithmetic; archived-data and native tests remain a separate unexecuted integration gate. Deterministic figures use only audited numerical inputs. All dates, protocol identity and script hashes are recorded; no post hoc positive family or synthetic empirical output was admitted.
''')
# Improve readability of deliberately compact scientific variable counts in supplement.
p=M/'SUPPLEMENTARY_METHODS_AND_RESULTS.md';s=p.read_text();s=s.replace('all8,465slots','all 8,465 slots').replace('allcandidate','all candidate').replace('all32','all 32');p.write_text(s)
put(M/'DECLARATIONS_AND_AUTHOR_FIELDS.md',r'''# Human completion fields — no approvals inferred

Authors / affiliations / corresponding author: UNVERIFIED; investigators complete.
CRediT roles: investigators assign Conceptualization, Data curation, Formal analysis, Methodology, Software, Validation, Visualization, Writing–original draft/review, Supervision, Funding acquisition. The automated agent is not an author.
Ethics approval and consent: investigator confirms the secondary use of the exact public/licensed summary statistics, original studies' consent/ethics, and whether a local exemption applies. No invented approval number.
Consent for publication: investigator determines applicability; no individual-level content is presented.
Funding and grant numbers: UNVERIFIED; human completion required.
Competing interests: UNVERIFIED; every author must declare.
Acknowledgments: confirm source consortium/reference acknowledgments and permissions; no invented collaborator.
Data/code availability: new branch contains protocols, derived numerical source tables, validators and source hashes; no restricted raw source redistribution. Exact licensed sources must be obtained from source providers under their terms. One canonical original provenance source and native inputs are presently inaccessible; full native reproduction is on hold. Investigators select a permanent archive/version DOI after review.
AI disclosure draft: OpenAI Codex was used for code/reproducibility audits, numerical table checks, literature metadata screening, manuscript drafting and figure generation. A responsible human investigator must verify all factual/statistical statements, references and final interpretations. Automated review is not independent peer review. No AI authorship is claimed.
Author approval/originality/exclusive submission: NOT PROVIDED. This package is not submission-approved and has not been sent to a journal.
''')
put(M/'COVER_LETTER_DRAFT.md',r'''# Unsigned cover letter for investigator adaptation

Dear Editors,

Please consider, after investigator approval, “Limits of locus validation in a six-disorder sleep genetic atlas: a reproducible secondary analysis” as a secondary-validation research article. We examine how evidence progresses from a 72-comparison global map and 25 pair-specific screening candidates through frozen local-QC gates, secondary covariance, exploratory colocalization and molecular/replication checks.

The manuscript reports the failed primary LAVA gate, all inapplicable secondary slots, three corrected covariance blocks that do not validate protected same-pair candidates, prior-sensitive single-signal results and the absence of independent two-trait locus confirmation. It also documents numerical diagnostics and prior regional reports. We do not claim a new mechanistic discovery or that a failed QC gate demonstrates biological absence.

For Genome Medicine, the subject relates to the Mental health and neuropsychiatric disorders collection; investigators must first establish whether its bounded contribution meets the journal's original-research standard. For a specialist fallback, emphasize the reproducibility and validation case study rather than expanding the discovery claim. Scientific review and the missing original-source verification remain required before sending this letter.

Authors must supply corresponding-author details, funding/COI/ethics declarations, author approval, exclusive-submission confirmation and any suggested reviewers. None of those confirmations is implied here.

Sincerely,
[Human corresponding author to complete; unsigned]
''')
put(M/'JOURNAL_SUBMISSION_CHECKLIST.md',r'''# Journal-specific submission preparation checklist

Genome Medicine: structured abstract≤350words;3–10keywords;title/authors/affiliations/correspondence;Background/Methods/Results/Discussion/Conclusions;abbreviations;Declarationswithethics,consent,data,COI,funding,contributions,acknowledgments. Main source follows this structure. Collection deadline12May2027 checked7Oct2026; topical fit does not imply editorial readiness.

Molecular Psychiatry: alternative unstructured abstract150–250words inMP_ABSTRACT.txt;conservative≤3,500mainwords,fivedisplayitemsand75refs due live-guideline/linkedPDFconflict. Fourmainfiguresandoneprincipalresults table. FiguresseparatePDF/SVG,legendsafterreferences; supplementseparate. Humanauthorsconfirmcurrentlimitsonsubmissionday. No author/ethics/COI claims are auto-filled.

BMC Medical Genomics fallback: scope medically relevant genomic analysis; confirm current research article requirements and investigator acceptance of secondary-validation framing. No automated submission.

All journals: verify raw-source audit hold; independent collaborator reviews covariance/correlation anomalies and model assumptions; allauthorsreviewandapprove;complete declarations;checksources/licensing/permanentarchive;removeinternalreviewstatusonlyafterapproval;verifyfigurelegibility/references/reportingchecklist;ensureallnegative/unestimatedstatesvisible. OriginaldiscoveryTrackAisNO-GO.
''')
# 22-item STROBE/STREGA mapping: adaptation rather than fabricated participant-level reporting.
items=[('1','Title/abstract','Title and structured abstract','Secondary table/provenance analysis explicit'),('2','Background/rationale','Background','Prior sharing established; validation question specified'),('3','Objectives','Background final paragraph','Evaluate validation and reproducibility, not causal discovery'),('4','Study design','Methods design','Retrospective continuation; historical frozen protocols separately dated'),('5','Setting','Methods sources; source ledger','Original GWAS setting/cohort dates must be verified by investigators'),('6','Participants','Source ledger; Supplement S1/S8','No new participants; original inclusion/selection and ancestry are source-specific'),('7','Variables','Methods sources; ledger','Phenotypes, build/effects and changed sleep thresholds explicit'),('8','Data sources/measurement','Methods; source ledger','Exact releases/hash fields; incomplete allele/N contracts disclosed'),('9','Bias','Discussion limitations; critical review','Overlap, power, selection, LD and available-source bias'),('10','Study size','Methods families; Results','Fixed counts and missingness; no new prospective power claim'),('11','Quantitative variables','Methods; S3–S5','N, covariance, posterior priors and boundaries retained'),('12','Statistical methods','Methods; Supplement S2–S7','Inherited FDR, full Bonferroni family, prior/matched-null sensitivity'),('13','Participants/results flow','Results; Figure2','Trait–locus/slot flow replaces participant flow; unavailable original individual flow explicit'),('14','Descriptive data','Global source table and source ledger','Age/sex/individual covariates not uniformly available; no imputation'),('15','Outcome data','All numerical source tables','All negative, missing and inapplicable outcomes'),('16','Main results','Results; principal table; Figures1–4','Effects/SE/P/q and posterior context; no clinical risk prediction'),('17','Other analyses','Sensitivity Results; S3–S8','Pilot, lower prior and positional subset separate'),('18','Key results','Discussion opening','Failed gate and limited validation, no discovery'),('19','Limitations','Discussion limitations','Missing raw replay, model validity, independence, novelty access'),('20','Interpretation','Discussion/Conclusions','No causal gene/null-from-QC/replication overclaim'),('21','Generalizability','Discussion limitations','EUR/meta-analysis/source availability; cannot compare disease power as architecture'),('22','Funding','Declarations template','UNVERIFIED; human completion required')]
write_tsv(M/'STROBE_STREGA_CHECKLIST.tsv',[dict(item=i,reporting_domain=d,manuscript_location=l,status_or_limitation=s) for i,d,l,s in items])
put(M/'STREGA_ADDITIONAL_ITEMS.md',r'''# STREGA-specific genetic reporting

State first-report/replication status: secondary validation; prior regional evidence, no independent two-trait replication.
Genotyping/imputation/calling and error: no new genotyping; original source/provider methods and exact QC fields must be checked. Missing native genotype/allele-order provenance prevents new multi-signal work.
Population stratification: European reference/source emphasis; original LDSC intercepts retained, not used to infer participant independence.
Relatedness: source-specific original exclusions; not assessed de novo from summary statistics.
Hardy–Weinberg equilibrium: no new participant genotypes; original studies/ref QC must be verified, not invented.
Haplotype/LD modeling: documented reference, build, SNP order, signed alleles and invalid-matrix gate; no silently repaired LD.
Replication and multiplicity: inherited396FDR, five-trackPLACO threshold, complete8,465secondaryslots and separate54/1,680BHfamilies; partial sensitivities cannot be independent replication.
Data availability: derived tables/code available; original native replay/source access holds disclosed. STREGA checklist is an adaptation for secondary summary analysis, not certification of missing original participant-level details.
''')
put(AREA/'REPRODUCE.md',r'''# Reproduction and execution boundary

From the new branch root, with Python3.11+:

```
python brain6/confirmatory_v2/scripts/audit_evidence.py
python -m pytest -q brain6/confirmatory_v2/tests
python brain6/confirmatory_v2/scripts/run_contract_tests.py --suite source-free
PYTHONPATH=extensions/brain6 python brain6/confirmatory_v2/scripts/run_brain6_tests.py
python brain6/confirmatory_v2/scripts/validate_manuscript.py
python brain6/confirmatory_v2/scripts/build_figures.py
```

Continuation dependencies: numpy,pandas,scipy,matplotlib,pytest,beautifulsoup4; recorded actual versions in qc/environment.json. Production requirements remain separately pinned in original requirements files. No local native R/PLINK production run is claimed. Original archived-data integration: `python brain6/confirmatory_v2/scripts/run_contract_tests.py --suite integration`; full unchanged suite: `python -m unittest discover -s tests -v`. Full Brain6 integration requires restored sources: `PYTHONPATH=extensions/brain6 python -m pytest brain6/scripts/tests extensions/brain6/tests`.

Optional bounded public acquisition:

```
python brain6/confirmatory_v2/scripts/acquire_literature.py
python brain6/confirmatory_v2/scripts/reacquire_manifests.py --manifest brain6/confirmatory_v2/review/focused_acquisition_manifest.tsv
# Explicit opt-in network replay; public metadata only, ≤manifest byte cap/object:
python brain6/confirmatory_v2/scripts/reacquire_manifests.py --execute --manifest brain6/confirmatory_v2/review/reference_acquisition_manifest.tsv
python brain6/confirmatory_v2/scripts/diagnostic_extension.py
```

The coordinate/metadata extension requires local raw public literature cache acquired under its manifests. Without it, use shipped derived tables and receipts; do not fabricate missing responses. Raw responses are excluded from Git and handoff. Checkpoint receipts are atomically updated after each acquisition; successful objects have SHA256. Downloads are bounded30MB orlower; querytruncationdeclared. Inventory hashes stream1MB blocks; numerical audits process small tables and avoid loading retained large GWAS files. Figures are deterministic, seeded nowhere because no stochastic operation is run. Matplotlib usesAgg; fixed metadata date removed fromPDF to permit byte-stable figures. Native analyses remain blocked rather than launching expensive downloads.

PDF: edit the standaloneMANUSCRIPT/brain6_manuscript.tex in the Codex built-in editor; use compile_latex_document for diagnostics. Terminal export via installedTectonic: `tectonic --outdir brain6/confirmatory_v2/MANUSCRIPT brain6/confirmatory_v2/MANUSCRIPT/brain6_manuscript.tex`. Figures are supplied separately as vectorPDF/SVG with numerical source tables. Original source package remains unchanged; manifests bind both old/new source hashes. Do not run inventory.py on unrelated mutable checkouts expecting identicaltimestamps; inventories are session evidence.
''')
put(AREA/'FINAL_HANDOFF.md',r'''# Final Brain6 handoff

**Track A: NO-GO. Track B: complete submission-preparation package, ready for human scientific review; not ready for submission approval.**

Completed: current repo/branch/archive inventory, historical table/hash/numerical audits, source/access and claim ledgers, literature metadata searches with targeted primary paper/supplement checks through7October2026, candidate novelty crosswalk, explicit integration/source-free CI partition and deterministic manuscript/figures/supplement/checklists. New branchbrain6/continuation-2026-10-07;main and all historical scientific outputs preserved.

Independent checks reproduce72/35global display counts,25/20candidate/region counts,failed17,465-cellLAVA decision,complete8,465secondaryfamily andzero same-paircandidateoverlaps,sixABFpriorreweightings,54/1,680enrichmentBH/permutation arithmetic and77availablehashes. One originalcanonical source is missing; the original auditor correctly fails. No nativeGWAS/LD replay or independentbiologicalreplication is claimed.

New audit findings: significantchr11ADHDderivedcorrelation1.1037requirescovariancewording;723/2,199numericcorrelationsareoutofbounds. Zu2026reportschr5rs77960sharedregionwithPP4.99;first-regionnoveltyrejected. Five of six ADHD intervals have coordinate-consistentpriorleads;twoexactleadrsIDsmatch. Closely overlapping INS–DEP priorpreprint evidence also weakensCandidateA. These are validation/novelty findings, notnewbiologicaldiscoveries.

LAVA remainsFAILED_QC_NOT_PROMOTED;SuSiE remainsfailedLD/blockedreference;single-signalABF remainsconditional;noeffector,independenttwo-traitreplication orFDR-significantpositionalenrichment. NOT_RUN/METHOD_INAPPLICABLE are not null findings.

Supported claims: an inheritedglobalmap, an explicitcandidateuniverse, a failedprimarygate, threecorrectedsecondarycovarianceblocks withoutsame-paircandidatevalidation, prior-sensitiveexploratorycoloc and boundednegative/subsetresults. No shared causal variant, new mechanism, cross-disease causal effect or universally failedmethod claim is supported.

GenomeMedicine/MolecularPsychiatry: NO-GO for a strong discoverysubmission; Genomecollectiontopicallyrelevantbutscientificbarunchanged. BMCMedicalGenomics is a conditional specialist fallback for the boundedstoryafterhumancontribution/validityassessment; publication is not assured.

Priority remaining: (1)restore exact canonicalsource/raw/reference provenance andrerunfullintegration,(2)independentstatisticalgeneticsreview includingcovariancediagnostics,(3)completeauthor/ethics/funding/COI/licensingfields,(4)authorizedcompetitorsupplement/build/allelecomparison,(5)onlythenconsider independenttwo-traitsourceacquisition and newprospective nativeprotocol. No contact,submission ormainmergeoccurred.

Start with SCIENTIFIC_STATUS.md, then PUBLICATION_READINESS.md and MANUSCRIPT/brain6_manuscript.pdf. Numerical provenance and execution logs are in qc/ and review/. REPRODUCE.md distinguishes arithmetic/software reproducibility from blocked native scientific reproduction. The final continuation receipt records commits, CI state, output hashes and immutable-file comparison.
''')
# Readability: use normal prose in final handoff rather than compact internal notation.
p=AREA/'FINAL_HANDOFF.md';s=p.read_text();replacements={'through7October2026':'through 7 October 2026','branchbrain6/':'branch brain6/','failed17,465-cellLAVA':'failed 17,465-cell LAVA','complete8,465secondaryfamily':'complete 8,465 secondary family','same-paircandidate':'same-pair candidate','single-signalABF':'single-signal ABF','No nativeGWAS/LD':'No native GWAS/LD','orFDR-significantpositionalenrichment':'or FDR-significant positional enrichment'}
for a,b in replacements.items():s=s.replace(a,b)
p.write_text(s)
print('Reports, source/claim ledgers and submission materials generated; final compiled manuscript/figures follow.')
# Apply final readable editorial text and completed public build checks.
import subprocess,sys
subprocess.run([sys.executable,str(AREA/'scripts/polish_handoff.py')],check=True)
