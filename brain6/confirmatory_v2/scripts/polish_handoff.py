#!/usr/bin/env python3
"""Readable final editorial material; no modification of scientific results."""
from pathlib import Path
from audit_evidence import AREA
M=AREA/'MANUSCRIPT'
def put(p,s):p.write_text(s.strip()+'\n')
put(AREA/'FINAL_HANDOFF.md','''# Final Brain6 handoff

**Track A: NO-GO. Track B: complete submission-preparation package for human scientific review; not ready for submission approval.**

Completed work includes the current repository, branch and local archive inventory; historical table, hash and numerical audits; source/access and claim ledgers; date-bounded literature searches and targeted primary paper/supplement checks through 7 October 2026; a candidate novelty crosswalk; explicit integration/source-free CI partitions; and a complete manuscript, vector figures, numerical sources, supplement and reporting/submission materials. The new branch is `brain6/continuation-2026-10-07`. Main and historical scientific outputs are preserved.

Independent checks reproduce the 72/35 global counts, 25 candidates in 20 regions, failed 17,465-cell LAVA decision, complete 8,465-slot secondary family, zero same-pair candidate overlaps, six ABF prior reweightings, 54-tissue/1,680-pathway arithmetic and 77 available hashes. One original canonical source is missing, and the original auditor correctly fails. No native GWAS/LD replay or independent biological replication is claimed.

New audit findings: the significant chr11 ADHD block has derived correlation 1.1037, requiring covariance-specific interpretation. Of 2,199 numeric correlations, 723 are outside bounds; 1,105 are missing. Zu 2026 reports the chr5 rs77960 shared region with PP4=0.99, defeating a first-region-discovery claim. Five of six ADHD intervals contain coordinate-consistent prior leads, and two exact Brain6 lead rsIDs recur. All 14 published lead positions match Ensembl GRCh37. A prior insomnia–depression preprint also has a closely overlapping chr11 regional signal, although full source context/build equivalence remains unresolved. These are validation and novelty audit findings, not new biological discoveries.

LAVA remains `FAILED_QC_NOT_PROMOTED`; SuSiE remains held by invalid LD and absent reference provenance. Single-signal ABF remains conditional. No effector, independent two-trait locus replication or FDR-significant positional enrichment was established. Unestimated and inapplicable outcomes are not biological nulls.

Supported wording describes an inherited global map, an explicit screening universe, a failed primary gate, three corrected secondary covariance blocks without same-pair candidate validation, prior-sensitive exploratory colocalization and bounded negative/subset results. No new shared causal variant, mechanism, cross-disease causal effect or general failure of these methods is supported.

Genome Medicine and Molecular Psychiatry are NO-GO for a strong original discovery submission. The Genome Medicine mental-health collection is topically relevant but does not change that verdict. BMC Medical Genomics is a conditional specialist fallback for the bounded paper, after human assessment of its scientific contribution. Publication is not assured.

Remaining priorities are: restore exact canonical source/raw/reference provenance and rerun full integration; obtain independent statistical-genetics review of covariance diagnostics and cohort/model assumptions; complete author/ethics/funding/COI/licensing fields; obtain authorized competing supplements and verify build/allele/LD equivalence; then consider independently sourced two-trait data under a new prospective native protocol. No third party was contacted, manuscript submitted or change merged to main.

Start with [Scientific status](SCIENTIFIC_STATUS.md), [Publication readiness](PUBLICATION_READINESS.md) and [Manuscript PDF](MANUSCRIPT/brain6_manuscript.pdf). The numerical provenance and execution logs are in qc/ and review/. REPRODUCE.md separates arithmetic/software validation from blocked native scientific reproduction. qc/final_receipt.json records final commits, CI state, output hashes and immutable-file comparisons.
''')
put(M/'SUPPLEMENTARY_METHODS_AND_RESULTS.md','''# Supplementary methods and results

This continuation audits a secondary table/provenance package. Historical native methods are described from their frozen protocols and receipts; new numerical checks replay arithmetic, not raw GWAS or native LD. Full native-input reproduction remains incomplete.

S1. The global source table contains all 72 estimates, standard errors, original q values, intercepts, sample sizes and source definitions. The original correction family is 396. Psychiatric and neurodegenerative detection counts are not a powered architecture comparison.

S2. Candidate and region tables preserve 25 pair identities and 20 display merges. Five tracks are insomnia–ADHD, insomnia–MDD, long sleep–bipolar, long sleep–PD and long sleep–SCZ. The PLACO_PLUS lineage includes protected legacy insomnia–ADHD. Lead P<1e-8 is the five-track rule (5e-8/5). Greedy clumping used LAVA UKB European v1.1 reference (100,000 samples), r²=0.1 and 500-kb windows with exact rsID/CHR/BP matching. Lead flanks ±500 kb merge within pair. The locus rule was frozen after initial screen outcomes and before annotation, not fully before discovery.

S3. Canonical LAVA includes 17,465 trait–locus evaluations across seven traits and 2,495 GRCh37 loci. The local-h² gate and complete family QC precede bivariate inference. NOT_RUN reasons are low local h² (3,564), fewer shared reference variants than minimum K (154) and other minimum-K failures (2). By trait: long sleep 1,291; insomnia 671; PD 636; MDD 589; ADHD 223; bipolar 168; SCZ 142. The native-N insomnia pilot has 47 tested/41 unestimated versus 44/44, a 3.41-percentage-point gain and 6.98% low-h² reduction, below alternatives of ten points or 25%. No full rescue was admitted. Hypothetical recovery bounds are not achieved results.

S4. Secondary SUPERGNOVA contains all 8,465 slots. The historical reference is 503-person 1000 Genomes EUR GRCh37, common unique nonpalindromic rsIDs with a genetic map. A recorded failed combined-map cM-reset attempt preceded per-chromosome reconstruction. Minimum shared SNPs is 120. Signed Z/study-total N were admitted only for two insomnia pairs as limited meta-analysis approximations. Three long-sleep pairs were inapplicable before outcomes. MHC is excluded; 1,693 blocks per pair and the complete five-pair Bonferroni denominator remain. Historical four-worker eligible-pair runs, implementation patches and source hashes are in run_provenance.json. The new Wald replay is P=2Φ(-|rho|/sqrt(var)); the retrospective correlation diagnostic does not change eligibility or the testing family.

S5. Fine-mapping has 50 candidate–trait units: 28 ABF estimates, 24 one-signal 95% credible sets, four insufficient posterior associations, nine input-QC holds and 13 long-sleep model holds. Conditional set sizes range 3–1,445 (median 53). These are historical model-conditional posteriors, not new SuSiE PIPs. All 32 assessed LD matrices failed. Six trait–trait coloc estimates retain full H0–H4 vectors, while 19 other pairs are unestimated. Historical coloc 5.2.3 priors are p1=p2=1e-4 and p12=1e-5 with 1e-6 and 5e-5 sensitivities. Analytic reweighting multiplies H4 by the prior ratio and renormalizes all five hypotheses; it checks serialized arithmetic, not raw Bayes factors.

S6. Molecular contexts are cortex/frontal cortex eQTL and LeafCutter sQTL in eQTL Catalogue release 7. GRCh37 windows were lifted to GRCh38; unique position/unordered alleles were matched. Admission required at least 500 shared variants and QTL plus component-GWAS P≤5e-8. QTL N=AN/2 and estimated sdY are model limitations. Four actual ADHD-component tests have maximum H4=0.0824. No pair-level effector, cell-specific mechanism or experimental validation is established. Published QTL credible sets and regulatory overlaps remain coordinate context.

S7. The full 20-region enrichment failed its matched-control requirement in two regions (four and seven anchors versus ten). A separately frozen 18-region/21-candidate sensitivity generated 10,000 accepted nonoverlapping matched null sets. All 54 GTEx v8 bulk tissues and 1,680 Reactome pathways were tested. The new check recomputed P=(extreme+1)/10001 and BH in separate complete families; it did not rerun permutations. Minimum q is 0.113388661134 for tissue and 0.2426424024267 for pathways. No test survives q<0.05. There is no full-set P value, causal-weighted enrichment or cell-resolved result.

S8. FinnGen ADHD supplies global directional context (rg=0.3817, SE=0.0521, P=2.327e-13), with original insomnia reused and complete participant independence unverified. The PGC MDD2025 no-UKB sensitivity has rg=0.4771, SE=0.0228, P=4.65e-97, with nine named PGC cohort overlaps. These historical results were not rerun and are not two-trait locus replications. Adjacent long-sleep lookups change thresholds and overlap UKB. No MR was run in this continuation.

S9. The novelty audit retains date-bounded search queries, truncation, access failures and review depth. Zu Table 3, actual supplement and all 14 insomnia–ADHD GRCh37 lead mappings were inspected. Catalog v2 release 2026-10-03 uses GRCh38.p14. Competitor PMID queries with zero indexed records and unmatched loci do not prove novelty. Full Schipper/Xue/Lin context remains incomplete where access failed.

S10. The final package has 77 available hash matches, one missing original canonical source and no mismatches. Native LD/GWAS/reference inputs were not replayed. Source-free tests exercise synthetic fixtures and numerical functions only; archived-source and native tests remain separate requirements. Exact commands, environments, residuals, claim wording and hash manifests are delivered. No hidden source replacement, threshold relaxation or synthetic empirical output was admitted.
''')
put(M/'JOURNAL_SUBMISSION_CHECKLIST.md','''# Journal-specific submission preparation checklist

Genome Medicine: structured abstract ≤350 words; three to ten keywords; title, authors, affiliations and correspondence; Background, Methods, Results, Discussion and Conclusions; abbreviations; full ethics, consent, data, COI, funding, contributions and acknowledgments declarations. The source follows this structure. Collection deadline is 12 May 2027 as checked 7 October 2026; topical fit does not establish editorial readiness.

Molecular Psychiatry: use MP_ABSTRACT.txt (150–250 words). Follow conservative limits of 3,500 main words, five display items and 75 references until investigators resolve the live guideline/linked PDF conflict. Four figures plus one table are supplied; vector figures are separate and legends follow references. Authors must confirm current limits on submission day. No author/ethics/COI approval is auto-filled.

BMC Medical Genomics fallback: medically relevant genomic scope is plausible. Confirm current article requirements and scientific merit of the bounded secondary-validation contribution. No automated submission is authorized.

Before any submission: resolve original-source audit hold; obtain independent statistical review; verify all numerical claims, references, source permissions and permanent code/data archive; complete declarations; obtain all author approvals; keep all negative and unestimated results visible; remove internal review status only after approval. Track A remains NO-GO.
''')
# Apply the completed independent build check without claiming paper allele identity.
p=AREA/'NOVELTY_AUDIT.md';s=p.read_text().replace('Ensembl GRCh37 maps rs77960 to the published coordinate','Ensembl GRCh37 maps all 14 published insomnia–ADHD leads to the published coordinates (14/14 exact checks), including rs77960').replace('GRCh37 consistency was independently checked.','GRCh37 consistency was independently checked for all 14 leads; effect-allele alignment was not verified.');p.write_text(s)
