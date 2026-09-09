# Discovery extension adversarial review

Generated: 2026-08-29T13:12:15Z

Verdict: **QUALIFIED PASS — GLOBAL DISCOVERY AND INDEPENDENT REPLICATION; DOWNSTREAM CLAIMS WITHHELD**

## Scope and immutable boundary

The core verifier reports: `CORE_CHECKPOINT_OK commit=22df92c8d6894f93c50561ae8a86acf2f48f50cb traits=45 pairs=396 result_sha256=161756ac61775ad3393572cf5fbdf51690eb47ff010b8a50dbcc5fabf79bc077`

The prospectively locked extension contains 100 traits and 1200 planned tests. It remains analytically and multiplicity-separated from the 45-trait/396-pair core. This audit accepts only the global discovery, literature, and independent-replication claims supported by real canonical artifacts; downstream claims are withheld.

Present real finding classes: extension h2, global rg, pair novelty, replication.

Absent real finding classes: local architecture, pleiotropic loci, fine-mapping/colocalization, mechanistic synthesis.

## Findings-level stress test

| Challenge | Evidence | Adversarial disposition |
|---|---|---|
| Core drift | 45 traits, 396 pairs, checkpoint verifier PASS | PASS; no core file is used as an extension output. |
| Post-result selection | 100-trait order/hash predates all 1,200 rg results | PASS; no trait was added or removed after result access. |
| Weak h2 / intercept inflation | 100/100 primary-pass; min h2 Z=4.8014; max intercept=1.1663 | PASS under the frozen Z>=4 and intercept<=1.2 gates. |
| Multiple testing | One isolated 1,200-test BH family; 603 FDR<0.05 | PASS; core P values were never pooled. |
| Discovery sample overlap | max absolute cross-trait intercept=0.1386 | DISCLOSE; UKB sample overlap remains plausible and discovery is not independent confirmation. |
| SNP overlap / ancestry | min valid overlap=1,128,043; discovery ancestry EUR | PASS for LDSC comparability, not ancestry generalizability. |
| Literature duplication | 603/603 hits audited; classes={'KNOWN_BUT_NEW_DATASET': 1, 'NO_DIRECT_RG_FOUND': 308, 'PARTIAL_EXTENSION': 294} | QUALIFIED PASS; one exact same-pair/new-dataset result retained, zero APPARENTLY_NOVEL, no first-ever claim. |
| Replication attrition | classes={'DIRECTIONALLY_CONCORDANT': 18, 'NO_INDEPENDENT_DATASET': 159, 'REPLICATED': 23, 'UNDERPOWERED': 17} | DISCLOSE; 159 unavailable and 17 underpowered cannot be interpreted as failures or successes. |
| Replication independence | 41 tested in FinnGen R13 or MVP EUR; non-overlap required in the locked manifest | PASS for cohort independence; discovery and replication phenotype definitions still require row-level qualification. |
| Direction stability | discordant=0/41 tested | PASS; absence of discordance does not rescue non-significant tests. |
| Phenotype-definition mismatch | exact=18, comparable-with-differences=5 among 23 replicated | QUALIFIED; five findings cannot be described as exact phenotype replication. |
| Effect heterogeneity | 7/23 replicated have heterogeneity P<0.05 | QUALIFIED; both estimates are retained and pooled-effect language is avoided. |
| Replication h2/intercepts | 4/13 sources failed prespecified h2/intercept gates | PASS fail-closed behavior; affected 17 pairs remain UNDERPOWERED. |
| Local/pleiotropy substitution | 23 replicated PLACO+ candidates locked; no real result artifact | PASS_BLOCKED; HapMap3 LDSC inputs are not relabeled as dense genome-wide inputs. |
| Fine-mapping/colocalization | signed LD, dense locus statistics, and QTL family absent | PASS_BLOCKED; physical overlap or global rg is not called colocalization. |
| Causal-language overreach | report claim limit forbids first-ever, causal, local, pleiotropic, colocalized, and mechanistic claims | PASS; strongest claim is independently replicated global genetic correlation. |

## Challenge ledger

| Risk | Required control | Current status | Evidence |
|---|---|---|---|
| Core drift | Re-run byte-level core checkpoint before and after every extension stage. | PASS | CORE_CHECKPOINT_OK commit=22df92c8d6894f93c50561ae8a86acf2f48f50cb traits=45 pairs=396 result_sha256=161756ac61775ad3393572cf5fbdf51690eb47ff010b8a50dbcc5fabf79bc077 |
| Post-result panel selection | Panel membership/order/hash must predate extension rg. | PASS | 2026-08-28T04:13:05Z |
| Multiple-testing leakage | Never combine the extension BH family with the immutable 396 core tests. | PASS_CONTRACT | extension_harmonization_policy.json |
| Weak h2 | Exclude rerun h2 Z<4 from primary rg without replacement. | PASS_REAL_H2_GATE_APPLIED | present:sha256=95b22a7da9924945643f43660080b9773496421d4ca7c184eac03a03961a626c |
| LDSC intercept inflation | Exclude rerun intercept>1.2 from primary rg and retain sensitivity status. | PASS_REAL_H2_GATE_APPLIED | present:sha256=95b22a7da9924945643f43660080b9773496421d4ca7c184eac03a03961a626c |
| UK Biobank sample overlap | Inspect cross-trait intercepts and disclose overlap; do not equate LDSC adjustment with independent replication. | PASS_DIAGNOSTICS_RECORDED_DISCLOSURE_REQUIRED | present:sha256=5089bfff7b607f356104f90aba017ff35e1841f66056242b765f05c6e7ba1166 |
| Sparse or proxy phenotypes | Retain exact phenotype definitions and distinguish medication/proxy traits from diagnoses. | PASS_METADATA | candidate_traits.tsv |
| Panel-level novelty inflation | Require pair-level direct/same-phenotype/same-direction/same-sleep-context audit. | PASS_PAIR_AUDIT_COMPLETE | present:sha256=af910799e1ff94de584239879d63495dda537632efd2f4877f1a238e1f2933d8 |
| Replication non-independence | Require non-overlapping participants and separately sourced summary statistics. | PASS_INDEPENDENT_REPLICATION_COMPLETE | present:sha256=c4b4da28640101102cc2784eeb772e8af0ab59f6fb3b3e3247015b9ae496cf84;present:sha256=d613c6ab4ebb33814621a6cb2e2d4e88b680153e32d95683695435ee442d2711;present:sha256=b91cd751d42431e94b2734782357f2f8e3d579344aa202522c8241a212e46dff |
| Global-to-local overreach | Do not call global rg evidence of a shared locus; retain a prespecified globally-null secondary local set. | PASS_QUEUE_LOCKED_INPUTS_BLOCKED | absent;present:sha256=55992eee1891497c42b6a381f9fee2d1ff3979bbfb6cdc294f02fe49951bbbc4;present:sha256=1fe6d3bd40a21423c93e9941e2de644e500e5cdf3fe902f49adcb999715715c7 |
| Pleiotropy or mediated effects | Evaluate horizontal, vertical/mediated, shared-factor, and sample-overlap alternatives; run PLACO+ on genome-wide data only. | PASS_CANDIDATE_FAMILY_LOCKED_INPUTS_BLOCKED | absent;present:sha256=f3508e6ed0ec22ed91a395ec2fe75e5cf7aec303f64674582c9ef83adde5d8fa;present:sha256=dfc107b45ecd92e166c8e4896eac554ec2ec397ae3d79676548e867d1d241069;present:sha256=b28bbc80736f4dab94e7681bceee2f6c5cd57d4f03edc6ca6b5bcfebcf19748c |
| Colocalization overclaim | Report H0-H4, priors, sensitivity, and claim guards; colocalization is not causality. | CODE_READY_INPUTS_BLOCKED | absent;present:sha256=10a13727ffbe9c7f9ca07a2d1a53fe2489ae17d36f86a6354002d58cbfad4328 |
| Synthetic/real result contamination | Synthetic tests stay under synthetic or temporary paths and carry explicit markers. | PASS | seven isolated synthetic workflows |
| Storage-driven partial acquisition | Do not silently analyze a result-selected subset of the locked panel. | PASS_LOCKED_PANEL_COMPLETE | streaming_receipts=100/100;mirror=BLOCKED_INSUFFICIENT_STORAGE |
| Local LD-reference mismatch | Require checksum-locked ancestry-matched LAVA/HDL-L references; do not fall back silently to a smaller panel. | PASS_BLOCKED | present:sha256=55992eee1891497c42b6a381f9fee2d1ff3979bbfb6cdc294f02fe49951bbbc4 |
| Local multiplicity or h2-gate leakage | Freeze the pair-by-locus family, LAVA local-h2 Bonferroni gate, and local-rg BH family before result access. | PASS_CONTRACT | config/local_architecture_contract.json |
| Replication candidate attrition | Lock the entire Tier A/B candidate family before source curation and preserve NO_INDEPENDENT_DATASET outcomes. | PASS_CONTRACT | scripts/21_prepare_replication_queue.py+22_lock_replication_manifest.py+23_collate_replication.py |
| LAVA simulation instability | Freeze a deterministic nonzero simulation seed per selected pair and retain it in the manifest. | PASS_CONTRACT | config/local_architecture_contract.json |
| Fine-mapping LD or allele mismatch | Require identical dense SNP order/alleles, checksum-locked signed LD, PSD/symmetry checks, RSS-LD s, and kriging outlier diagnostics. | PASS_CONTRACT_INPUTS_BLOCKED | config/fine_mapping_colocalization_contract.json |
| Molecular-QTL source attrition | Complete eQTL/sQTL/pQTL searches per locked locus and preserve evidence-backed NO_SUITABLE_DATASET rows. | PASS_CONTRACT_INPUTS_BLOCKED | scripts/33_prepare_finemapping_queue.py+34_lock_finemapping_manifest.py |
| Coloc prior cherry-picking | Lock p1/p2/p12 and the full p12 sensitivity grid before result access; label prior robustness separately. | PASS_CONTRACT_INPUTS_BLOCKED | config/fine_mapping_colocalization_contract.json |
| Mechanistic source or release drift | Require exact release/accession, local query/result snapshot, SHA-256, access date, license/terms, and primary citation for every supported evidence row. | PASS_CONTRACT_INPUTS_BLOCKED | config/mechanistic_annotation_contract.json+config/mechanistic_sources.tsv |
| Mechanistic chain gap concealment | Emit every required chain edge, mark unsupported edges MISSING, and forbid a complete causal narrative from annotation-only evidence. | PASS_CONTRACT_INPUTS_BLOCKED | scripts/39_validate_mechanistic_evidence.py+40_synthesize_mechanisms.py |

## Major residual risks

- The discovery screen is UKB-centered and restricted to EUR summary statistics; transferability is untested.
- 603/1,200 FDR hits reflect pervasive polygenic correlation and shared-cohort structure as well as biology; the LDSC intercept is a diagnostic, not proof that overlap is harmless.
- Replication availability is highly selective: only 58/217 candidates had a suitable source, and only 41 survived source h2/intercept gates.
- Five replicated comparisons use comparable rather than identical phenotype definitions, and seven show nominal effect heterogeneity.
- The literature audit found no pair that met the stringent APPARENTLY_NOVEL bar; 'underreported' is the maximum defensible novelty language.
- Global rg cannot identify loci, distinguish horizontal from vertical pleiotropy, nominate genes, or establish causality.

## Final adversarial conclusion

The real global extension and independent replication family pass their locked computational, multiplicity, ancestry, overlap-documentation, and provenance checks. Twenty-three pairs support a qualified replicated-global-rg claim. The audit rejects stronger novelty language and rejects all local, pleiotropic, fine-mapped, colocalized, gene, mechanism, or causal interpretations because the required real artifacts do not exist.

This is a qualified findings-level pass, not a declaration that every planned downstream stage succeeded.
