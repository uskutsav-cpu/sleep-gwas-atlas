# Independent integrated scientific review of the V8 interpretation candidate

Reviewed 11 October 2026. Finite scope: actual V8 candidate rows/schema/hold/manifest, four direct scientific lookups, the rendered two new PNG figures and their embedded labels, and the biological-follow-up decision. No native fit, historical rerun, broad unchanged-source audit or download was performed. Separate manuscript captions were not yet present at review time; the embedded figure captions were inspected visually.

## Decision

**Scientifically defensible as a qualified interpretation and manuscript-preparation candidate. It is not a cleared public deposition or a claim of journal readiness.** The release correctly separates historical reproduction, a current source-direction hold, existing QC exclusions and the one new exploratory transport estimate. No new independent two-trait replication, source contrast, causal or molecular result should be inferred.

Two documentation items should accompany manuscript/release use: state that the core heatmap includes all QC-eligible point estimates irrespective of BH significance, and supply a package-level link to the new FinnGen/MVP source contract and complete limitations. These require bounded documentation, not new fits or infrastructure. The frozen seven legacy tables need not be rewritten to achieve that.

## Actual quarantine and testing-family treatment

`release_candidate/core_rg.tsv` has 396 records, all labelled `family_size=396` and `adjustment_method=frozen_BH_396`. Exactly 36 use `delange_2017_ibd`, `delange_2017_crohn` or `delange_2017_uc`: twelve sleep constructs × three outcome sources. All 36 have empty rg, z and both 95% interval fields. Each row explicitly withholds biological interpretation and limits retained SE/P/q/threshold/class fields to historical reproduction metadata. Neither empty values nor retained historical P/q signify biological absence or cleared evidence of sharing.

The uncertainty concerns source direction. Published control variants support a concern, while the exact source-wide mechanism remains unresolved. The release appropriately does not describe an authoritative correction, silently reverse signed estimates, claim every source row is reversed, or remove the records from the testing family. The current hold is a scientific interpretation hold, not a claim that numerical reproduction failed. Conditional sign-reexpression diagnostics, if retained separately, cannot clear that hold without a defensible source-wide contract.

The new heatmap's total is consistent with actual candidate fields: 36 source-held records, 63 additional unheld records failing one or both existing pairwise h²-Z/intercept conditions, and 297 displayed point estimates. H and x are visually distinct from the near-zero color. Removing estimates from biological interpretation does not redefine the historical BH denominator. No FDR recomputation or new testing family was introduced by this view.

The manifest and schema distinguish the one-row exploratory table from the core396, extension1200 and historical validation217 families. The manifest reports the seven other legacy TSV outputs with the same recorded historical digests. This finite review did not rehash or reanalyse those unchanged tables; it inspected the manifest's family-preservation claim and the relevant current rows, without claiming a new independent numerical audit of historical estimates.

## Actual scientific lookup behavior

Four independent direct calls to `resource_tools/query_atlas.py`'s existing reader were examined:

| Query | Actual returned behavior | Scientific assessment |
| --- | --- | --- |
| Core insomnia–IBD | `HISTORICAL_RECORD_WITH_INTERPRETATION_HOLD`; blank rg, z and CI; hold basis and historical scope supplied | Correct abstention from signed biological interpretation |
| Core insomnia–frailty | Estimate 0.6405056712964976 plus source/coding/QC/marginal uncertainty metadata | Valid retrieval example; historical association, not new discovery or a frailty-project modification |
| Validation insomnia–MVP abdominal pain | Estimate 0.40982718740242974; links PanUKB discovery; reports sleep reuse and unavailable calibrated difference test | Correctly limits this to external outcome-side validation |
| Core insomnia–unknown outcome | `NO_ESTIMATE_IN_THIS_RELEASE_TABLE`; explicit no biological-null inference | Correct distinction between missing record and absent association |

The required V8 source hold is read from the candidate directory; the reader checks missing-hold failure for a V8 manifest. Its saved examples disclose an ineligible melanoma record's QC fields; retrieving a historical estimate is not scientific clearance. Users should apply the exposed QC and interpretation fields before prioritizing biology. These examples demonstrate transparent evidence retrieval, not completed human usability validation. The accompanying pilot remains necessary for any empirical usability claim.

## Actual new result and source-transport figure

`exploratory_source_transport.tsv` contains exactly one separately labelled record: FinnGen registry insomnia/MVP abdominal pain, rg=0.39912075217856924, SE=0.050695938374874, CI approximately 0.29975854–0.49848297, nominal P=3.4669408448782324×10⁻¹⁵. Its scope states new sleep source, reused historical MVP outcome, CI95/log-Wald and constant-effective-N assumptions, no independent replication and no cross-source difference test. Values agree with the independently reviewed actual native result; no estimator rerun was needed.

The visually inspected `insomnia_abdominal_pain_source_transport.png` shows three source-defined estimates with marginal 95% intervals:

- UKB complaints/PanUKB pain: 0.559 [0.489, 0.629].
- UKB complaints/MVP pain: 0.410 [0.345, 0.475].
- FinnGen registry insomnia/MVP pain: 0.399 [0.300, 0.498].

These match the current release's historical source rows and reviewed new result at displayed precision. The new row is clearly exploratory, and the visible footer discloses the shared MVP outcome, lack of a calibrated between-source difference test and source/CI/effective-N assumptions. The axis presents correlation rather than a causal effect, and the points are not connected into a biological trajectory. Figure placement beside historical estimates cannot itself create replication or equivalence.

Manuscript caption should retain: historical selection informed the new pair; the two MVP rows share the same outcome GWAS; clinical/source constructs differ with cohort and ascertainment; intervals are marginal, assumption-conditional and not simultaneous; the clinical-source result supports positive transport without identifying a pure measurement effect, causal mechanism or a source difference. The visually close two MVP estimates do not establish statistical equivalence. The larger PanUKB point does not establish stronger sharing after accounting for joint uncertainty.

The one-row TSV's brevity is acceptable if the distributed candidate includes an explicit source-contract supplement. The seven preserved legacy tables still leave MVP sample counts/build unknown in their older crosswalk; that is appropriately historical provenance, but readers must be directed to the newly verified accession/count/build evidence and the remaining assumptions. The supplement should include unknown individual overlap, Finnish LD-reference mismatch, the HM3-without-INFO route and assumed rather than verified SNP-wise N. It should not describe FinnGen diagnostic coding as clinically adjudicated chronic insomnia.

## Actual core figure

The rendered `core_interpretive_quarantine.png` is legible and distinguishes source holds from pairwise QC failures. It states that colors are point estimates, uncertainty is in the release table, H/x are withheld interpretations rather than zero effects, there is no causal/between-trait contrast claim, and the complete original 396-pair family is unchanged. The symmetric color scale labelled −0.7/+0.7 conveys direction and magnitude; it is not a significance or reliability scale.

Caption should explicitly state that all 297 QC-eligible unheld estimates are colored irrespective of BH significance. No significance symbols appear in the figure, so darker color must not be described as statistically stronger evidence without consulting uncertainty and frozen q-values. Unheld status does not establish perfect harmonization, participant disjointness or validated measurement reliability; each displayed estimate retains its source-specific qualifications. This heatmap is a source-defined descriptive map, not a discovery panel.

## Biological use case and claims to retain

The new evidence is positive project-level sleep-side source transport under qualified assumptions. It is useful for prioritizing a sleep–pain genomic question and demonstrating transparent source/uncertainty interpretation. Field-level novelty is unestablished, and the already identified clinical-versus-habitual sleep prior art narrows broad novelty claims.

The independent `biological_followup_admissibility_v8.md` appropriately retains this global use case and declines to manufacture a gene, locus or tissue claim from it. No dense regional variance/LD/molecular contract has been admitted for the selected pair. This is not proof that a mechanism is absent or impossible; it is an evidence-specific limit on current claims. Repeated global fits and broad local/QTL scans would not close that precise scientific gap automatically.

Essential remaining manuscript/release work is bounded: include the figures, qualified new row, source-contract references and selection-aware text; preserve the bowel-source hold everywhere; complete actual human quality/usability assessment before asserting it occurred; and obtain source-specific release clearance before public deposition. Optional contingent work is an independently designed regional or cohort validation once compatible inputs and adequate precision are established. No unchanged-file audit expansion, controller version or historical native rerun is justified by this review.

## Suggested captions

**Core interpretation map.** “Source-defined marginal genetic correlations across the historical 396-pair core. Colors show 297 unheld point estimates passing the existing pairwise h²-Z/intercept conditions, irrespective of BH significance. H marks 36 records with unresolved bowel-source effect direction; x marks 63 additional records failing existing pairwise QC. H/x are withheld interpretations, not zero effects. Marginal uncertainty and original complete-family BH q-values are provided in the candidate release. No causal or between-trait difference inference is made.”

**Sleep–pain source transport.** “Three source-defined insomnia–abdominal-pain estimates and marginal 95% intervals. The FinnGen registry-insomnia/MVP result is a separately specified exploratory analysis selected using historical evidence; the two MVP rows share an outcome GWAS. Cohort and ascertainment differ with the sleep source. The positive new estimate supports qualified common-variant source/construct transport, while effective-N, CI-reconstruction and LD-reference assumptions remain. These intervals do not test between-source differences or establish independent two-trait replication.”

A public deposit remains explicitly uncleared. This review establishes scientific interpretation quality for the current candidate and finite documentation actions; it does not guarantee Cell Genomics acceptance.

## Documentation closure, 11 October 2026

Both bounded documentation actions are complete. Inspected only the new `figures/CAPTIONS.md`, expanded exploratory TSV fields, corresponding schema entry and manifest status. The core caption explicitly states that all 297 QC-eligible estimates are colored regardless of BH significance, preserves the complete 396-pair denominator and distinguishes holds/QC exclusions from biological nulls. The transport caption preserves historical selection, shared MVP outcome, marginal conditional intervals, CI95/log-Wald and effective-N assumptions, HM3 fallback, approximate Finnish/EUR LD and unknown cross-enrollment/relatives. It excludes discovery, independent replication, mechanism, causality and calibrated source contrasts.

The exploratory TSV now has expanded `inference_scope`, `source_contract_ref` and `independent_actual_result_review_ref`; both relative references resolve within the private V8 package and are included in the updated schema. The manifest still marks public deposition uncleared and historical families unchanged. No repeat numerical review, rendered-figure audit, historical-input hashing or estimator fit was performed for this closure. **The two documentation conditions in this finite scientific review are resolved; retain the qualified private manuscript-preparation candidate verdict.** Public code availability is distinct from permission to deposit source data or candidate scientific outputs.
