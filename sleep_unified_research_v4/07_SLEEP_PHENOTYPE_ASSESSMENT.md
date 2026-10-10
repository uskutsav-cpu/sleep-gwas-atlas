# Sleep-phenotype evidence assessment

Status: COMPLETE_WITH_QUALIFICATIONS for the documented phenotype crosswalk; BLOCKED_HUMAN_REVIEW for clinical interpretation; NOT_ADMISSIBLE for a new measurement-difference claim. This is a technical evidence inventory, not manuscript text.

The 12-row `tables/sleep_clinician_core_crosswalk_v4.tsv` records exact source, construct, questionnaire/device coding, window, effect units and direction, N/case/control counts, ancestry, build, reliability limits, overlap and outstanding interpretation question. Historical rounded heritability fields are explicitly identified as historical. Newly validated native estimates reside in `tables/core_native_full_precision_h2_v4.tsv`; the crosswalk does not replace those estimates or certify comparability across scales.

| Measurements | Required interpretation boundary |
|---|---|
| Insomnia | Frequent initiation/maintenance complaints; chronic insomnia disorder criteria were not required. |
| Reported duration and selected short/long categories | Twenty-four-hour recalled duration includes naps; both tails share normal-duration controls. No diagnosis, independent mechanism or causal U-shape follows. |
| Chronotype | Morning-person preference, with direction opposite to later device midpoint; not measured biological circadian phase. |
| Sleepiness and napping | Ordinal propensity to doze and nap frequency are distinct constructs. Neither is a physiological subtype or an objective sleepiness test. |
| Snoring | Participant-reported complaint by a close contact; no respiratory-event or apnea diagnosis is established. |
| FinnGen R9 sleep apnea | Registry endpoint with healthcare ascertainment; an ICD9 3472/3472A labeling anomaly needs Finnish coding expertise. The reported code-related count is not a count of uniquely misclassified cases. |
| Device efficiency, duration and midpoint | Main-period activity-derived measurements in a selected later UKB subset; quiet wake, shared boundaries, transforms and assessment lag limit comparison with baseline questionnaires. |

The two-row new-source crosswalk and five-row objective-source screen preserve clinical definitions and availability limits. Public MVP GIA Phe_327_4 rules were resolved in v2 and remain qualified: at least two mapped ICD instances versus zero for controls. Spacing, chronicity, observation horizon and accession-specific effect/CI semantics are still unresolved. FinnGen R13 F5_INSOMNIA is a related registry construct, with missing per-variant INFO/N and no clinical signoff; acquisition does not establish independent replication.

No direct difference in genetic correlations across measurement modes is admitted. Such a test requires a prespecified disease target, calibrated joint uncertainty including heritability denominators, construct and overlap qualification, prospective precision and a valid independent source. Counts of significant rows do not answer that question. No exact-source reliability coefficient or unestimated genetic absence is invented.

Evidence: `reviews/sleep_clinician_review_v4.md`, its immutable receipt and SHA256 inventory; `tables/sleep_clinician_core_crosswalk_v4.tsv`, `tables/sleep_clinician_new_sleep_crosswalk_v4.tsv`, `tables/sleep_clinician_objective_source_screen_v4.tsv`, and `tables/sleep_clinician_human_source_questions_v4.tsv`. Outstanding qualified-human questions are retained verbatim in those tables.
