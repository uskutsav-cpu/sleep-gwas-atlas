# Researcher utility study: integration benefit

Date: 2026-10-11. Status: preparation for independent review; no participant recruited, answer collected, or benefit established. This is a new human-researcher study, distinct from the deferred V7 automated-answer experiment. It cannot be presented as execution of that experiment.

## Scientific question

Does explicit source/phenotype linkage and contextualized result retrieval reduce interpretation errors when researchers receive the same numerical results and source evidence?

The intervention is the existing candidate resource plus its exact-key query interface. The comparator receives the same eight historical-family candidate TSV views, SCHEMA.md, schema.json and source_direction_hold.json, without the query interface or prejoined answer. Both groups can use their normal local text/table tools. Both receive identical external-resource links. No task requires additional primary-paper content: the source dictionary already provides the source-defined information being tested. This comparison evaluates integration plus contextual interpretive guardrails; it does not establish superiority over every public platform, biological novelty, or the benefit of access to additional data.

## Finite pilot design

Invite independent statistical-genetics/sleep researchers only through investigators after the applicable institutional determination and source-use decisions. No invitation is authorized by this protocol. Target 8–12 participants for feasibility, not a powered confirmatory efficacy claim. Exclude atlas builders and anyone who has inspected the scoring key; record relevant experience without identifiable demographic data. Do not replace participants based on scores.

Each participant receives the same 12 task IDs, each once. Assign six tasks to plain-table workflow and six to integrated-query workflow. Use complementary allocations in successive participant pairs, balancing the condition for every task; randomize task order with a frozen seed before any responses. Distinct estimate records avoid direct reuse of numerical answers; general source and interpretation learning across tasks remains a limitation. Record order and prior familiarity. Provide equal orientation time and a common answer template; do not give one group extra source knowledge.

A pilot may stop at the prespecified recruitment window whether scores favor the atlas or not. Investigators must record the recruitment window and intended participant count before collecting answers. No outcome-dependent task/participant substitution is allowed. If recruitment is unavailable, retain the packet as unexecuted; automated agents do not count as human participants.

## Evidence equality and tasks

Both conditions use the exact same frozen V8 interpretation-view files, source dictionary, schema and source-direction hold in shared_evidence_manifest.json. Availability failures affect both conditions and must be retained. No task depends on undisclosed scoring-key information. The 12 tasks cover phenotype meaning, source identification, fixed testing families, marginal versus contrast uncertainty, validation independence, sensitivity interpretation, QC failure, and absent-row interpretation. Tasks are predefined in tasks.json. T04 and T05 use different outcome records to avoid answer transfer for the same estimate. General source/interpretation learning across tasks remains a limitation. They are purposive safety/interpretation tasks, not a representative sample of all genetic research.

Do not run the study on any record with a newly identified unresolved correctness issue. Record any affected task exclusion or replacement before responses, with its reason, and preserve the original task set. Historical estimates remain immutable; any corrected release becomes a separately frozen input version.

## Scoring and analysis

Remove condition labels and participant identity before independent scoring. Two qualified human scorers use the stored source-supported key. Score the required factual units correct/incorrect/not assessable and flag unsupported causality, clinical equivalence, replication or difference claims. Accept semantically equivalent source-supported phrasing; do not require exact strings. Numeric rg, SE and CI answers must agree within 0.0005 absolute (approximately three-decimal reporting), family sizes/threshold classifications exactly. No unstated numerical rounding test is added. A fully supported task requires all required units correct and no material unsupported claim. Scorer disagreements remain visible; conservative agreement is the primary pilot summary, with adjudicated scores separately reported.

Primary descriptive endpoint: fully supported tasks / six assigned tasks in each workflow for each participant. Unanswered, incomplete, withdrawn-task and scorer-unresolved responses count as not fully supported; keep separate reasons. Pre-response instrument failures affecting both workflows are reported separately and excluded only through a documented outcome-blind amendment. Entire participant withdrawals are reported and retained in the flow accounting, with consent governing data retention; do not invent their answers. Report every participant's two workflow totals and task-level completion, errors, order, and missingness. Report paired participant differences and their distribution; do not treat 12 tasks per person as 12 independent people. Task and participant clustering, purposive selection, small pilot sample and learning prohibit a broad confirmatory user-benefit claim. Time-to-answer is a secondary feasibility measure recorded by the participant interface or observer, not model/tool latency. Report incomplete-task censoring explicitly.

No significance threshold or favorable minimum effect is a pilot success gate. Success means the study is executable and reveals measurable failure modes; unfavorable and null results are retained. A later confirmatory design would require a separate prospectively justified sample size and protocol informed by the pilot, with participant/task dependence handled explicitly.

## Required human decisions before execution

Investigators: applicable institutional/consent determination, permitted private sharing of candidate data/documents, independent participant and scorer recruitment, fixed count/window, and any compensation. These are actual prerequisites, not computation tasks. No identities, consent records or contact details belong in Git.

## Deliverables

Task packet and shared evidence inventory; complementary randomized allocations; deidentified responses/times; blinded human scores/disagreement; complete pilot summaries; independent scientific review. Until responses exist, there are no utility results.

## Input freeze amendment before responses

The shared-evidence manifest now freezes V8_INTERPRETIVE_QUARANTINE_CANDIDATE. No response has been collected. Its eight historical-family TSV views preserve unaffected task records; 36 bowel-source core records have signed fields withheld and mandatory interpretation holds. None of the 12 tasks requests an affected bowel-source estimate. Both workflows receive identical hold/schema information. The new exploratory FinnGen estimate is outside the task packet. The original V4 manifest is preserved in input_manifests/preserved_V4_shared_evidence_manifest.json. Input freeze does not satisfy human execution prerequisites.
