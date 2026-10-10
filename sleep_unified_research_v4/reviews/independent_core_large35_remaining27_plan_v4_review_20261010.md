# Independent review: remaining27 core replay v4

Read-only structural review of scripts 87/88 and the workspace draft plan. No workers or scientific analyses were run.

SHA-256 identities reviewed:

- `87_prepare_core_large35_remaining27_v4.py`: `366c8b7d080ddbea8ab122a689fdb675df892419c70f028768de8d63e57125ff`
- `88_run_core_large35_remaining27_v4.py`: `17a696ba13d7c2444c3e35f015fb80dfe38995fc3d9299768f03357d69e59976`
- `work/core_large35_remaining27_plan_v4.review.json`: `abe31a477fdfc4c1c28d0ef6253dd62b74d15119d0add3407b5ff30bc88c05da`

The draft contains exactly the ordered v3 members 9–35 (27 traits), with all successor-owned harmonized, munged, receipt, and spool routes in the fresh `large35_bounded_replay_v4_remaining27` SSD namespace. The seven completed v3 members remain bound by their retained output and worker hashes. The snoring evidence now includes hashes for its QC, source gate, harmonized and munged outputs, failed comparison, successful and failed worker receipts, logs, and journals; the comparison and sole step-14 label difference are checked against the exact 0-dropped / 7,168,629-retained values. The runner rechecks this excluded evidence before work and in final identity checks. The original v3 failed attempt remains separate.

Controller accounting is internally consistent: 27 ordered members, five monitored commands per member (135 total), one worker, the original shared mutex, and the inherited 3/5 GiB free-space floors, 2 GiB RSS ceiling, 16 GiB core namespace cap, 300 GiB campaign cap, and 2-hour/96-hour deadlines. Estimator calls remain zero. The v3 validator and cleanup helper are reused unchanged; the cleanup helper’s plan-parent requirement is met when the frozen production plan is emitted inside the v4 namespace. The workspace draft is for review and must not be passed as the production plan path.

I found no concrete correctness, data-integrity, or safety defect in the reviewed candidate. This review supports separate root execution admission; it is not evidence that the 27 production replays have run or passed.
