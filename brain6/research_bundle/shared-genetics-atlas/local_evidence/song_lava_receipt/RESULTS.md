# LAVA receipt overlapping the published 3p21.31 interval

Of 2,495 saved insomnia–FI LAVA receipts, one has numerical coordinates overlapping the [Song et al. ρ-HESS interval](../LOCAL_FRAILTY_AUDIT.md): LAVA locus **464**, chromosome 3:47,588,462–50,387,742, contains the reported 47,727,212–49,316,972 interval. Coordinate-build equivalence and the published partition settings remain unverified; the two different windows and methods are not interchangeable tests.

The LAVA receipt reports `process_status=PROCESSED` but `univ_status=PARTIAL`: insomnia was tested, **frailty was dropped**, and the bivariate test was `NOT_ELIGIBLE`. The accompanying log specifically reports a negative variance estimate for frailty. Consequently the saved local correlation and its p-value are absent. A missing bivariate estimate is **not a negative or null finding**, and the receipt does not contradict or independently validate the published ρ-HESS result. The log does not establish the cause of the negative variance estimate or a correction for it.

This locus-level inspection is compatible with, but does not resolve, the failed aggregate insomnia LAVA quality gate. It is a post-hoc coordinate lookup, not a rerun. The saved insomnia-only univariate result is not promoted as a validated local finding because the family gate failed and relevant source/reference and calibration issues remain. Coordinate overlap does not establish identical SNP support, sample overlap treatment or statistical estimand.

[Audit code](audit.py) reads all 2,495 receipts and retains the matching receipt, its SHA-256 and the matching log's SHA-256 in [receipt.json](receipt.json). Source files remain on the SSD; their original process output is not modified.

A [cross-trait receipt comparison](cross_trait/RESULTS.md) finds the same FI input in all 12 runs at this locus. Eleven dropped FI; the sleep-apnea run retained it with different marker and component counts but did not qualify for a bivariate test. This is an observation about the archived analysis, not a verified explanation of the negative estimate.
