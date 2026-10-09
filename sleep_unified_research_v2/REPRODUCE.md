# Addendum reproduction and integrity

Run from repository root. Completed receipts and destinations are immutable; use a newly versioned experiment directory for another execution. Do not rerun archival recovery into the sealed v1 package.

The actual v2 sequence was: one read-only sequential 56.99 GB archive scan; metadata filename walk over readable SSD namespaces; HEAD checks; a preserved slow sequential transfer failure; a verified 1 MiB byte-range probe; separately planned exact ranged acquisition; existing-file hashes; exact DIRECT_TSV decompression; and one original native prefilter at a time. The inaccessible OS `.TemporaryItems` namespace and excluded directories are explicit coverage qualifications. No filename match alone was admitted.

Scripts `11`–`17` record their plans and comparisons. The ranged downloader caps eight workers and 64 MiB response bodies, requires exact HTTP 206/Content-Range/ETag/length, rechecks every chunk, and admits only exact historical complete-file SHA-256. Chunk bodies, completed bodies and the old failed partial remain on the new SSD namespace. The frozen plan budgets their entire 13,916,278,250-byte footprint; additional direct-source/prefilter plans have separate bounds. Network acquisition does not alter allele coding or source versions.

The native prefilter wrapper invokes unchanged `scripts/21_prefilter_hm3.py`, whose Git blob matches original main `659d01cf`, under the recovered Python 3.9.23 environment. It binds before/after source, allowlist, native script, Python binary, QC and acquisition receipt hashes. Its own version, command and resource limits are frozen in each plan. It compares total/retained rows and deterministic gzip bytes to historical QC records. Historical and newly calculated gzip hashes agree; this does not complete downstream allele/MAF/INFO/MHC/sample-size harmonization, munging or LDSC.

Portable validation needs only the two committed evidence packages:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s sleep_unified_research_v1/tests -v
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s sleep_unified_research_v2/tests -v
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/10_seal_package.py --verify
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v2/scripts/18_seal_addendum.py --verify
```

Clean-checkout validation is recorded separately from scientific validation. Portable tests inspect data/proof invariants; the independent reviewer separately hashes actual SSD files. Neither check proves full native analysis, achieved power, calibrated heterogeneity or biological novelty.

Full historical replay remains governed by the v1 native runner and 3 GiB internal / 5 GiB SSD guards. Its 190 planned jobs reproduce estimators from processed inputs; complete source harmonization is additional work. Because v1 is sealed, future native outputs must be redirected in a new reviewed execution version rather than written into v1. The clinical-source v3 amendment authorizes no new pair outcome until all remaining source, h2 and power gates pass.
