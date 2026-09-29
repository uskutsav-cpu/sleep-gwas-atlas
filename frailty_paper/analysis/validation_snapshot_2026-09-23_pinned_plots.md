# Selected plot regeneration and regression validation

**Date:** 2026-09-23

The FI global-rg forest, sleep × latent-frailty heatmap, and insomnia/sleep-apnea frailty-dimension forest were regenerated with Python 3.11.11 and Matplotlib 3.9.4, matching the repository pins. Their provenance sidecars record the source and output SHA-256 values. PNG dimensions are 3174×2294, 2963×2194, and 4223×2376, respectively; each PDF contains one page. The figures remain provisional and do not establish replication or support downstream locus/mechanism claims.

The package suite passes 64/64 tests in the existing Python 3.13 environment. Five focused plot/path regression tests pass in the isolated Python 3.11.11 environment. The complete test suite in that isolated environment remains incomplete because `lxml` is unavailable in the offline package cache (`test_pubmed_acquisition` cannot import). The existing project `.venv` was not replaced. Figure 1–2 runtime evidence remains as previously recorded.
