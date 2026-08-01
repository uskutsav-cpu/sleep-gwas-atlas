# Environment lock — Phase 2/3 toolchain

Recorded 2026-08-01T21:33:35Z on Darwin arm64, 8 GB RAM.

## Installed and verified

| component | version / commit | how |
|---|---|---|
| R | 4.6.1 (2026-06-24) | `brew install r` |
| Homebrew | 6.0.14 | preinstalled |
| remotes | 2.5.0 | CRAN |
| data.table | 1.18.4 | CRAN |
| dplyr | 1.2.1 | CRAN |
| ggplot2 | 4.0.3 | CRAN |
| **LAVA** | **0.1.5** | `remotes::install_github("josefin-werme/LAVA")` |
| **PLACO / PLACO+** | **commit 3ba3cae1d323ad117fb4540e620bcefa79f70663** | vendored source, not a package |
| LDSC | CBIIT fork, branch `ldsc39` | symlink to shared install |
| Python (LDSC env) | 3.9.23 | conda env |

R packages live in the project-local library `.r-lib/` (gitignored). Set
`export R_LIBS_USER="$PWD/.r-lib"` before any R invocation.

## Vendored PLACO

PLACO is distributed as plain R source, not an R package — `install_github`
fails on it because the repository has no DESCRIPTION. Vendored at a pinned
commit instead:

| file | bytes | sha256 (first 16) |
|---|---|---|
| `vendor/placo/PLACO_v0.2.0.R` | 7515 | fb684a8ed88f27dd |
| `vendor/placo/PLACO_v0.2.0_example.R` | 1754 | 2b46946d354c2f51 |

`placo.plus` — the sample-overlap-robust variant the project requires — is
present and exported.

## Verification result

`bash scripts/21_verify_phase23_tools.sh` → **TOOLCHAIN_VERIFIED**

- LAVA: `process.input`, `run.univ`, `run.bivar`, `read.loci` all present.
- PLACO+ calibration on 5,000 simulated null Z pairs (a SOFTWARE test, not a
  scientific result): `var.placo` = 1.0541, 0.9958 (expected ≈1 under the
  null); 200/200 p-values finite and within [0,1]; median p = 0.506
  (expected ≈0.5). Correctly calibrated.

## Not yet installed

**Octave** — required for pleioFDR (conjunctional FDR). Not installed in this
run. The pleioFDR repository's Octave support was last tested on Octave 4.0.2,
so a controlled reproduction of its chromosome-21 demonstration is required
before any output is trusted. Until that test passes, conjunctional FDR is
**pending**, not failed, and must not be reported as run.

## LD reference

`eur_w_ld_chr` (1000G Phase 3 EUR, HapMap3 SNPs), 48 files, plus
`w_hm3.snplist`. LAVA requires its own LD-block reference, which is **not yet
downloaded** — see BLOCKERS.md.
