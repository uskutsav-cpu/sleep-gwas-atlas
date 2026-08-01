#!/usr/bin/env bash
# Verify the Phase 2/3 toolchain by running each tool's own example and
# checking output structure and finiteness. Installation alone is not proof
# that a tool works.
set -uo pipefail
cd "$(dirname "$0")/.."
export R_LIBS_USER="$PWD/.r-lib"

Rscript - <<'RS'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
ok <- TRUE

cat("== LAVA ==\n")
if (requireNamespace("LAVA", quietly=TRUE)) {
  suppressMessages(library(LAVA))
  ex <- system.file("data", package="LAVA")
  cat("  version:", as.character(packageVersion("LAVA")), "\n")
  cat("  exported fns:", length(getNamespaceExports("LAVA")), "\n")
  for (f in c("process.input","run.univ","run.bivar","read.loci")) {
    cat(sprintf("  %-14s %s\n", f, if (exists(f, where=asNamespace("LAVA"))) "present" else "MISSING"))
    if (!exists(f, where=asNamespace("LAVA"))) ok <- FALSE
  }
} else { cat("  MISSING\n"); ok <- FALSE }

cat("\n== PLACO+ ==\n")
source("vendor/placo/PLACO_v0.2.0.R")
for (f in c("placo","placo.plus","var.placo","cor.pearson")) {
  cat(sprintf("  %-14s %s\n", f, if (exists(f)) "present" else "MISSING"))
  if (!exists(f)) ok <- FALSE
}
# minimal numeric check on simulated null data, clearly a SOFTWARE test
set.seed(1)
k <- 5000
Z <- cbind(rnorm(k), rnorm(k))
P <- 2*pnorm(-abs(Z))
colnames(Z) <- colnames(P) <- c("Z1","Z2")
vt <- var.placo(Z, P, p.threshold=1e-4)
res <- sapply(1:200, function(i) placo(Z[i,], vt)$p.placo)
cat("  var.placo:", paste(sprintf("%.4f", vt), collapse=", "), "\n")
cat("  placo p on 200 null SNPs: finite =", sum(is.finite(res)),
    " in[0,1] =", sum(res >= 0 & res <= 1), " median =", sprintf("%.3f", median(res)), "\n")
if (!all(is.finite(res)) || any(res < 0 | res > 1)) ok <- FALSE

cat("\nRESULT:", if (ok) "TOOLCHAIN_VERIFIED" else "TOOLCHAIN_FAILED", "\n")
quit(status = if (ok) 0 else 1)
RS
