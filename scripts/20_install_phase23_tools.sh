#!/usr/bin/env bash
# Install and pin the Phase 2/3 R toolchain into a project-local library.
#
# Absence of software is not an external blocker: R installs via Homebrew and
# LAVA/PLACO install from their official GitHub repositories. This script is
# idempotent -- rerunning it skips what is already present.
set -uo pipefail
cd "$(dirname "$0")/.."
export R_LIBS_USER="$PWD/.r-lib"
mkdir -p "$R_LIBS_USER" environment

Rscript - <<'RS'
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
repos <- "https://cloud.r-project.org"
need <- function(p) !requireNamespace(p, quietly=TRUE)
for (p in c("remotes","data.table","dplyr","ggplot2","RcppArmadillo")) {
  if (need(p)) { cat("installing", p, "\n"); install.packages(p, lib=lib, repos=repos) }
}
if (need("LAVA"))  remotes::install_github("josefin-werme/LAVA", lib=lib, upgrade="never")
if (need("PLACO")) remotes::install_github("RayDebashree/PLACO", lib=lib, upgrade="never")
cat("\n--- installed ---\n")
for (p in c("remotes","data.table","dplyr","ggplot2","LAVA","PLACO")) {
  cat(sprintf("%-16s %s\n", p,
      if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING"))
}
RS
