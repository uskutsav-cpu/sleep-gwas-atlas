#!/usr/bin/env Rscript
# Phase 3: PLACO+ cross-trait pleiotropic scan.
#
# Runs on the LDSC-munged sumstats already produced in Phase 1, so no new
# reference download is required. Input per trait is SNP / A1 / A2 / Z / N.
#
#   Rscript scripts/30_phase3_placo.R <sleep_trait> <disease_trait>
#
# PLACO+ (placo.plus) is used rather than legacy PLACO because these traits are
# correlated and UK Biobank sample overlap is plausible throughout this project.
#
# Effect alleles are aligned explicitly: Z for the disease trait is flipped
# where its A1 matches the sleep trait's A2. Strand-ambiguous variants (A/T,
# C/G) cannot be resolved by allele codes alone and are removed.
suppressMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  library(data.table)
})
source("vendor/placo/PLACO_v0.2.0.R")

args <- commandArgs(trailingOnly = TRUE)
t1 <- args[1]; t2 <- args[2]
# PRE-DECLARED SCAN SCOPE. Genome-wide PLACO+ is 337 min/pair on this machine
# (measured; see results/blockers.tsv). The scan is therefore restricted to
# variants reaching P < SCREEN_P in AT LEAST ONE trait -- declared BEFORE
# looking at any result, and recorded in the provenance row. This is a
# pleiotropy SCREEN, not a genome-wide scan, and is labelled as such
# everywhere it is reported.
SCREEN_P <- as.numeric(Sys.getenv("PLACO_SCREEN_P", "1e-3"))
outdir <- "results/phase3"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tag <- paste0(t1, "__", t2)

rd <- function(t) {
  d <- fread(cmd = paste("gunzip -c", shQuote(sprintf("data/munged/%s.sumstats.gz", t))))
  d <- d[!is.na(Z) & !is.na(N)]
  setnames(d, c("A1", "A2", "Z", "N"), c("a1", "a2", "z", "n"))
  unique(d, by = "SNP")
}

d1 <- rd(t1); d2 <- rd(t2)
m <- merge(d1, d2, by = "SNP", suffixes = c(".1", ".2"))
n_merged <- nrow(m)

# align effect alleles
same <- m$a1.1 == m$a1.2 & m$a2.1 == m$a2.2
flip <- m$a1.1 == m$a2.2 & m$a2.1 == m$a1.2
m <- m[same | flip]
m[flip[same | flip], z.2 := -z.2]
n_aligned <- nrow(m)

# strand-ambiguous cannot be resolved from allele codes alone
amb <- paste0(m$a1.1, m$a2.1) %in% c("AT", "TA", "CG", "GC")
m <- m[!amb]
n_unamb <- nrow(m)

# PLACO's documented guidance: remove variants with extreme squared Z
zsq_cut <- qchisq(1e-4, df = 1, lower.tail = FALSE)
keep <- !(m$z.1^2 > 80 | m$z.2^2 > 80)
n_extreme <- sum(!keep)
m <- m[keep]

# variance and correlation nuisance parameters are estimated on the FULL
# merged set, as the method requires -- only the tested set is restricted
Zfull <- as.matrix(m[, .(z.1, z.2)]); colnames(Zfull) <- c("Z1", "Z2")
Pfull <- 2 * pnorm(-abs(Zfull)); colnames(Pfull) <- c("P1", "P2")
vt  <- var.placo(Zfull, Pfull, p.threshold = 1e-4)
rho <- cor.pearson(Zfull, Pfull, p.threshold = 1e-4, returnMatrix = FALSE)

n_full <- nrow(m)
sel <- (Pfull[, 1] < SCREEN_P) | (Pfull[, 2] < SCREEN_P)
m <- m[sel]
Z <- as.matrix(m[, .(z.1, z.2)]); colnames(Z) <- c("Z1", "Z2")
P <- 2 * pnorm(-abs(Z)); colnames(P) <- c("P1", "P2")
vt_unused <- var.placo(Z, P, p.threshold = 1e-4)
# placo.plus is per-SNP (Z is a length-2 vector), so it must be applied
# row-wise. It performs numerical integration per variant, which dominates
# runtime -- benchmark first, then run genome-wide.
pv <- vapply(seq_len(nrow(Z)),
             function(i) placo.plus(Z[i, ], VarZ = vt, CorZ = rho)$p.placo,
             numeric(1))
out <- data.table(SNP = m$SNP, CHR = NA_integer_, A1 = m$a1.1, A2 = m$a2.1,
                  Z1 = Z[, 1], Z2 = Z[, 2], N1 = m$n.1, N2 = m$n.2,
                  P_placo = as.numeric(pv))
out <- out[is.finite(P_placo)]
out[, direction := ifelse(sign(Z1) == sign(Z2), "same", "opposite")]
fwrite(out[order(P_placo)], sprintf("%s/placo_%s.tsv.gz", outdir, tag), sep = "\t")

gw <- out[P_placo < 5e-8]
lam <- median(qchisq(out$P_placo, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)

prov <- data.table(
  pair = tag, sleep_trait = t1, disease_trait = t2,
  n_merged = n_merged, n_after_qc = n_full, screen_p = SCREEN_P, scan_scope = "SCREEN_not_genome_wide", n_allele_aligned = n_aligned, n_unambiguous = n_unamb,
  n_extreme_removed = n_extreme, n_tested = nrow(out),
  varZ1 = vt[1], varZ2 = vt[2], cor_Z = rho,
  lambda_gc = lam, n_genome_wide = nrow(gw),
  n_same_dir = sum(gw$direction == "same"), n_opp_dir = sum(gw$direction == "opposite"))
fwrite(prov, sprintf("%s/prov_%s.tsv", outdir, tag), sep = "\t")

cat(sprintf("%-34s tested=%8d  lambda=%.3f  gw=%4d (same %d / opp %d)\n",
            tag, nrow(out), lam, nrow(gw),
            sum(gw$direction == "same"), sum(gw$direction == "opposite")))
