#!/usr/bin/env Rscript
# Exploratory single-signal GWAS–QTL coloc for pre-gated candidate contexts.

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
out <- file.path(repo, "brain6/results/brain6_exploratory_functional_v1")
gate <- read.delim(file.path(out, "coloc_gate.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
eligible <- gate[gate$status == "ELIGIBLE_SINGLE_SIGNAL_ABF", , drop = FALSE]
if (nrow(eligible) == 0L) stop("No preregistered eligible coloc contexts")
lib <- Sys.getenv("BRAIN6_R_LIB", "/private/tmp/brain6-r-lib")
.libPaths(c(lib, .libPaths()))
if (!requireNamespace("coloc", quietly = TRUE)) stop("Pinned coloc package unavailable")
if (as.character(packageVersion("coloc")) != "5.2.3") stop("Unexpected coloc version")
if (file.exists(file.path(out, "coloc_results.tsv"))) stop("Completed coloc result exists")

rows <- vector("list", nrow(eligible))
for (i in seq_len(nrow(eligible))) {
  x <- eligible[i, ]
  path <- file.path(repo, x$input_path)
  d <- read.delim(gzfile(path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(d) != x$n_unique_shared_variants || anyDuplicated(d$variant)) stop("Input count or variant identity mismatch")
  if (any(!is.finite(d$gwas_beta) | !is.finite(d$gwas_se) | !is.finite(d$qtl_beta) |
          !is.finite(d$qtl_se) | !is.finite(d$qtl_maf) | !is.finite(d$qtl_n))) stop("Invalid source statistics")
  gwas <- list(snp = d$variant, beta = d$gwas_beta, varbeta = d$gwas_se^2, type = "cc")
  qtl <- list(snp = d$variant, beta = d$qtl_beta, varbeta = d$qtl_se^2, type = "quant",
              MAF = d$qtl_maf, N = d$qtl_n)
  fit <- coloc::coloc.abf(gwas, qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  posterior <- fit$summary
  top <- fit$results[which.max(fit$results$SNP.PP.H4), , drop = FALSE]
  rows[[i]] <- data.frame(candidate_locus_id = x$candidate_locus_id, region_grch37 = x$region_grch37,
                         gwas_trait = x$gwas_trait, qtl_dataset_id = x$qtl_dataset_id,
                         qtl_context = x$qtl_context, molecular_trait_id = x$molecular_trait_id,
                         gene_id = x$gene_id, n_shared = posterior[["nsnps"]],
                         pp_h0 = posterior[["PP.H0.abf"]], pp_h1 = posterior[["PP.H1.abf"]],
                         pp_h2 = posterior[["PP.H2.abf"]], pp_h3 = posterior[["PP.H3.abf"]],
                         pp_h4 = posterior[["PP.H4.abf"]], top_shared_variant = top$snp,
                         top_shared_variant_conditional_h4_posterior = top$SNP.PP.H4,
                         status = "EXPLORATORY_SINGLE_SIGNAL_ABF_ONLY",
                         input_path = x$input_path, stringsAsFactors = FALSE)
}
result <- do.call(rbind, rows)
write.table(result, file.path(out, "coloc_results.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat("coloc version", as.character(packageVersion("coloc")), "eligible tests", nrow(result), "\n")
