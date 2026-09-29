#!/usr/bin/env Rscript
# Numerical-only audit of the pre-outcome UKB BCOR-derived regional matrices.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: script <LD manifest TSV> <output TSV>")
tasks <- utils::read.delim(args[[1]], check.names=FALSE)
if (!nrow(tasks) || any(tasks$status != "READY")) stop("Expected complete LD extraction manifest")
rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i,]
  x <- readRDS(task$ld_rds)$ld
  if (nrow(x) != ncol(x) || !all(is.finite(x))) stop("Malformed matrix")
  y <- .995*x
  diag(y) <- 1
  chol_ok <- !inherits(try(chol(y), silent=TRUE), "try-error")
  data.frame(analysis_id=task$analysis_id, chromosome=task$chromosome,
             n_snps=nrow(x), raw_max_abs_r=max(abs(x)),
             raw_above_one=sum(abs(x[upper.tri(x)]) > 1),
             raw_asymmetry=max(abs(x-t(x))),
             raw_diag_error=max(abs(diag(x)-1)),
             lambda=0.005, regularized_max_abs_r=max(abs(y)),
             regularized_cholesky=if (chol_ok) "PASS" else "FAIL",
             valid_for_v1_susie=if (chol_ok && max(abs(y)) <= 1+1e-8) "YES" else "NO",
             stringsAsFactors=FALSE)
})
result <- do.call(rbind, rows)
utils::write.table(result, args[[2]], sep="\t", row.names=FALSE,
                   col.names=TRUE, quote=FALSE)
cat("LD matrices",nrow(result),"valid under v1 lock",sum(result$valid_for_v1_susie=="YES"),"\n")
