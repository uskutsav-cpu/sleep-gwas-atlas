#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: 58_extract_lava_ld.R REFERENCE_PREFIX CHROMOSOME VARIANT_ORDER OUT MIN_REFERENCE_N")
}

reference_prefix <- args[[1]]
chromosome <- as.integer(args[[2]])
order_path <- args[[3]]
out_path <- args[[4]]
minimum_reference_n <- as.integer(args[[5]])
if (!is.finite(chromosome) || chromosome < 1 || chromosome > 22 ||
    !is.finite(minimum_reference_n) || minimum_reference_n < 1) {
  stop("invalid chromosome or reference sample-size threshold")
}

suppressPackageStartupMessages(library(LAVA))
if (as.character(packageVersion("LAVA")) != "0.1.5") {
  stop("LAVA runtime differs from locked 0.1.5")
}
order <- read.delim(order_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(order), "SNP") || nrow(order) < 1 || anyDuplicated(order$SNP)) {
  stop("invalid fine-mapping variant-order file")
}
order$SNP <- tolower(order$SNP)
chromosome_prefix <- paste0(reference_prefix, "_chr", chromosome)
reference <- LAVA:::load.reference(LAVA:::check.reference(chromosome_prefix))
if (!is.finite(reference$sample.size) || reference$sample.size < minimum_reference_n) {
  stop("LAVA reference sample size is below the locked minimum")
}
loaded <- LAVA:::read.ld(reference, order$SNP, require.freq = TRUE)
if (loaded$mode != "ld" || !identical(rownames(loaded$ld), order$SNP) ||
    !identical(colnames(loaded$ld), order$SNP)) {
  stop("LAVA signed-LD order differs from the locked locus order")
}
matrix <- loaded$ld
if (any(!is.finite(matrix)) || max(abs(matrix - t(matrix))) > 1e-8 ||
    max(abs(diag(matrix) - 1)) > 1e-6 || max(abs(matrix)) > 1 + 1e-6) {
  stop("LAVA signed-LD matrix fails finite/symmetry/unit-diagonal/range QC")
}
temporary <- paste0(out_path, ".tmp")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(temporary, open = "wt", compression = 9)
write.table(
  data.frame(SNP = rownames(matrix), matrix, check.names = FALSE),
  connection, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA"
)
close(connection)
if (!file.rename(temporary, out_path)) stop("failed to atomically install LAVA signed-LD matrix")
cat(sprintf("LAVA_SIGNED_LD_OK chromosome=%d variants=%d reference_n=%d\n",
            chromosome, nrow(matrix), as.integer(reference$sample.size)))
