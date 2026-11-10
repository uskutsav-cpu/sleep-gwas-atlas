# Independent base-R coefficient join; read only the two required columns.
folder <- "outputs/sleep_frailty_ml_neuronal_activity"
metadata <- file.path(folder, "ngn2_element_activity_metadata_202602_bbmap_only_mpralib095_normalized_bcalm_dna_ge_1.tsv.gz")
fitted <- file.path(folder, "cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_bc10_toptreat_results_scrambled_NP_MK_result_095_dev.tsv.gz")
read_required <- function(path) {
  handle <- gzfile(path, "rt")
  header <- strsplit(readLines(handle, n=1), "\t", fixed=TRUE)[[1]]
  close(handle)
  classes <- rep("NULL", length(header))
  classes[header == "name"] <- "character"
  classes[header == "logFC"] <- "numeric"
  stopifnot(sum(header == "name") == 1, sum(header == "logFC") == 1)
  read.table(gzfile(path), header=TRUE, sep="\t", quote="", comment.char="",
             check.names=FALSE, colClasses=classes, na.strings=c("NA", "NaN", ""))
}
a <- read_required(metadata)
b <- read_required(fitted)
stopifnot(!anyDuplicated(b$name))
good <- is.finite(a$logFC)
j <- match(a$name[good], b$name)
stopifnot(!anyNA(j))
error <- max(abs(a$logFC[good] - b$logFC[j]))
stopifnot(error < 1e-10)
receipt <- data.frame(passed=TRUE, metadata_rows=nrow(a), finite_original_activity_rows=sum(good),
    missing_activity_rows=sum(!good), all_finite_coefficients_matched=TRUE,
    maximum_absolute_coefficient_error=error, statistical_significance_filter_used=FALSE,
    original_counts_or_BCalm_model_refitted=FALSE, new_model_fitted=FALSE,
    reserved_study_outcomes_read=FALSE)
write.table(receipt, file.path(folder, "neuronal_activity_independent_R_receipt.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
print(receipt)
