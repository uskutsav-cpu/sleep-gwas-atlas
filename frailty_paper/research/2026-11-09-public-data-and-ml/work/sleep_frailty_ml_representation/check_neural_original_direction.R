# Independent direction algebra and full original-coefficient check in base R.
folder <- '/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/outputs/sleep_frailty_ml_representation/neural_source_trace'
code<-readLines(file.path(folder,'perform_bcalm_variants_cli.R'),warn=FALSE)
stopifnot(any(grepl('alt = grepl("alt", colnames(mpra))',code,fixed=TRUE)),
          any(grepl('topTable(mpralm_fit, coef = 2, number = Inf)',code,fixed=TRUE)))
names<-c('sample1_bc1_ref','sample1_bc1_alt','sample2_bc1_ref','sample2_bc1_alt','sample3_bc1_ref','sample3_bc1_alt')
design<-as.matrix(data.frame(intcpt=1,alt=grepl('alt',names)))
errors<-numeric()
for (baseline in c(-1,0,2)) for (delta in c(-.7,0,.5)) {
  response<-baseline+as.numeric(grepl('alt',names))*delta
  fitted<-qr.solve(design,response)
  errors<-c(errors,abs(fitted[2]-delta))
}
stopifnot(max(errors)<1e-12)
labels<-read.csv(gzfile(file.path(folder,'neural_original_signed_development.csv.gz')),check.names=FALSE)
original<-read.delim(file.path(folder,'NGN2_variant_bcalm_80k_bbmap_251128_bcalm_normalized_variant_map_names_dev.tsv'),check.names=FALSE)
stopifnot(!anyDuplicated(original$variant_id),!anyDuplicated(labels$observation_index),all(labels$split %in% c(0,1)))
index<-match(labels$source_variant_id,original$variant_id);stopifnot(!anyNA(index))
coefficient_error<-max(abs(original$logFC[index]-labels$signed_ALT_minus_REF_log2_effect))
stopifnot(coefficient_error<1e-12)
write.table(data.frame(passed=TRUE,all_recovered_rows_checked=nrow(labels),
          max_original_coefficient_error=coefficient_error,max_direction_algebra_error=max(errors),
          original_counts_refitted=FALSE),file.path(folder,'neural_direction_independent_R_receipt.tsv'),
          sep='\t',row.names=FALSE,quote=FALSE)
cat('NEURONAL ORIGINAL DIRECTION INDEPENDENT R CHECK PASS',nrow(labels),'\n')
