# Independent base-R reconstruction directly from the original GEO count rows.
# Prints aggregate agreement only; no effect distributions or model performance.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root <- normalizePath(args[[1]])
source_dir <- file.path(root,"outputs/sleep_frailty_ml_channels/source_candidates")
dest <- file.path(root,"outputs/sleep_frailty_ml_measurement")
raw <- read.delim(gzfile(file.path(source_dir,"GSE87711_RBC_MPRA_minP_raw.txt.gz")),stringsAsFactors=FALSE,check.names=FALSE)
cases <- read.csv(file.path(source_dir,"Ulirsch_candidate_strict_no_history_or_published_pretraining.csv"),stringsAsFactors=FALSE,check.names=FALSE)
released <- read.csv(gzfile(file.path(dest,"Ulirsch_sealed_count_measurements.csv.gz")),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(cases)==nrow(released))
columns <- c("K562_minP_DNA1","K562_minP_DNA2",paste0("K562_CTRL_minP_RNA",1:6))
key_fields <- c("chr","pos","ref","alt","type","bot","top","clean")
keys <- do.call(paste,c(raw[key_fields],sep="|"))
values <- as.matrix(raw[columns])
stopifnot(all(is.finite(values)),all(values>=0),all(values==round(values)))
aggregated <- rowsum(values,keys,reorder=FALSE)
counts <- table(keys)
construct_key <- function(id) {
  match <- regexec("^([^:]+):([0-9]+):([^:]+):([^_]+)_(Ref|Mut)_Bot_([^_]+)_Top_([^_]+)(_mod)?$",id)
  parts <- regmatches(id,match)[[1]]
  stopifnot(length(parts)==9L)
  paste(c(parts[2:8],if (parts[9]=="_mod") "mod" else "var"),collapse="|")
}
truth <- function(value) tolower(as.character(value))=="true"
maximum_beta <- maximum_se <- maximum_empirical <- maximum_floor <- 0
eligible_n <- 0L
for (i in seq_len(nrow(cases))) {
  akey <- construct_key(cases$source_reference_construct[i])
  bkey <- construct_key(cases$source_alternate_construct[i])
  stopifnot(counts[[akey]]==14L,counts[[bkey]]==14L)
  ref <- aggregated[akey,];alt <- aggregated[bkey,]
  if (!truth(cases$GRCh38_reference_is_source_Ref[i])) { temp <- ref;ref <- alt;alt <- temp }
  stopifnot(all(ref==unlist(released[i,paste0("genomic_REF_",columns)],use.names=FALSE)),all(alt==unlist(released[i,paste0("genomic_ALT_",columns)],use.names=FALSE)))
  dna_ok <- all(ref[1:2]>=10) && all(alt[1:2]>=10)
  rna_ok <- sum(ref[3:8])+sum(alt[3:8])>0
  eligible <- dna_ok && rna_ok
  stopifnot(eligible==truth(released$eligible_point_contrast[i]))
  if (!eligible) {stopifnot(is.nan(released$signed_genomic_ALT_minus_REF_log2_ratio[i]));next}
  eligible_n <- eligible_n+1L
  ratio <- log2(alt+0.5)-log2(ref+0.5)
  beta <- mean(ratio[3:8])-mean(ratio[1:2])
  empirical <- var(ratio[3:8])/6+var(ratio[1:2])/2
  variance <- trigamma(alt+0.5)+trigamma(ref+0.5)
  floor <- (sum(variance[3:8])/36+sum(variance[1:2])/4)/log(2)^2
  se <- sqrt(max(empirical,floor))
  maximum_beta <- max(maximum_beta,abs(beta-released$signed_genomic_ALT_minus_REF_log2_ratio[i]))
  maximum_se <- max(maximum_se,abs(se-released$working_SE[i]))
  maximum_empirical <- max(maximum_empirical,abs(empirical-released$empirical_variance[i]))
  maximum_floor <- max(maximum_floor,abs(floor-released$count_variance_floor[i]))
}
stopifnot(maximum_beta<1e-9,maximum_se<1e-8,maximum_empirical<1e-9,maximum_floor<1e-8)
receipt <- data.frame(metric=c("passed","candidate_designs","eligible_point_contrasts","maximum_beta_difference","maximum_SE_difference","maximum_empirical_variance_difference","maximum_count_variance_floor_difference"),value=c(1,nrow(cases),eligible_n,maximum_beta,maximum_se,maximum_empirical,maximum_floor))
write.table(receipt,file.path(dest,"Ulirsch_independent_R_reconstruction_receipt.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("INDEPENDENT R COUNT RECONSTRUCTION PASSED; designs",nrow(cases),"eligible",eligible_n,"maximum beta difference",maximum_beta,"maximum SE difference",maximum_se,"\n")
