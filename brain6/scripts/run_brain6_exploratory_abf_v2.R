#!/usr/bin/env Rscript
# Pinned single-signal ABF inference for all source-admitted Brain6 loci.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("usage: script <repository-output-dir> <SSD-detail-dir> <private-R-library>")
.libPaths(c(args[[3]], .libPaths()))
suppressPackageStartupMessages(library(coloc))
if (as.character(packageVersion("coloc")) != "5.2.3") stop("Expected coloc 5.2.3")
out <- normalizePath(args[[1]], mustWork=TRUE)
detail <- normalizePath(args[[2]], mustWork=TRUE)
trait_input <- utils::read.delim(file.path(out,"fine_mapping_input_manifest.tsv"),check.names=FALSE)
pair_input <- utils::read.delim(file.path(out,"trait_coloc_input_manifest.tsv"),check.names=FALSE)
if (nrow(trait_input)!=50L || nrow(pair_input)!=25L ||
    anyDuplicated(trait_input$analysis_id) || anyDuplicated(pair_input$candidate_locus_id))
  stop("Expected exactly 50 unique candidate-trait and 25 pair inputs")

sha256 <- function(path) {
  value <- system2("/usr/bin/shasum", c("-a","256",shQuote(path)), stdout=TRUE)
  if (length(value)!=1L || !grepl("^[a-f0-9]{64} ",value)) stop("SHA-256 failed: ",path)
  substr(value,1,64)
}
read_input <- function(path, expected_sha) {
  if (!file.exists(path) || sha256(path)!=expected_sha) stop("Input SHA-256 mismatch: ",path)
  x <- utils::read.delim(path,check.names=FALSE)
  if (!identical(names(x),c("SNP","CHR","POS","A1","A2","BETA","SE","N_AUDIT")) ||
      nrow(x)<500L || anyDuplicated(x$SNP) || any(!is.finite(x$BETA)) ||
      any(!is.finite(x$SE) | x$SE<=0)) stop("Malformed ABF input: ",path)
  x
}
dataset <- function(x) list(beta=x$BETA,varbeta=x$SE^2,snp=as.character(x$SNP),type="cc")
finemap <- function(x) {
  f <- coloc::finemap.abf(dataset(x),p1=1e-4)
  null <- f$SNP.PP[f$snp=="null"]
  if (length(null)!=1L || !is.finite(null) || null<0 || null>1) stop("Invalid ABF null posterior")
  variants <- f[f$snp!="null",c("snp","SNP.PP"),drop=FALSE]
  if (!setequal(variants$snp,x$SNP)) stop("ABF SNP mismatch")
  variants <- variants[match(x$SNP,variants$snp),,drop=FALSE]
  association <- 1-null
  variants$conditional_pip <- if (association>0) variants$SNP.PP/association else NA_real_
  if (association>0 && abs(sum(variants$conditional_pip)-1)>1e-6)
    stop("Conditional ABF PIP does not sum to one")
  list(null=null,association=association,pip=variants)
}
credible_set <- function(f) {
  if (f$association<0.8) return(data.frame())
  order <- order(-f$pip$conditional_pip,f$pip$snp)
  ordered <- f$pip[order,,drop=FALSE]
  cumulative <- cumsum(ordered$conditional_pip)
  cut <- which(cumulative>=0.95)[1]
  if (is.na(cut)) stop("95% credible set coverage unavailable")
  ans <- ordered[seq_len(cut),,drop=FALSE]
  ans$rank <- seq_len(cut)
  ans$cumulative_pip <- cumulative[seq_len(cut)]
  ans
}
posterior <- function(result,key) {
  value <- unname(result$summary[[key]])
  if (length(value)!=1L || !is.finite(value)) stop("Missing coloc posterior ",key)
  as.numeric(value)
}

trait_rows <- list(); cs_rows <- list(); pip_rows <- list(); errors <- character()
for (i in seq_len(nrow(trait_input))) {
  m <- trait_input[i,]
  row <- list(candidate_locus_id=m$candidate_locus_id,pair_id=m$pair_id,
              trait_id=m$trait_id,status=m$status,reason=m$reason,
              n_snps=NA_integer_,PP_null=NA_real_,PP_association=NA_real_,
              top_snp="",top_snp_conditional_pip=NA_real_,credible_set_size=NA_integer_,
              credible_set_cumulative_pip=NA_real_,
              analysis_label="EXPLORATORY_SINGLE_SIGNAL_ABF",
              input_sha256=m$input_sha256)
  if (m$status=="READY") {
    row <- tryCatch({
      x <- read_input(m$input_tsv,m$input_sha256)
      f <- finemap(x)
      cs <- credible_set(f)
      top <- which.max(f$pip$conditional_pip)
      row$n_snps <- nrow(x)
      row$PP_null <- f$null
      row$PP_association <- f$association
      row$top_snp <- f$pip$snp[top]
      row$top_snp_conditional_pip <- f$pip$conditional_pip[top]
      row$credible_set_size <- nrow(cs)
      row$credible_set_cumulative_pip <- if(nrow(cs)) tail(cs$cumulative_pip,1) else NA_real_
      row$status <- if(nrow(cs)) "EXPLORATORY_ABF_CREDIBLE_SET" else "NO_ASSOCIATION_SUPPORTED_CREDIBLE_SET"
      row$reason <- if(nrow(cs)) "" else "PP_ASSOCIATION_LT_0.80"
      if(nrow(cs)) {
        cs$candidate_locus_id <- m$candidate_locus_id
        cs$pair_id <- m$pair_id
        cs$trait_id <- m$trait_id
        cs_rows[[length(cs_rows)+1L]] <- cs
      }
      f$pip$candidate_locus_id <- m$candidate_locus_id
      f$pip$pair_id <- m$pair_id
      f$pip$trait_id <- m$trait_id
      pip_rows[[length(pip_rows)+1L]] <- f$pip
      row
    },error=function(e) {
      row$status <- "RUNTIME_FAILURE"
      row$reason <- conditionMessage(e)
      errors <<- c(errors,paste(m$analysis_id,conditionMessage(e)))
      row
    })
  }
  trait_rows[[length(trait_rows)+1L]] <- as.data.frame(row,stringsAsFactors=FALSE)
  cat("fine_map",m$analysis_id,row$status,"\n")
}

pair_rows <- list(); sensitivity_rows <- list()
for (i in seq_len(nrow(pair_input))) {
  m <- pair_input[i,]
  row <- list(candidate_locus_id=m$candidate_locus_id,pair_id=m$pair_id,
              status=m$status,reason=m$reason,n_shared=NA_integer_,
              PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_,
              PP.H4_conditional_H3H4=NA_real_,
              trait1_PP_association_shared=NA_real_,trait2_PP_association_shared=NA_real_,
              descriptive_support="NOT_ESTIMABLE",analysis_label="EXPLORATORY_SINGLE_SIGNAL_ABF",
              p12=1e-5,trait1_sha256=m$trait1_sha256,trait2_sha256=m$trait2_sha256)
  if(m$status=="READY") {
    row <- tryCatch({
      a <- read_input(m$trait1_input_tsv,m$trait1_sha256)
      b <- read_input(m$trait2_input_tsv,m$trait2_sha256)
      if(!identical(a$SNP,b$SNP) || !identical(a$A1,b$A1) || !identical(a$A2,b$A2))
        stop("Shared allele-orientation mismatch")
      fa <- finemap(a); fb <- finemap(b)
      result <- coloc::coloc.abf(dataset(a),dataset(b),p1=1e-4,p2=1e-4,p12=1e-5)
      ps <- sapply(paste0("PP.H",0:4,".abf"),function(k) posterior(result,k))
      if(abs(sum(ps)-1)>1e-6) stop("Coloc H0-H4 do not sum to one")
      row$n_shared <- nrow(a)
      for(j in 0:4) row[[paste0("PP.H",j)]] <- ps[j+1]
      row$PP.H4_conditional_H3H4 <- if(ps[4]+ps[5]>0) ps[5]/(ps[4]+ps[5]) else NA_real_
      row$trait1_PP_association_shared <- fa$association
      row$trait2_PP_association_shared <- fb$association
      supported <- ps[5]>=0.8 && is.finite(row$PP.H4_conditional_H3H4) &&
                   row$PP.H4_conditional_H3H4>=0.8 &&
                   fa$association>=0.8 && fb$association>=0.8
      row$descriptive_support <- if(supported) "ABF_MODEL_SUPPORT" else "ABF_MODEL_SUPPORT_NOT_MET"
      row$status <- "EXPLORATORY_TRAIT_COLOC_ABF"
      row$reason <- "SINGLE_CAUSAL_VARIANT_PER_TRAIT_ASSUMPTION_UNVERIFIED"
      for(prior in c(1e-6,5e-5)) {
        alt <- coloc::coloc.abf(dataset(a),dataset(b),p1=1e-4,p2=1e-4,p12=prior)
        sensitivity_rows[[length(sensitivity_rows)+1L]] <- data.frame(
          candidate_locus_id=m$candidate_locus_id,pair_id=m$pair_id,p12=prior,
          PP.H0=posterior(alt,"PP.H0.abf"),PP.H1=posterior(alt,"PP.H1.abf"),
          PP.H2=posterior(alt,"PP.H2.abf"),PP.H3=posterior(alt,"PP.H3.abf"),
          PP.H4=posterior(alt,"PP.H4.abf"))
      }
      row
    },error=function(e) {
      row$status <- "RUNTIME_FAILURE"
      row$reason <- conditionMessage(e)
      errors <<- c(errors,paste(m$candidate_locus_id,conditionMessage(e)))
      row
    })
  }
  pair_rows[[length(pair_rows)+1L]] <- as.data.frame(row,stringsAsFactors=FALSE)
  cat("coloc",m$candidate_locus_id,row$status,"\n")
}

write_tsv <- function(x,path) utils::write.table(x,path,sep="\t",row.names=FALSE,
                                                  col.names=TRUE,quote=FALSE,na="")
trait_result <- do.call(rbind,trait_rows)
pair_result <- do.call(rbind,pair_rows)
if(nrow(trait_result)!=50L || nrow(pair_result)!=25L) stop("Incomplete result family")
write_tsv(trait_result,file.path(out,"fine_mapping_50.tsv"))
write_tsv(pair_result,file.path(out,"trait_coloc_25.tsv"))
cs <- if(length(cs_rows)) do.call(rbind,cs_rows) else
  data.frame(snp=character(),SNP.PP=numeric(),conditional_pip=numeric(),rank=integer(),
             cumulative_pip=numeric(),candidate_locus_id=character(),pair_id=character(),trait_id=character())
write_tsv(cs,file.path(out,"fine_mapping_credible_sets.tsv"))
pip <- if(length(pip_rows)) do.call(rbind,pip_rows) else
  data.frame(snp=character(),SNP.PP=numeric(),conditional_pip=numeric(),
             candidate_locus_id=character(),pair_id=character(),trait_id=character())
write_tsv(pip,file.path(detail,"variant_pip.tsv"))
sens <- if(length(sensitivity_rows)) do.call(rbind,sensitivity_rows) else
  data.frame(candidate_locus_id=character(),pair_id=character(),p12=numeric(),
             PP.H0=numeric(),PP.H1=numeric(),PP.H2=numeric(),PP.H3=numeric(),PP.H4=numeric())
write_tsv(sens,file.path(out,"trait_coloc_prior_sensitivity.tsv"))
cat("fine_map_rows",nrow(trait_result),"coloc_rows",nrow(pair_result),
    "credible_set_variants",nrow(cs),"pip_variants",nrow(pip),"\n")
if(length(errors)) stop(length(errors)," runtime failures: ",paste(errors,collapse="; "))
