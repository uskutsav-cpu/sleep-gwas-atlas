#!/usr/bin/env Rscript
# Calls upstream PLACO functions; no rewritten tail approximation.
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need(file.exists(cfg$placo_source), "Pinned upstream PLACO R source is missing")
source(cfg$placo_source)
need(cfg$method %in% c("PLACO_PLUS","PLACO"), "Unknown PLACO method")
need(is.numeric(cfg$extreme_z2) && cfg$extreme_z2>0, "Freeze the extreme-Z threshold")
need(is.numeric(cfg$p_threshold) && cfg$p_threshold>0 && cfg$p_threshold<1,
     "Freeze null parameter-estimation threshold")

# The statistic is Z1*Z2. When either Z is exactly zero, the observed
# statistic is exactly zero and its two-sided tail probability is exactly 1.
# The upstream density has an integrable K0 singularity at zero; numerical
# integration at that boundary may otherwise return NA or an out-of-range
# value. This exact boundary case is handled analytically, without changing
# the upstream PLACO+ calculation for nonzero statistics.
placo_plus_pvalue <- function(z, VarZ, CorZ, AbsTol) {
  if (any(z == 0)) return(1.0)
  placo.plus(z, VarZ=VarZ, CorZ=CorZ, AbsTol=AbsTol)$p.placo.plus
}

if (cfg$mode == "estimate") {
  # Use ALL genome-wide eligible variants, not just significant candidates.
  d <- read_dt(cfg$pair_file, c("Z1","Z2","P1","P2"))
  finite_cols(d,c("Z1","Z2","P1","P2"))
  keep <- d$Z1^2 <= cfg$extreme_z2 & d$Z2^2 <= cfg$extreme_z2
  d <- d[keep]
  need(nrow(d) >= cfg$min_variants, "Insufficient genome-wide parameter-estimation data")
  need(all(d$P1>=0 & d$P1<=1 & d$P2>=0 & d$P2<=1), "Invalid marginal P")
  null_var <- !(d$P1<cfg$p_threshold & d$P2<cfg$p_threshold)
  null_cor <- d$P1>=cfg$p_threshold & d$P2>=cfg$p_threshold
  need(sum(null_var)>=1000 && sum(null_cor)>=1000, "Too few null variants for stable global parameters")
  Z <- as.matrix(d[,.(Z1,Z2)]); P <- as.matrix(d[,.(P1,P2)])
  vz <- var.placo(Z,P,p.threshold=cfg$p_threshold)
  rho <- cor.pearson(Z,P,p.threshold=cfg$p_threshold,returnMatrix=FALSE)
  need(all(is.finite(vz)) && all(vz>0) && is.finite(rho) && abs(rho)<1,
       "Invalid PLACO global variance/correlation parameters")
  if (cfg$method == "PLACO") {
    need(isTRUE(cfg$uncorrelated_inputs_reviewed), "Original PLACO requires uncorrelated inputs review")
  }
  write_json(list(VarZ=unname(vz), CorZ=unname(rho), method=cfg$method,
                 p_threshold=cfg$p_threshold,extreme_z2=cfg$extreme_z2,
                 n_genomewide=nrow(d),n_null_var=sum(null_var),n_null_cor=sum(null_cor)),
             opath("parameters.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  finish(extra=list(mode="estimate",n_genomewide=nrow(d)))
} else if (cfg$mode == "chunk") {
  par <- fromJSON(cfg$parameters)
  need(par$method == cfg$method && par$extreme_z2 == cfg$extreme_z2 &&
       par$p_threshold == cfg$p_threshold, "Chunk/global parameter policy mismatch")
  d <- read_dt(cfg$chunk_file, c("SNP","CHR","BP","Z1","Z2"))
  finite_cols(d,c("CHR","BP","Z1","Z2"))
  need(!anyDuplicated(d$SNP), "Duplicate SNPs in chunk")
  d[, `:=`(P_PLACO=NA_real_,status="NUMERICAL_FAILURE")]
  exclusions <- d$Z1^2 > cfg$extreme_z2 | d$Z2^2 > cfg$extreme_z2
  d[exclusions,status:="EXCLUDED_EXTREME_Z"]
  errors <- list()
  for (i in which(!exclusions)) {
    ans <- tryCatch({
      if (cfg$method == "PLACO_PLUS") {
        placo_plus_pvalue(c(d$Z1[i],d$Z2[i]),VarZ=par$VarZ,CorZ=par$CorZ,AbsTol=cfg$absolute_tolerance)
      } else {
        placo(c(d$Z1[i],d$Z2[i]),VarZ=par$VarZ,AbsTol=cfg$absolute_tolerance)$p.placo
      }
    }, error=function(e) { errors[[length(errors)+1L]] <<- list(SNP=d$SNP[i],message=conditionMessage(e)); NA_real_ })
    # Negative/out-of-range numerical values are not clipped into significance.
    if (length(ans)==1L && is.finite(ans) && ans>=0 && ans<=1) {
      set(d,i,"P_PLACO",ans); set(d,i,"status","TESTED")
    }
  }
  fwrite(d,opath("chunk.tsv"),sep="\t",na="NA")
  write_json(errors,opath("numerical_errors.json"),pretty=TRUE,auto_unbox=TRUE)
  finish(extra=list(mode="chunk",n_rows=nrow(d),numerical_failures=sum(d$status=="NUMERICAL_FAILURE"),
                    note="Per-chunk completion is not pair-level QC acceptance; run collate-placo."))
} else stop("Mode must be estimate or chunk")
