#!/usr/bin/env Rscript
# Install into a NEW isolated library only. This is setup, not statistical validation.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=1L) stop("Usage: Rscript bootstrap_native.R NEW_LIBRARY_DIRECTORY")
lib <- path.expand(a[[1]])
if(dir.exists(lib) && length(list.files(lib,all.files=TRUE,no..=TRUE))) stop("Library is not empty; do not change a running project's packages")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lib <- normalizePath(lib)
.libPaths(c(lib,.Library))
Sys.setenv(R_LIBS_USER=lib)
options(repos=c(CRAN="https://cloud.r-project.org"),Ncpus=1L)
cran <- c("remotes","jsonlite","data.table","R.utils","susieR","coloc","renv")
install.packages(cran,lib=lib,dependencies=NA)
if(!all(vapply(cran,requireNamespace,logical(1),quietly=TRUE))) stop("Missing CRAN dependency; setup failed")
pins <- c(LAVA="josefin-werme/LAVA@e729a245f7b6923967a96804fbf5246eadf2d6c6",
          GenomicSEM="GenomicSEM/GenomicSEM@6b65ca5db39fdade08b0d811477be1cdd57b5039",
          TwoSampleMR="MRCIEU/TwoSampleMR@c14776b89056e88422d9dff92b517b870a1339a4")
for(p in names(pins)) {
  remotes::install_github(pins[[p]],lib=lib,dependencies=NA,upgrade="never",build_vignettes=FALSE)
  if(!requireNamespace(p,quietly=TRUE)) stop(paste("Install failed:",p))
  desc <- packageDescription(p,lib.loc=lib)
  expected <- sub(".*@","",pins[[p]])
  if(is.null(desc$RemoteSha) || desc$RemoteSha!=expected) stop(paste("Wrong GitHub commit:",p))
}
ip <- installed.packages(lib.loc=lib)
write.table(ip[,intersect(c("Package","Version","Built"),colnames(ip)),drop=FALSE],
            file.path(lib,"INSTALLED_VERSIONS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(lib,"SESSION_INFO.txt"))
cat("Isolated R packages installed. Run native smoke tests and freeze an environment lock BEFORE real analysis.\n")
cat("No GWAS, PLINK, LDSC, LD reference, or MATLAB installation was performed.\n")
