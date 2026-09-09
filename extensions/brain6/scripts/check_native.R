#!/usr/bin/env Rscript
# Does not install anything. Records native dependencies and parses all adapters.
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(normalizePath(script))
for(f in list.files(root,pattern="[.]R$",full.names=TRUE)) { parse(f); cat("PARSE PASS",basename(f),"\n") }
packages <- c("data.table","jsonlite","R.utils","susieR","coloc","LAVA","GenomicSEM","TwoSampleMR")
missing <- character()
for(p in packages) {
 if(requireNamespace(p,quietly=TRUE)) cat(p,as.character(packageVersion(p)),"\n")
 else { cat(p,"MISSING\n");missing <- c(missing,p) }
}
if(length(missing)) quit(status=2)
