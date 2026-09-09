#!/usr/bin/env Rscript
# Run on the intended production host, not assumed to exist in this delivery.
if (!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite in the reviewed environment first")
packages <- c("data.table","jsonlite","R.utils","susieR","coloc","LAVA","GenomicSEM","TwoSampleMR")
versions <- setNames(lapply(packages,function(p) {
  if (requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NULL
}),packages)
cat(jsonlite::toJSON(list(R=as.character(getRversion()),platform=R.version$platform,
    packages=versions,note="Record actual installed versions; this is not a tested environment lock."),
    auto_unbox=TRUE,pretty=TRUE,null="null"))
cat("\n")
