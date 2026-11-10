# Independent source-table replay, using a deterministic sample per source.
root <- '/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo'
out <- file.path(root,'outputs/sleep_frailty_ml_direction')
path <- file.path(out,'original_paired_signed_development.csv.gz')
con <- gzfile(path,'rt'); header <- strsplit(readLines(con,n=1),',',fixed=TRUE)[[1]];close(con)
classes <- rep('NULL',length(header));names(classes)<-header
for (key in c('observation_index','signed_ALT_minus_REF_log2_effect','split')) classes[key]<-'numeric'
for (key in c('original_source_path','source_row_ID','source_variant_GRCh37','construct_genomic_orientation')) classes[key]<-'character'
data <- read.csv(gzfile(path),colClasses=classes,check.names=FALSE)
stopifnot(all(data$split %in% c(0,1)),!anyDuplicated(data$observation_index))
records <- list()
for (source in unique(data$original_source_path)) {
  rows <- data[data$original_source_path==source,,drop=FALSE]
  selected <- unique(c(round(seq(1,nrow(rows),length.out=min(256,nrow(rows)))),
                       which(rows$construct_genomic_orientation=='-')[1]))
  selected <- selected[!is.na(selected)]
  rows <- rows[selected,,drop=FALSE]
  stopifnot(!anyDuplicated(rows$source_row_ID))
  found <- rep(FALSE,nrow(rows)); maxerror<-0
  handle <- gzfile(source,'rt')
  fields <- strsplit(readLines(handle,n=1),'\t',fixed=TRUE)[[1]]
  at <- match(c('ID','SNP','Log2Skew','window','strand','haplotype'),fields)
  stopifnot(!anyNA(at))
  repeat {
    lines <- readLines(handle,n=8192,warn=FALSE)
    if (!length(lines)) break
    parts <- strsplit(lines,'\t',fixed=TRUE)
    ids <- vapply(parts,function(x) x[at[1]],character(1))
    lookup <- match(ids,rows$source_row_ID)
    for (j in which(!is.na(lookup))) {
      i<-lookup[j];x<-parts[[j]]
      # Source IDs can recur in other windows: match the exact retained design.
      if (x[at[4]]!='center'||x[at[5]]!='fwd'||x[at[6]]!='ref') next
      stopifnot(!found[i],x[at[2]]==rows$source_variant_GRCh37[i])
      error <- abs(as.numeric(x[at[3]])-rows$signed_ALT_minus_REF_log2_effect[i])
      stopifnot(error<1e-12);maxerror<-max(maxerror,error);found[i]<-TRUE
    }
  }
  close(handle);stopifnot(all(found))
  records[[length(records)+1]]<-data.frame(source=source,checked_rows=nrow(rows),max_coefficient_error=maxerror,passed=TRUE)
}
write.table(do.call(rbind,records),file.path(out,'signed_source_independent_R_receipt.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('SIGNED SOURCE ROWS INDEPENDENT R CHECK PASS',sum(vapply(records,function(r) r$checked_rows,numeric(1))),'\n')
