# Read public, aggregate corrected-p matrix only; no individual records.
data <- readRDS('work/hpp_prs_corrected.Rds')
dict <- readRDS('work/hpp_prs_names.Rds')
scores <- c('Sleeplessness / insomnia', 'Non-cancer illness code, self-reported: insomnia',
            'Sleep duration', 'Nap during day', 'Morning/evening person (chronotype)')
stopifnot(all(scores %in% names(data)), nrow(data) == 730,
          data$Phenotype[1] == 'prs')
out <- do.call(rbind, lapply(scores, function(score) {
  v <- suppressWarnings(as.numeric(data[[score]][-1]))
  code <- unique(dict$code[dict$name == score])
  data.frame(score=score, score_code=paste(code,collapse=';'),
             rows=length(v), nonmissing=sum(!is.na(v)),
             less_than_005=sum(v<.05,na.rm=TRUE), minimum=min(v,na.rm=TRUE),
             maximum=max(v,na.rm=TRUE))
}))
write.csv(out, 'outputs/stage6_hpp_prs_matrix_summary.csv', row.names=FALSE)
print(out)
