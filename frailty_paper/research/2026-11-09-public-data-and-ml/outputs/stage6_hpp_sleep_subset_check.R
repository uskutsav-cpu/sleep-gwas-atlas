# Explore public HPP association-file coverage; no participant data or prediction.
models <- list(
  `15` = read.csv('outputs/unvalidated_15_variant_score_prototype.csv'),
  `22` = read.csv('outputs/unvalidated_22_shared_direction_score.csv')
)
files <- c(SleepLatency='work/hpp_sleep_latency.Rds', TotalSleepTime='work/hpp_total_sleep_time.Rds')
results <- list()
for (trait in names(files)) {
  x <- readRDS(files[[trait]])
  stopifnot(all(c('ID','P','OBS_CT') %in% names(x)), nrow(x)>0)
  for (model in names(models)) {
    ids <- models[[model]]$variant
    matches <- x[x$ID %in% ids, , drop=FALSE]
    results[[length(results)+1]] <- data.frame(trait=trait, model=model,
      rows=nrow(x), minimum_p=min(x$P, na.rm=TRUE), maximum_p=max(x$P, na.rm=TRUE),
      obs_min=min(x$OBS_CT, na.rm=TRUE), obs_max=max(x$OBS_CT, na.rm=TRUE),
      candidate_ids=length(unique(ids)), matched_ids=length(unique(matches$ID)))
  }
}
result <- do.call(rbind, results)
write.csv(result, 'outputs/stage6_hpp_sleep_subset_inventory.csv', row.names=FALSE)
print(result)
