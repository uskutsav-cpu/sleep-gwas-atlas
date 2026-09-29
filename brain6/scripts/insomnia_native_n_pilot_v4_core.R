# Input triage for the receipt-bound insomnia pilot. Zero rows are a biological
# no-data state; a missing or malformed shard remains an execution error.
check_insomnia_shard <- function(path, expected_rows) {
  if (!file.exists(path)) stop('MISSING_SUMSTATS_FILE', call. = FALSE)
  if (length(expected_rows) != 1L || is.na(expected_rows) ||
      expected_rows < 0 || expected_rows != as.integer(expected_rows))
    stop('INVALID_EXPECTED_ROW_COUNT', call. = FALSE)
  con <- gzfile(path, open = 'rt')
  on.exit(close(con))
  lines <- tryCatch(readLines(con, n = 2L, warn = FALSE),
                    error = function(e) stop(paste0('UNREADABLE_SUMSTATS_FILE: ', conditionMessage(e)), call. = FALSE))
  if (length(lines) < 1L || lines[[1L]] != 'SNP\tA1\tA2\tZ\tN')
    stop('INVALID_SUMSTATS_HEADER', call. = FALSE)
  if (expected_rows == 0L) {
    if (length(lines) != 1L) stop('ZERO_ROW_RECEIPT_HAS_DATA', call. = FALSE)
    return('NOT_RUN_NO_SNPS')
  }
  if (length(lines) < 2L || length(strsplit(lines[[2L]], '\t', fixed = TRUE)[[1L]]) != 5L)
    stop('NONZERO_ROW_RECEIPT_HAS_NO_VALID_DATA', call. = FALSE)
  'PROCEED'
}

capture_lava_result <- function(fn) {
  tryCatch(list(ok = TRUE, value = fn()),
           error = function(e) list(ok = FALSE, reason = substr(conditionMessage(e), 1L, 500L)))
}
