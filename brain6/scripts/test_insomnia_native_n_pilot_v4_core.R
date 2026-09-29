#!/usr/bin/env Rscript
source('brain6/scripts/insomnia_native_n_pilot_v4_core.R')
d <- tempfile('insomnia-v4-test-'); dir.create(d)
on.exit(unlink(d, recursive = TRUE), add = TRUE)
write_gz <- function(name, lines) {
  path <- file.path(d, name)
  con <- gzfile(path, open = 'wt'); writeLines(lines, con); close(con)
  path
}
expect_error <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = function(e) conditionMessage(e))
  stopifnot(!is.null(error), grepl(pattern, error, fixed = TRUE))
}
header <- 'SNP\tA1\tA2\tZ\tN'
empty <- write_gz('empty.gz', header)
valid <- write_gz('valid.gz', c(header, 'rs1\tA\tC\t1.2\t386533'))
bad_header <- write_gz('bad-header.gz', c('SNP\tA1\tA2\tZ', 'rs1\tA\tC\t1.2'))
stopifnot(check_insomnia_shard(empty, 0L) == 'NOT_RUN_NO_SNPS')
stopifnot(check_insomnia_shard(valid, 1L) == 'PROCEED')
expect_error(check_insomnia_shard(file.path(d, 'missing.gz'), 0L), 'MISSING_SUMSTATS_FILE')
expect_error(check_insomnia_shard(bad_header, 1L), 'INVALID_SUMSTATS_HEADER')
expect_error(check_insomnia_shard(empty, 1L), 'NONZERO_ROW_RECEIPT_HAS_NO_VALID_DATA')
expect_error(check_insomnia_shard(valid, 0L), 'ZERO_ROW_RECEIPT_HAS_DATA')
expect_error(check_insomnia_shard(valid, -1L), 'INVALID_EXPECTED_ROW_COUNT')
ok <- capture_lava_result(function() 3L)
bad <- capture_lava_result(function() stop('real LAVA failure'))
stopifnot(ok$ok, ok$value == 3L, !bad$ok, bad$reason == 'real LAVA failure')
cat('INSOMNIA_V4_CORE_TESTS_PASS\n')
