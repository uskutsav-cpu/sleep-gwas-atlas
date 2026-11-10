# Synthetic input only; no NHATS person-level records.
grid <- expand.grid(walking = c(1L, 2L, -1L), vigorous = c(1L, 2L, -1L))
grid$official_r <- with(grid, ifelse((vigorous == 2) & (walking == 2), 1, ifelse((vigorous < 0) & (walking < 0), NA, 0)))
grid$known_answer_rule <- with(grid, ifelse(walking == 1 | vigorous == 1, 0, ifelse(walking == 2 & vigorous == 2, 1, NA)))
grid$description <- ifelse(is.na(grid$official_r), 'missing', ifelse(grid$official_r == 1, 'low activity', 'not low activity'))
grid$differs_from_known_answer <- xor(is.na(grid$official_r), is.na(grid$known_answer_rule)) | (!is.na(grid$official_r) & !is.na(grid$known_answer_rule) & grid$official_r != grid$known_answer_rule)
write.csv(grid, 'outputs/stage6_nhats_activity_synthetic_truth_table.csv', row.names=FALSE, na='NA')
stopifnot(sum(grid$differs_from_known_answer) == 2, grid$official_r[grid$walking == 2 & grid$vigorous == -1] == 0)
print(grid)
