# Illustrative fixed-design linear-regression power for a single additional score.
# ncp is an approximation for a specified population partial R-squared; it is not empirical score performance.
alpha <- 0.05; k <- 10L
power_at <- function(n,r2) {df <- n-k-2; lambda <- (n-k-1)*r2/(1-r2); crit <- qf(1-alpha,1,df); pf(crit,1,df,ncp=lambda,lower.tail=FALSE)}
grid <- expand.grid(n=c(401L,922L,4006L),partial_r2=c(.001,.005,.01,.02))
grid$power <- mapply(power_at,grid$n,grid$partial_r2)
grid$min_r2_for_80pct <- vapply(grid$n,function(n) uniroot(function(r) power_at(n,r)-.8,c(.00001,.2))$root,numeric(1))
write.csv(grid,'outputs/stage6_score_validation_idealized_power.csv',row.names=FALSE)
print(grid)
