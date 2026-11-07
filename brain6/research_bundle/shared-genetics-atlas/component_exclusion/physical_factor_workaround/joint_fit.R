args <- commandArgs(trailingOnly=TRUE)
folder <- if (length(args)) args[1] else 'joint_fit'
y <- scan(file.path(folder,'moments.txt'),quiet=TRUE)
V <- as.matrix(read.table(file.path(folder,'covariance.txt')))
deletes <- as.matrix(read.table(file.path(folder,'deletes.txt')))
stopifnot(length(y)==9,nrow(deletes)==200,ncol(deletes)==9,min(eigen(V,symmetric=TRUE)$values)>0)
W <- solve(V)
prediction <- function(b) c(b[1]*b[2],b[1]*b[3],b[2]*b[3],b[1]*b[4],b[2]*b[4],b[3]*b[4],b[1]*b[5],b[2]*b[5],b[3]*b[5])
jacobian <- function(b) rbind(c(b[2],b[1],0,0,0),c(b[3],0,b[1],0,0),c(0,b[3],b[2],0,0),c(b[4],0,0,b[1],0),c(0,b[4],0,b[2],0),c(0,0,b[4],b[3],0),c(b[5],0,0,0,b[1]),c(0,b[5],0,0,b[2]),c(0,0,b[5],0,b[3]))
objective <- function(b,z) { e<-z-prediction(b);as.numeric(t(e)%*%W%*%e) }
gradient <- function(b,z) as.vector(-2*t(jacobian(b))%*%W%*%(z-prediction(b)))
lower<-c(0,0,0,-1,-1);upper<-rep(1,5)
starts<-rbind(c(.67,.44,.19,-.11,-.41),c(.3,.6,.5,-.3,-.3),c(.8,.3,.1,0,-.7),c(.5,.5,.5,.2,.2),c(.9,.9,.9,-.7,-.7))
fit <- function(z,initial) optim(initial,objective,gradient,z=z,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=2000,factr=1e5,pgtol=1e-8))
fits<-lapply(seq_len(nrow(starts)),function(i)fit(y,starts[i,])); stopifnot(all(vapply(fits,function(x)x$convergence==0,logical(1))))
best<-fits[[which.min(vapply(fits,function(x)x$value,numeric(1)))]]
checks<-lapply(seq_len(nrow(starts)),function(i)nlminb(starts[i,],objective,gradient,z=y,lower=lower,upper=upper,control=list(iter.max=2000,eval.max=5000,rel.tol=1e-12)))
stopifnot(all(vapply(checks,function(x)x$convergence==0,logical(1))))
check<-checks[[which.min(vapply(checks,function(x)x$objective,numeric(1)))]]
stopifnot(abs(check$objective-best$value)<1e-6)
start_results<-do.call(rbind,lapply(seq_len(nrow(starts)),function(i)rbind(data.frame(optimizer='LBFGSB',start=i,objective=fits[[i]]$value,t(fits[[i]]$par)),data.frame(optimizer='nlminb',start=i,objective=checks[[i]]$objective,t(checks[[i]]$par)))))
write.table(start_results,file.path(folder,'multistart.txt'),row.names=FALSE,quote=FALSE)
write.table(data.frame(observed=y,predicted=prediction(best$par),residual=y-prediction(best$par)),file.path(folder,'moment_residuals.txt'),row.names=FALSE,quote=FALSE)
# Check the analytic gradient against symmetric finite differences at an interior point.
btest<-c(.61,.47,.23,-.14,-.38); step<-1e-6
numeric_gradient<-sapply(seq_len(5),function(i){d<-rep(0,5);d[i]<-step;(objective(btest+d,y)-objective(btest-d,y))/(2*step)})
stopifnot(max(abs(numeric_gradient-gradient(btest,y)))<1e-5)
block<-matrix(NA_real_,200,5);convergence<-numeric(200)
for (i in seq_len(200)) {
 a<-fit(deletes[i,],best$par);b<-fit(deletes[i,],starts[2,]);
 candidates<-list(a,b);ok<-vapply(candidates,function(x)x$convergence==0,logical(1));stopifnot(any(ok));valid<-candidates[ok]
 chosen<-valid[[which.min(vapply(valid,function(x)x$value,numeric(1)))]];block[i,]<-chosen$par;convergence[i]<-chosen$convergence
}
covariance<-199/200*crossprod(sweep(block,2,colMeans(block)))
contrast<-block[,5]-block[,4];contrast_se<-sqrt(199/200*sum((contrast-mean(contrast))^2));delta<-best$par[5]-best$par[4]
labels<-c('walking_loading','grip_loading','activity_loading','registry_factor_correlation','questionnaire_factor_correlation')
write.table(data.frame(parameter=labels,estimate=best$par,jackknife_se=sqrt(diag(covariance))),file.path(folder,'parameters.txt'),row.names=FALSE,quote=FALSE)
write.table(data.frame(block=seq_len(200),block),file.path(folder,'parameter_deletes.txt'),row.names=FALSE,quote=FALSE)
write.table(data.frame(metric=c('gls_Q','approximate_df','approximate_fit_p','questionnaire_minus_registry','paired_jackknife_se','unadjusted_normal_contrast_p','boundary_full','boundary_delete_count','optimizer_objective_difference','optimizer_parameter_max_difference','multistart_objective_range','gradient_check_max_difference'),value=c(best$value,4,pchisq(best$value,4,lower.tail=FALSE),delta,contrast_se,2*pnorm(-abs(delta/contrast_se)),any(best$par<lower+1e-6 | best$par>upper-1e-6),sum(apply(block,1,function(b)any(b<lower+1e-6|b>upper-1e-6))),abs(check$objective-best$value),max(abs(check$par-best$par)),diff(range(vapply(fits,function(x)x$value,numeric(1)))),max(abs(numeric_gradient-gradient(btest,y))))),file.path(folder,'diagnostics.txt'),row.names=FALSE,quote=FALSE)
print(read.table(file.path(folder,'parameters.txt'),header=TRUE));print(read.table(file.path(folder,'diagnostics.txt'),header=TRUE))
