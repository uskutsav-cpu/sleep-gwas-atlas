# Independent base-R matrix replay and finite differences on actual data.
out <- '/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/outputs/sleep_frailty_ml_representation'
readmat <- function(name) as.matrix(read.csv(file.path(out,name),header=FALSE))
x<-readmat('audit_magnitude_design.csv');parameters<-readmat('audit_magnitude_parameters.csv')
labels<-readmat('audit_actual_labels.csv');beta<-parameters[,1];penalty<-parameters[,2]
objective<-function(b) {
  residual<-labels[,1]+drop(x%*%b)-labels[,2]
  sum(labels[,3]*sqrt(residual^2+.01^2))/sum(labels[,3])+.5*sum(penalty*b^2)
}
residual<-labels[,1]+drop(x%*%beta)-labels[,2]
gradient<-drop(crossprod(x,labels[,3]*residual/sqrt(residual^2+.01^2)))/sum(labels[,3])+penalty*beta
error<-max(abs(gradient-parameters[,3]));stopifnot(error<1e-10)
selected<-unique(c(order(abs(gradient),decreasing=TRUE)[1:8],which(penalty==0)[1:2]))
relative<-numeric()
for (index in selected) {
  bplus<-beta;bminus<-beta;bplus[index]<-bplus[index]+1e-6;bminus[index]<-bminus[index]-1e-6
  finite<-(objective(bplus)-objective(bminus))/2e-6
  relative<-c(relative,abs(finite-gradient[index])/max(abs(finite),abs(gradient[index]),1e-10))
}
stopifnot(max(relative)<1e-4)
ridge<-readmat('audit_signed_ridge_system.csv');p<-nrow(ridge)
independent_beta<-solve(ridge[,seq_len(p),drop=FALSE],ridge[,p+1])
ridge_error<-max(abs(independent_beta-ridge[,p+2]));stopifnot(ridge_error<1e-10)
a<-labels[,3]/sum(labels[,3]);u<-residual/sqrt(residual^2+.01^2)
unpen<-which(penalty==0);pen<-which(penalty>0)
stopifnot(all(rowSums(x[,unpen,drop=FALSE])==1))
line<-max.col(x[,unpen,drop=FALSE],ties.method='first')
for (value in unique(line)) {
  ix<-line==value;u[ix]<-u[ix]-sum(a[ix]*u[ix])/sum(a[ix])
  u[ix]<-u[ix]/max(1,max(abs(u[ix])))
}
dual_constraint_error<-max(abs(drop(crossprod(x[,unpen,drop=FALSE],a*u))))
stopifnot(dual_constraint_error<1e-12,max(abs(u))<=1+1e-12)
dual_gradient<-drop(crossprod(x[,pen,drop=FALSE],a*u))
dual<-sum(a*(u*(labels[,1]-labels[,2])+.01*sqrt(pmax(1-u*u,0))))-sum(dual_gradient^2)/.2
stopifnot(objective(beta)>=dual-1e-12)
write.table(data.frame(passed=TRUE,actual_rows=nrow(x),magnitude_batch_objective=objective(beta),
                       maximum_gradient_disagreement=error,max_finite_difference_relative_error=max(relative),
                       signed_ridge_independent_solve_error=ridge_error,
                       dual_line_constraint_error=dual_constraint_error,
                       independently_computed_batch_dual_bound=dual,
                       independently_computed_batch_primal_dual_gap=objective(beta)-dual),
            file.path(out,'independent_R_readout_numerical_receipt.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('READOUT NUMERICS INDEPENDENT R CHECK PASS\n')
