folder <- 'joint_fit'
y<-scan(file.path(folder,'moments.txt'),quiet=TRUE)
V<-as.matrix(read.table(file.path(folder,'covariance.txt')));W<-solve(V)
deletes<-as.matrix(read.table(file.path(folder,'deletes.txt')))
original<-read.table(file.path(folder,'parameters.txt'),header=TRUE)$estimate
lower<-c(rep(1e-6,3),rep(-1+1e-6,2),-1,-1)
upper<-c(rep(1-1e-6,5),1,1)
results<-list();diagnostics<-list()
for (indicator in seq_len(3)) {
 label<-c('walking','grip','activity')[indicator]
 prediction<-function(b) {
  m<-c(b[1]*b[2],b[1]*b[3],b[2]*b[3],b[1]*b[4],b[2]*b[4],b[3]*b[4],b[1]*b[5],b[2]*b[5],b[3]*b[5])
  m[3+indicator]<-m[3+indicator]+sqrt((1-b[indicator]^2)*(1-b[4]^2))*b[6]
  m[6+indicator]<-m[6+indicator]+sqrt((1-b[indicator]^2)*(1-b[5]^2))*b[7]
  m
 }
 objective<-function(b,z){e<-z-prediction(b);as.numeric(t(e)%*%W%*%e)}
 gradient<-function(b,z){
  J<-cbind(rbind(c(b[2],b[1],0,0,0),c(b[3],0,b[1],0,0),c(0,b[3],b[2],0,0),c(b[4],0,0,b[1],0),c(0,b[4],0,b[2],0),c(0,0,b[4],b[3],0),c(b[5],0,0,0,b[1]),c(0,b[5],0,0,b[2]),c(0,0,b[5],0,b[3])),matrix(0,9,2))
  for(k in seq_len(2)){
   row<-3*k+indicator;rho<-b[3+k];psi<-b[5+k];loading<-b[indicator]
   J[row,indicator]<-J[row,indicator]-loading*sqrt((1-rho^2)/(1-loading^2))*psi
   J[row,3+k]<-J[row,3+k]-rho*sqrt((1-loading^2)/(1-rho^2))*psi
   J[row,5+k]<-sqrt((1-loading^2)*(1-rho^2))
  }
  as.vector(-2*t(J)%*%W%*%(z-prediction(b)))
 }
 starts<-rbind(c(original,0,0),c(original,.3,.3),c(original,-.3,-.3),c(original,.3,-.3),c(original,-.3,.3))
 btest<-c(original,.1,-.1);step<-1e-6
 numeric_gradient<-sapply(seq_len(7),function(i){d<-rep(0,7);d[i]<-step;(objective(btest+d,y)-objective(btest-d,y))/(2*step)})
 stopifnot(max(abs(numeric_gradient-gradient(btest,y)))<1e-4)
 fit<-function(z,initial)nlminb(initial,objective,gradient,z=z,lower=lower,upper=upper,control=list(iter.max=3000,eval.max=10000,rel.tol=1e-10))
 fits<-lapply(seq_len(nrow(starts)),function(i)fit(y,starts[i,]));ok<-vapply(fits,function(x)x$convergence==0,logical(1));stopifnot(any(ok))
 valid<-fits[ok];best<-valid[[which.min(vapply(valid,function(x)x$objective,numeric(1)))]]
 alternative<-optim(best$par,objective,gradient,z=y,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=3000,factr=1e5,pgtol=1e-8))
 stopifnot(alternative$convergence==0,abs(alternative$value-best$objective)<1e-5)
 block<-matrix(NA_real_,200,7);failed<-integer(0)
 for(i in seq_len(200)) {
  attempts<-list(fit(deletes[i,],best$par),fit(deletes[i,],starts[1,]))
  good<-vapply(attempts,function(x)x$convergence==0,logical(1))
  if(!any(good)){failed<-c(failed,i);next}
  eligible<-attempts[good];chosen<-eligible[[which.min(vapply(eligible,function(x)x$objective,numeric(1)))]];block[i,]<-chosen$par
 }
 ses<-rep(NA_real_,7);dse<-NA_real_;boundaries<-NA_real_
 if(!length(failed)) {
  ses<-sqrt(diag(199/200*crossprod(sweep(block,2,colMeans(block)))))
  d<-block[,5]-block[,4];dse<-sqrt(199/200*sum((d-mean(d))^2))
  boundaries<-sum(apply(block,1,function(b)any(b<lower+1e-5|b>upper-1e-5)))
 }
 results[[indicator]]<-data.frame(scenario=label,parameter=c('walking_loading','grip_loading','activity_loading','registry_factor_correlation','questionnaire_factor_correlation','registry_residual_correlation','questionnaire_residual_correlation'),estimate=best$par,jackknife_se=ses)
 diagnostics[[indicator]]<-data.frame(scenario=label,gls_Q=best$objective,approximate_df=2,approximate_fit_p=pchisq(best$objective,2,lower.tail=FALSE),questionnaire_minus_registry=best$par[5]-best$par[4],contrast_se=dse,full_boundary=any(best$par<lower+1e-5|best$par>upper-1e-5),delete_boundaries=boundaries,delete_failures=length(failed),converged_starts=sum(ok),start_objective_range=diff(range(vapply(valid,function(x)x$objective,numeric(1)))),optimizer_difference=abs(alternative$value-best$objective))
 write.table(data.frame(block=seq_len(200),block),file.path(folder,paste0('residual_',label,'_deletes.txt')),row.names=FALSE,quote=FALSE)
 write.table(data.frame(observed=y,predicted=prediction(best$par),residual=y-prediction(best$par)),file.path(folder,paste0('residual_',label,'_moments.txt')),row.names=FALSE,quote=FALSE)
}
write.table(do.call(rbind,results),file.path(folder,'residual_parameters.txt'),row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,diagnostics),file.path(folder,'residual_diagnostics.txt'),row.names=FALSE,quote=FALSE)
print(do.call(rbind,results));print(do.call(rbind,diagnostics))
