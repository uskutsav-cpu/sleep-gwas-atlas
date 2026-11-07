p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/excitatory_coloc/finngen_sensitivity/statistic_sensitivity/qtl_t_distribution'
f <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/excitatory_coloc/finngen_sensitivity/matched.tsv'
d <- read.delim(f);stopifnot(nrow(d)==2441,all(d$qtl_se>0),all(d$qtl_p>0),all(d$qtl_p<=1))
t <- abs(d$qtl_beta/d$qtl_se); lead <- which.max(t); logp <- log(d$qtl_p)
objective <- function(df) log(2)+pt(-t[lead],df=df,log.p=TRUE)-logp[lead]
fit <- uniroot(objective,c(20,1000000),tol=1e-9)$root;df <- round(fit)
others <- seq_len(nrow(d))!=lead
tdiff <- abs((log(2)+pt(-t,df=df,log.p=TRUE)-logp)/log(10))
ndiff <- abs((log(2)+pnorm(-t,log.p=TRUE)-logp)/log(10))
q <- function(x) c(median=median(x),p95=unname(quantile(x,.95)),max=max(x))
summ <- rbind(data.frame(method='student_t',as.list(q(tdiff[others]))),data.frame(method='normal',as.list(q(ndiff[others]))))
write.table(summ,file.path(p,'discrepancies.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(n=nrow(d),held_out=sum(others),fitting_variant=d$variant[lead],t_stat=t[lead],source_p=d$qtl_p[lead],df_fitted=fit,df_rounded=df),file.path(p,'fit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(summ);cat('fit:',d$variant[lead],fit,'rounded',df,'\n')
