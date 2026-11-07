d<-read.csv('outputs/shared-genetics-atlas/phenotype_comparison/ukb_cross_tab.csv')
d$fried_category<-factor(d$fried_category,levels=c('robust','prefrail','frail'));d$fi_category<-factor(d$fi_category,levels=c('robust','mild','moderate_severe'));m<-xtabs(n~fried_category+fi_category,d);q<-m/sum(m);e<-outer(rowSums(q),colSums(q));dist<-abs(outer(0:2,0:2,'-'));weights<-list(unweighted=1*(dist==0),linear=1-dist/2,quadratic=1-(dist/2)^2)
x<-read.delim('outputs/shared-genetics-atlas/phenotype_comparison/category_agreement/agreement.tsv')
for(n in names(weights)){w<-weights[[n]];v<-(sum(w*q)-sum(w*e))/(1-sum(w*e));stopifnot(abs(v-x$kappa[x$weight==n])<1e-12)}
cat('PASS: independent R computation matches all three kappas\n')
