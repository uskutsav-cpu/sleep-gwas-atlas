folder<-'joint_fit';y<-scan(file.path(folder,'moments.txt'),quiet=TRUE);d<-as.matrix(read.table(file.path(folder,'deletes.txt')))
index<-function(z)c(sum(z[4:6]),sum(z[7:9]))/sqrt(3+2*sum(z[1:3]))
idx<-index(y);id<-t(apply(d,1,index));se<-function(a)sqrt(199/200*sum((a-mean(a))^2))
draw<-function(){
 old<-par(no.readonly=TRUE);on.exit(par(old));par(mfrow=c(1,2),mar=c(4.5,10,3,1),oma=c(3.5,0,2.5,0),family='sans',cex=1)
 labels<-c('Walking pace','Left grip strength','Moderate activity','Equal-weight genetic index')
 for(s in 1:2){
  offset<-if(s==1)3 else 6;point<-c(y[(offset+1):(offset+3)],idx[s]);err<-c(apply(d[,(offset+1):(offset+3)],2,se),se(id[,s]));pos<-4:1
  color<-if(s==1)'#147D92' else '#7546A0'
  plot(NA,xlim=c(-.4,.15),ylim=c(.5,4.5),xlab='Genetic correlation',ylab='',yaxt='n',bty='n',main=if(s==1)'Registry insomnia' else 'Questionnaire insomnia',cex.main=1.1)
  abline(v=0,col='#777777',lty=2);abline(h=pos,col='#ECECEC')
  axis(2,at=pos,labels=labels,las=1,tick=FALSE,cex.axis=.85)
  segments(point-1.95996398454*err,pos,point+1.95996398454*err,pos,col=color,lwd=2)
  segments(point-1.95996398454*err,pos-.07,point-1.95996398454*err,pos+.07,col=color)
  segments(point+1.95996398454*err,pos-.07,point+1.95996398454*err,pos+.07,col=color)
  points(point,pos,pch=c(16,16,16,18),cex=c(1.15,1.15,1.15,1.6),col=color)
 }
 mtext('Insomnia and physical-function genetics',outer=TRUE,side=3,line=.8,cex=1.1,font=2)
 mtext('Intervals: unadjusted normal 95%, using 200-block jackknife SEs. Same 1,053,614 SNPs.',outer=TRUE,side=1,line=.8,cex=.75)
 mtext('Post-hoc physical genetic index; shared UKB data and source/reference assumptions limit interpretation.',outer=TRUE,side=1,line=2,cex=.75)
}
png('explicit_index.png',width=1800,height=760,res=160);draw();dev.off()
svg('explicit_index.svg',width=11.25,height=4.75);draw();dev.off()
