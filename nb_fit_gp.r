#!/software/biosoft/software/R-3.2.0/bin/Rscript
options(warn = -1)
Args<- commandArgs(trailing = TRUE)
library(MASS)
cov_file <- Args[1]
ploidy <- as.numeric(Args[2])

sample.cov<-read.delim(cov_file,header=F)
colnames(sample.cov)<-c("chr","dep","count","genomeLen","perc")

chr_s="genome"
chr.cov<-subset(sample.cov,chr==chr_s)
median_all<-median(rep(chr.cov$dep,chr.cov$count))
peak_select<-subset(chr.cov,dep>1/2*median_all&dep<3/2*median_all)
peak_fit<-fitdistr(rep(peak_select$dep,peak_select$count),densfun="negative binomial")
likelihood_matrix<-matrix(0,nrow=nrow(chr.cov),ncol=(2*ploidy+1))
for(i in 1:(2*ploidy+1)){
  if(i==1){
    likelihood_matrix[,i]<-dnbinom(chr.cov$dep,size=as.numeric(peak_fit$estimate[1]),mu=1)
  }else{
    likelihood_matrix[,i]<-dnbinom(chr.cov$dep,size=as.numeric(peak_fit$estimate[1]),mu=(i-1)/ploidy*as.numeric(peak_fit$estimate[2]))
  }
}
bootstrap_abs=1
i=1
prior_old<-data.frame(Var1=1:(2*ploidy+1),Freq_p=rep(1/(2*ploidy+1),(2*ploidy+1)))
prior_matrix<-data.frame(row.names = 1:(2*ploidy+1))

while( bootstrap_abs > 0){
chr.cov$gt<-apply(likelihood_matrix,1,function(a){b=a*prior_old$Freq_p;which.max(b)})
chr.cov$gt<-factor(chr.cov$gt,levels = 1:(2*ploidy+1))
prior_new<-data.frame(Freq=tapply(chr.cov$count,chr.cov$gt,sum))
prior_new[is.na(prior_new$Freq),]=0
prior_new$Freq_p<-prior_new$Freq/sum(prior_new$Freq)
if(ploidy>1){
  prior_new[2:ploidy,"Freq_p"]<-sum(prior_new[2:ploidy,"Freq_p"])/(ploidy-1)
  prior_new[(ploidy+2):(2*ploidy),"Freq_p"]<-sum(prior_new[(ploidy+2):(2*ploidy),"Freq_p"])/(ploidy-1)
}

bootstrap_abs<-sum(abs(prior_new$Freq_p-prior_old$Freq_p))
prior_old=prior_new
prior_matrix[,i]<-prior_new$Freq_p
i=i+1
}
colnames(prior_matrix)<-1:ncol(prior_matrix)
result_cov<-data.frame()
for (chr_s in levels(sample.cov$chr)){
if(chr_s!="genome"){
chr.cov<-subset(sample.cov,chr==chr_s)
likelihood_matrix<-matrix(0,nrow=nrow(chr.cov),ncol=(2*ploidy+1))
for(i in 1:(2*ploidy+1)){
  if(i==1){
    likelihood_matrix[,i]<-dnbinom(chr.cov$dep,size=as.numeric(peak_fit$estimate[1]),mu=1)
  }else{
    likelihood_matrix[,i]<-dnbinom(chr.cov$dep,size=as.numeric(peak_fit$estimate[1]),mu=(i-1)/ploidy*as.numeric(peak_fit$estimate[2]))
  }
}
chr.cov$gt<-apply(likelihood_matrix,1,function(a){b=a*prior_new$Freq_p;which.max(b)})
chr.cov$gt<-as.numeric(as.character(chr.cov$gt))
chr.cov$gt<-chr.cov$gt-1
chr.cov$gq<-apply(likelihood_matrix,1,function(a){b=a*prior_new$Freq_p;b=b[order(b,decreasing=T)];(log10(b[1])-log10(b[2]))*10})
result_cov<-rbind(result_cov,chr.cov)
}
}
write.table(result_cov,row.names = F,col.names = F,quote=F,sep="\t")
