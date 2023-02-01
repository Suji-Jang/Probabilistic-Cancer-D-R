library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(tidyr)
figuresfolder <- "figures"
resultsfolder <- "results"

pdf(file.path(figuresfolder,"Supple Figure - RSD Eta-squared.pdf"))
for (k in 1:3) {
  if (k==1) {
    risklevel <- 1e-4
    RSD.label <- "RSD04"
  } else if (k==2) {
    risklevel <- 1e-5
    RSD.label <- "RSD05"
  } else if (k==3) {
    risklevel <- 1e-6
    RSD.label <- "RSD06"
  }
  RSD.etasq.df<-fread(file.path(resultsfolder,paste0(RSD.label,".etasq.csv")))
  names(RSD.etasq.df)[7]<-"BMD Model"
  names(RSD.etasq.df)[11]<-"BMD Parameters"
  print(ggpairs(RSD.etasq.df[,c(5,7:11)])+ggtitle(RSD.label)+theme_bw())
  print(ggplot(RSD.etasq.df)+geom_histogram(aes(RSD.log10ci90))+ggtitle(RSD.label)+theme_classic())
  RSD.etasq.df.long <- pivot_longer(RSD.etasq.df[,c(1,7:11)],cols=2:6)
  print(ggplot(RSD.etasq.df.long)+geom_histogram(aes(value))+facet_wrap(~name)+theme_bw())
  RSD.etasq.df$eta_max <- aggregate(value~Dataset,data=RSD.etasq.df.long,max)[[2]]
  RSD.etasq.df$var_max.indx <- aggregate(value~Dataset,data=RSD.etasq.df.long,which.max)[[2]]
  RSD.etasq.df$var_max <- (names(RSD.etasq.df)[7:11])[RSD.etasq.df$var_max.indx]
  RSD.etasq.df$var_max <- factor(RSD.etasq.df$var_max,
                                 levels=c("BMD Model","BMD Parameters","DAF","AHU","sigmaH"))
  RSD.etasq.df$RSD.log10ci90.bin <- cut(RSD.etasq.df$RSD.log10ci90,breaks=seq(1,8))
  print(ggplot(RSD.etasq.df,aes(x="All Chemicals",fill=var_max))+geom_bar(position="fill")+
          xlab("")+ylab("Proportion of chemicals")+theme_classic()+
          scale_fill_viridis_d("Largest contribution\nto uncertainty",
                               direction=-1,begin=0.2,drop=FALSE)+ggtitle(RSD.label))
  print(ggplot(RSD.etasq.df,aes(x=RSD.log10ci90.bin,fill=var_max))+geom_bar(position="fill")+
    xlab("90% Confidence Width (log10 units)")+ylab("Proportion of chemicals")+theme_classic()+
    scale_fill_viridis_d("Largest contribution\nto uncertainty",
                         direction=-1,begin=0.2,drop=FALSE)+ggtitle(RSD.label))
  print(ggplot(RSD.etasq.df,aes(x=RSD.log10ci90.bin,fill=var_max))+geom_bar()+
          xlab("90% Uncertainty Width (log10 units)")+ylab("Number of chemicals")+theme_classic()+
          scale_fill_viridis_d("Largest contribution\nto uncertainty",
                               direction=-1,begin=0.2,drop=FALSE)+ggtitle(RSD.label))
}
dev.off()
