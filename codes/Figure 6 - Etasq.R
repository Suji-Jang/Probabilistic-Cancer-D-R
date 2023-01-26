library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(tidyr)
figuresfolder <- "figures"
resultsfolder <- "results"

risklevel <- 1e-6
RSD.label <- "RSD06"
RSD.etasq.df<-fread(file.path(resultsfolder,paste0(RSD.label,".etasq.csv")))
names(RSD.etasq.df)[7]<-"BMD Model"
names(RSD.etasq.df)[11]<-"BMD Parameters"
RSD.etasq.df.long <- pivot_longer(RSD.etasq.df[,c(1,7:11)],cols=2:6)
RSD.etasq.df$eta_max <- aggregate(value~Dataset,data=RSD.etasq.df.long,max)[[2]]
RSD.etasq.df$var_max.indx <- aggregate(value~Dataset,data=RSD.etasq.df.long,which.max)[[2]]
RSD.etasq.df$var_max <- (names(RSD.etasq.df)[7:11])[RSD.etasq.df$var_max.indx]
RSD.etasq.df$var_max <- factor(RSD.etasq.df$var_max,
                               levels=c("BMD Model","BMD Parameters","DAF","AHU","sigmaH"))
RSD.etasq.df$RSD.log10ci90.bin <- cut(RSD.etasq.df$RSD.log10ci90,breaks=seq(1,8))

pdf(file.path(figuresfolder,"Eta-squared - Resized(6,4).pdf"),width=6,height=4)
print(ggplot(RSD.etasq.df,aes(x=RSD.log10ci90.bin,fill=var_max))+geom_bar()+
        xlab("90% Uncertainty Width (log10 units)")+ylab("Number of chemicals")+theme_classic()+
        scale_fill_viridis_d("Largest contribution\nto uncertainty",
                             direction=-1,begin=0.2,drop=FALSE))
dev.off()
