library(viridis)
library(ggplot2)
library(stats)
library(ggpubr)
library(lsr)
library(scales)
library(data.table)

functionfolder <- "functions"
source(file.path(functionfolder,"Pop_incidence_functions.R"))
source(file.path(functionfolder,"BMD_functions.R"))
source(file.path(functionfolder,"HDMI_functions.R"))
source(file.path(functionfolder,"Extra_risk_functions.R"))
source(file.path(functionfolder,"RSD_functions.R"))
bmdfolder <- "BMD-data"
total.data <- fread(file.path(bmdfolder,"Total_data_092022.csv"))
bw.a.df <- fread(file.path(bmdfolder,"BW_data_final_092022.csv"))
bbmd.df <- fread(file.path(bmdfolder,"Summary_Output_data_092022.csv"))
datasets_dich <- fread(file.path(bmdfolder,"BMD_data_092022.csv"))
resultsfolder <- "results"
hdmi.df <- fread(file.path(resultsfolder,"HDMI_data_092922.csv"))
rsd.06.df <- fread(file.path(resultsfolder,"RSD06.quants.csv"))
figuresfolder <- "figures"

# Data loading and Parameter setting
study.index <- unique(datasets_dich$Study_Index)
num.study <- length(study.index)
study.index.all <- datasets_dich$Study_Index
dose.all <- datasets_dich$Dose
subjnum.all <- datasets_dich$Number
casenum.all <- datasets_dich$Effect
zeroish <- 1e-8  # Globle Variable

bw.a.df <- bw.a.df[,c("Index","bw.a")]
colnames(bw.a.df) <- c("Dataset","bw.a")
bw.h <- 70

bbmd.bmdl <- bbmd.df[,c("V1","MA-ext10L")]
colnames(bbmd.bmdl) <- c("Index","BMDL")
bbmd.bmdl$DAF <- (bw.a.df$bw.a/70)^0.25
bbmd.bmdl$HED <- bbmd.bmdl$BMDL*bbmd.bmdl$DAF

rsd.06.5th <- rsd.06.df$X5.

for (s in 1:length(study.index)){
  s <- 1
  dose.temp<-dose.all[study.index.all==s]
  nsub<-subjnum.all[study.index.all==s]
  ncase<-casenum.all[study.index.all==s]
  dose.max<-max(dose.temp)
  min.resp.bmr <- min(ncase/nsub)+0.1
  
  ma.bmdl <- bbmd.bmdl$BMDL[s]
  bbmd.hed <- bbmd.bmdl$HED[s]
  
  # Figure 5A - length(study.index)
  temp.df <- fread(file.path(resultsfolder,"D-R Samples",paste0("Index ",s,".csv")))
  
  temp.df <- temp.df[,2:52]
  plot.df <- as.data.frame(t(apply(temp.df,2,quantile,probs=c(0.05,0.5,0.95))))
  plot.df$Dose <- rownames(plot.df)
  plot.df$Dose <- as.numeric(gsub("X","",plot.df$Dose))
  
  act.data <- data.frame(Dose=dose.temp,Response=ncase/nsub)
  
  bmr.approx <- approx(plot.df$Dose,plot.df$`95%`,xout=ma.bmdl)[2][[1]]
  
  p1 <- ggplot(plot.df,aes(x=Dose,y=`50%`)) + labs(tag="A") +
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#FDE725FF") +
    geom_line() + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,by=0.2)) + 
    geom_point(data=act.data,aes(x=Dose,y=Response)) +
    geom_point(aes(x=ma.bmdl,y=bmr.approx),shape=2,size=3) +
    geom_segment(aes(x=-Inf,y=bmr.approx,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    geom_segment(aes(x=ma.bmdl,y=-Inf,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    theme_classic() + ylab("Incidence") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA),
          plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.18, 0.93))
  

  # Figure 5B
  rsd.06 <- rsd.06.5th[s]

  samp.parms <- subset(total.data,Dataset==s)[,-1] # Removing first column (dataset number) #
  bw.a <- subset(bw.a.df,Dataset==s)$bw.a[1] # Animal body weight
  set.seed(3.14159)
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a,
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  Dose.vec <- 10^seq(rsd.06,3,0.5)
  NormDose.vec <- Dose.vec/dose.max
  pop.incidence.df <- get.pop.incidence(NormDose.vec,
                                        uncertainty.samples,
                                        n.pop = 1e4,
                                        prob = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975))
  pop.incidence.df$OrigDose <- pop.incidence.df$Dose * dose.max
  
  hdmi.med.plot.df <- subset(hdmi.df,Incidence==0.5&Index==s)[,c("indivrisk","Incidence","X5.","X50.","X95.")] # Extracting I=50%

  p2 <- ggplot(hdmi.med.plot.df) + labs(tag="B") +
    geom_line(aes(x=`X50.`,y=indivrisk)) +
    geom_ribbon(aes(xmin=`X5.`,xmax=`X95.`,y=indivrisk),fill="#541352FF",alpha=0.5) +
    scale_x_log10(breaks=10^(-7:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + 
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed, # Old x=bbmd.hed
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    ylab("Individual Risk (I = 50%)") + annotation_logticks(sides="b") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA),
          plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.18, 0.93))
  
  
  # Figure 5C
  hdmi.01.plot.df <- subset(hdmi.df,Incidence==0.01&Index==s)[,c("indivrisk","Incidence","X5.","X50.","X95.")] # Extracting I=1%

  p3 <- ggplot(hdmi.01.plot.df) + labs(tag="C") +
    geom_line(aes(x=`X50.`,y=indivrisk)) +
    geom_ribbon(aes(xmin=`X5.`,xmax=`X95.`,y=indivrisk),fill="#2A788EFF",alpha=0.5) +
    scale_x_log10(breaks=10^(-7:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed,
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    ylab("Individual Risk (I = 1%)") + annotation_logticks(sides="b") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA),
          plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.18, 0.93))
  
  
  # Figure 5D
  p4 <- ggplot(pop.incidence.df,aes(x=OrigDose,y=`50%`)) + labs(tag="D") +
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#7AD151FF") + # 90% CI of population incidence
    geom_line() +  # Median population incidence
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed,
                     y = 1e-6,yend = 0.1),linetype="dashed") +
    geom_point(aes(x=10^rsd.06,y=1e-6),shape=18,size=3) +
    scale_x_log10(breaks=10^(-7:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    xlab("Dose (mg/kg-d)") + ylab("Population Risk") + annotation_logticks(sides="b") +
    theme(legend.position="none",panel.border = element_rect(colour = "black",fill=NA),
          plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.18, 0.93))
  
  # ggarrange
  plot <- ggarrange(p1,p2,p3,p4,ncol=1)
  annotate_figure(plot,top=text_grob(paste("Index",s), face="bold",size=14))
  ggsave(file.path(figuresfolder,"Figure 5 - Summary Figrues",paste("Fig 5 - Index",s,"_103122.pdf")),width=4,height=12)
  
}
