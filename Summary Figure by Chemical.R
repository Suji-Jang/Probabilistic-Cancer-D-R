library(viridis)
library(ggplot2)
library(stats)
library(ggpubr)
library(data.table)
library(lsr)
library(scales)

# Data loading and Parameter setting
functionfolder <- "functions"
source(file.path(functionfolder,"BMD_functions.R"))
source(file.path(functionfolder,"Pop_incidence_functions.R"))
source(file.path(functionfolder,"Extra_risk_functions.R"))

bmdfolder <- "BMD-data"
total.data <- fread(file.path(bmdfolder,"Total_data_121821.csv"))
datasets_dich <- fread(file.path(bmdfolder,"BMD_data_121821.csv"))
bw.a.df <- fread(file.path(bmdfolder,"BW_data_final_121821.csv"))
bmds.df <- fread(file.path(bmdfolder,"BMDS_RSD_020422.csv")) 
bbmd.df <- fread(file.path(bmdfolder,"Summary_Output_data.csv"))

resultsfolder <- "results"
rsd.06.df <- fread(file.path(resultsfolder,"RSD06.quants.csv"))
rsd.05.df <- fread(file.path(resultsfolder,"RSD05.quants.csv"))
rsd.04.df <- fread(file.path(resultsfolder,"RSD04.quants.csv"))
hdmi.output.df <- fread(file.path(resultsfolder,"HDMI_filtered.csv"))

study.index <- unique(datasets_dich[,1])
num.study <- length(study.index)
study.index.all <- datasets_dich[[1]]
dose.all <- datasets_dich[[2]]
subjnum.all <- datasets_dich[[3]]
casenum.all <- datasets_dich[[4]]
zeroish <- 1e-8  # Globle Variable

bw.a.df <- bw.a.df[c(1:98,100:197,199:273),c("Index","bw.a")]
bw.a.df$Index <- 1:271
colnames(bw.a.df) <- c("Dataset","bw.a")
bw.h <- 70

bbmd.bmdl <- bbmd.df[,c("V1","MA-ext10L")]
colnames(bbmd.bmdl) <- c("Index","BMDL")
bbmd.bmdl$DAF <- (bw.a.df$bw.a/70)^0.25
bbmd.bmdl$HED <- bbmd.bmdl$BMDL*bbmd.bmdl$DAF

rsd.06.5th <- rsd.06.df$X5.
rsd.05.5th <- rsd.05.df$X5.
rsd.04.5th <- rsd.04.df$X5.

for (s in 1:271) {
  dose.temp<-dose.all[study.index.all==s]
  nsub<-subjnum.all[study.index.all==s]
  ncase<-casenum.all[study.index.all==s]
  dose.max<-max(dose.temp)
  min.resp.bmr <- min(ncase/nsub)+0.1
  
  ma.bmdl <- bbmd.bmdl$BMDL[s]
  bbmd.hed <- bbmd.bmdl$HED[s]
  
  # # Figure 5A - length(study.index)
  temp.df <- read.csv(file.path(bmdfolder,"BBMD-outputs",paste0("Index ",s,".csv")),as.is=TRUE)
  temp.df <- temp.df[,2:52]
  plot.df <- as.data.frame(t(apply(temp.df,2,quantile,probs=c(0.05,0.5,0.95))))
  plot.df$Dose <- rownames(plot.df)
  plot.df$Dose <- as.numeric(gsub("X","",plot.df$Dose))

  act.data <- data.frame(Dose=dose.temp,Response=ncase/nsub)

  bmr.approx <- approx(plot.df$Dose,plot.df$`95%`,xout=ma.bmdl)[2][[1]]

  p1 <- ggplot(plot.df,aes(x=Dose,y=`50%`)) +
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#FDE725FF") +
    geom_line() + ylim(0,1) +
    geom_point(data=act.data,aes(x=Dose,y=Response)) +
    geom_point(aes(x=ma.bmdl,y=bmr.approx),shape=6,size=3) +
    geom_segment(aes(x=-Inf,y=bmr.approx,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    geom_segment(aes(x=ma.bmdl,y=-Inf,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    theme_classic() + xlab("Dose (mg/kg-day)") + ylab("Incidence") + ggtitle("D-R 90% CI")

  # Figure 5B
  rsd.06 <- rsd.06.5th[s]
  rsd.05 <- rsd.05.5th[s]
  rsd.04 <- rsd.04.5th[s]
  
  samp.parms <- subset(total.data,Dataset==s)[,-1] # Removing first column (dataset number) #
  bw.a <- subset(bw.a.df,Dataset==s)$bw.a[1] # Animal body weight
  set.seed(3.14159)
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a,
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  Dose.vec <- 10^seq(floor(min(rsd.06,log10(bbmd.hed*1e-5))),
                     ceiling(log10(dose.max)),0.5)
  NormDose.vec <- Dose.vec/dose.max
  pop.incidence.df <- get.pop.incidence(NormDose.vec,
                                        uncertainty.samples,
                                        n.pop = 1e4,
                                        prob = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975))
  pop.incidence.df$OrigDose <- pop.incidence.df$Dose * dose.max
  
  hdmi.med.plot.df <- subset(hdmi.output.df,Incidence==0.5&Index==s) # Extracting I=50%

  p2 <- ggplot(hdmi.med.plot.df) +
    geom_line(aes(x=`X50.`,y=indivrisk)) +
    geom_ribbon(aes(xmin=`X5.`,xmax=`X95.`,y=indivrisk),fill="#541352FF",alpha=0.5) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    annotation_logticks(  short = unit(0.05, "cm"),
                          mid = unit(0.1, "cm"),
                          long = unit(0.15, "cm"),
                          size = 0.25)+
    theme_classic() +
    coord_cartesian(xlim=range(pop.incidence.df$OrigDose),
                    ylim=c(1e-6,0.1),expand=FALSE)+
    geom_point(aes(x=bbmd.hed,y=0.1),shape=25,size=5,fill="black") + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed, 
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    xlab("Dose (mg/kg-day)") + ylab("Individual Risk (I = 50%)") + ggtitle("Individual Risk(I = 50%)")
  
  # Figure 5C
  hdmi.01.plot.df <- subset(hdmi.output.df,Incidence==0.01&Index==s) # Extracting I=1%

  p3 <- ggplot(hdmi.01.plot.df) +
    geom_line(aes(x=`X50.`,y=indivrisk)) +
    geom_ribbon(aes(xmin=`X5.`,xmax=`X95.`,y=indivrisk),fill="#2A788EFF",alpha=0.5) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    annotation_logticks(  short = unit(0.05, "cm"),
                          mid = unit(0.1, "cm"),
                          long = unit(0.15, "cm"),
                          size = 0.25)+
    theme_classic() +
    coord_cartesian(xlim=range(pop.incidence.df$OrigDose),
                    ylim=c(1e-6,0.1),expand=FALSE)+
    geom_point(aes(x=bbmd.hed,y=0.1),shape=25,size=5,fill="black") + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed,
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    xlab("Dose (mg/kg-day)") + ylab("Individual Risk (I = 1%)") + ggtitle("Individual Risk(I = 1%)")
  
  # Figure 5D
  p4 <- ggplot(pop.incidence.df,aes(x=OrigDose,y=`50%`)) +
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#7AD151FF") + # 90% CI of population incidence
    geom_line() +  # Median population incidence
    geom_point(aes(x=bbmd.hed,y=0.1),shape=25,size=5,fill="black") + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed, 
                     y = 1e-6,yend = 0.1),linetype="dashed") +
    geom_point(aes(x=10^rsd.06,y=1e-6),shape=17,size=5) +
    # geom_point(aes(x=10^rsd.05,y=1e-5),shape=18,size=3) +
    # geom_point(aes(x=10^rsd.04,y=1e-4),shape=18,size=3) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", scales::math_format(10^.x))) + 
    annotation_logticks(  short = unit(0.05, "cm"),
                          mid = unit(0.1, "cm"),
                          long = unit(0.15, "cm"),
                          size = 0.25)+
    theme_classic() +
    coord_cartesian(xlim=range(pop.incidence.df$OrigDose),
                    ylim=c(1e-6,0.1),expand=FALSE)+
    xlab("Dose (mg/kg-day)") + ylab("Population Risk") + ggtitle("Population Risk") +
    theme(legend.position="none")
  
  # ggarrange
  sumplot <- ggarrange(p1,p2,p3,p4,ncol=1) 
  annotate_figure(sumplot,top=text_grob(paste("Index",s), face="bold",size=14))
  ggsave(file.path(resultsfolder,"Summary Figures",paste("SummaryFigure - Index",s,".pdf")),width=4,height=12)
  
}
