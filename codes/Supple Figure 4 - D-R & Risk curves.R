  library(viridis)
  library(ggplot2)
  library(stats)
  library(ggpubr)
  library(lsr)
  library(scales)
  library(data.table)
  
  # Load data and functions
  functionfolder <- "functions"
  source(file.path(functionfolder,"Pop_incidence_functions.R"))
  source(file.path(functionfolder,"BMD_functions.R"))
  source(file.path(functionfolder,"HDMI_functions.R"))
  source(file.path(functionfolder,"Extra_risk_functions.R"))
  source(file.path(functionfolder,"RSD_functions.R"))
  bmdfolder <- "BMD-data"
  total.data <- fread(file.path(bmdfolder,"Sample_parameters.csv"))
  bw.a.df <- fread(file.path(bmdfolder,"BW_data.csv"))
  bbmd.df <- fread(file.path(bmdfolder,"Summary_Output_data.csv"))
  datasets_dich <- fread(file.path(bmdfolder,"BMD_data.csv"))
  chemname.df <- fread(file.path(bmdfolder,"Supple Table 1 Cancer Dose Response Data.csv"))
  resultsfolder <- "results"
  hdmi.df <- fread(file.path(resultsfolder,"Supple Table 5 BBMD MA HDMI results.csv"))
  rsd.06.df <- fread(file.path(resultsfolder,"RSD06_quants.csv"))
  figuresfolder <- "figures"
  
  # Set Parameters
  study.index <- unique(datasets_dich$Study_Index)
  num.study <- length(study.index)
  study.index.all <- datasets_dich$Study_Index
  dose.all <- datasets_dich$Dose
  subjnum.all <- datasets_dich$Number
  casenum.all <- datasets_dich$Effect
  zeroish <- 1e-8  # Globle Variable
  
  chem.name <- unique(chemname.df[,c("Index","Chemical.Name")])
  
  bw.a.df <- bw.a.df[,c("Index","bw.a")]
  colnames(bw.a.df) <- c("Dataset","bw.a")
  bw.h <- 70
  
  bbmd.bmdl <- bbmd.df[,c("V1","MA-ext10L")]
  colnames(bbmd.bmdl) <- c("Index","BMDL")
  bbmd.bmdl$DAF <- (bw.a.df$bw.a/70)^0.25
  bbmd.bmdl$HED <- bbmd.bmdl$BMDL*bbmd.bmdl$DAF
  
  rsd.06.5th <- rsd.06.df$X5.
  
# Plot
for (s in 1:length(study.index)){
  # Extract data for each dataset
  dose.temp<-dose.all[study.index.all==s]
  nsub<-subjnum.all[study.index.all==s]
  ncase<-casenum.all[study.index.all==s]
  dose.max<-max(dose.temp)
  min.resp.bmr <- min(ncase/nsub)+0.1
  
  ma.bmdl <- bbmd.bmdl$BMDL[s]
  bbmd.hed <- bbmd.bmdl$HED[s]
  
  # Load D-R data
  drsampfolder <- "results/D-R-Samples"
  temp.df <- fread(file.path(drsampfolder,paste0("Index ",s,".csv")))
  
  temp.df <- temp.df[,2:52]
  plot.df <- as.data.frame(t(apply(temp.df,2,quantile,probs=c(0.05,0.5,0.95))))
  plot.df$Dose <- rownames(plot.df)
  plot.df$Dose <- as.numeric(gsub("X","",plot.df$Dose))
  
  act.data <- data.frame(Dose=dose.temp,Response=ncase/nsub)
  
  # Approximate BMDL
  bmr.approx <- approx(plot.df$Dose,plot.df$`95%`,xout=ma.bmdl)[2][[1]]
  
  # Figure 5 - Example Panel 1 - Dose-response curve (90% CI)
  p1 <- ggplot(plot.df,aes(x=Dose,y=`50%`)) + 
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#FDE725FF") +
    geom_line() + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,by=0.2)) + 
    geom_point(data=act.data,aes(x=Dose,y=Response)) +
    geom_point(aes(x=ma.bmdl,y=bmr.approx),shape=2,size=3) +
    geom_segment(aes(x=-Inf,y=bmr.approx,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    geom_segment(aes(x=ma.bmdl,y=-Inf,xend=ma.bmdl,yend=bmr.approx),linetype="dotted") +
    theme_classic() + ylab("Incidence") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA))
  
  # Population incidence with RSD
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
  
  # Extract HDMI at I = 50%
  hdmi.med.plot.df <- subset(hdmi.df,I_Incidence==0.5&Index==s)[,c("M_indivrisk","I_Incidence",
                                                                   "HDMI.p05","HDMI.p50","HDMI.p95")]
  
  # Figure 5 - Example Panel 2 - Dose-Individual risk of median population (90% CI)
  p2 <- ggplot(hdmi.med.plot.df) + 
    geom_line(aes(x=HDMI.p50,y=M_indivrisk)) +
    geom_ribbon(aes(xmin=HDMI.p05,xmax=HDMI.p95,y=M_indivrisk),fill="#541352FF",alpha=0.5) +
    scale_x_log10(breaks=10^(-10:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + 
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed, # Old x=bbmd.hed
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    ylab("Individual Risk (I = 50%)") + annotation_logticks(sides="b") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA))
  
  # Extract HDMI at I = 1%
  hdmi.01.plot.df <- subset(hdmi.df,I_Incidence==0.01&Index==s)[,c("M_indivrisk","I_Incidence",
                                                                   "HDMI.p05","HDMI.p50","HDMI.p95")]
  
  # Figure 5 - Example Panel 3 - Dose-Individual risk of 1% sensitive population (90% CI)
  p3 <- ggplot(hdmi.01.plot.df) + 
    geom_line(aes(x=HDMI.p50,y=M_indivrisk)) +
    geom_ribbon(aes(xmin=HDMI.p05,xmax=HDMI.p95,y=M_indivrisk),fill="#2A788EFF",alpha=0.5) +
    scale_x_log10(breaks=10^(-10:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed,
                     y = 1e-6,yend = 0.1),linetype="dashed") + 
    ylab("Individual Risk (I = 1%)") + annotation_logticks(sides="b") +
    theme(axis.title.x=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA))

    # Figure 5 - Example Panel 4 - Dose-Population risk (90% CI)
  p4 <- ggplot(pop.incidence.df,aes(x=OrigDose,y=`50%`)) +
    geom_ribbon(aes(ymin=`5%`,ymax=`95%`),fill="#7AD151FF") + # 90% CI of population incidence
    geom_line() +  # Median population incidence
    geom_point(aes(x=bbmd.hed,y=0.1),size=3) + # BBMD - HED
    geom_segment(aes(x = bbmd.hed*1e-5,xend = bbmd.hed,
                     y = 1e-6,yend = 0.1),linetype="dashed") +
    geom_point(aes(x=10^rsd.06,y=1e-6),shape=18,size=3) +
    scale_x_log10(breaks=10^(-10:2),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + 
    scale_y_log10(breaks=10^(-6:-1),labels=trans_format("log10",math_format(10^.x)),expand=c(0,0)) + theme_classic() +
    coord_cartesian(xlim=c(min(pop.incidence.df$OrigDose),max(pop.incidence.df$OrigDose)),ylim=c(1e-6,0.1))+
    xlab("Dose (mg/kg-d)") + ylab("Population Risk") + annotation_logticks(sides="b") +
    theme(legend.position="none",panel.border = element_rect(colour = "black",fill=NA))
  
  plot <- ggarrange(p1,p2,p3,p4,ncol=1)
  annotate_figure(plot,top=text_grob(chem.name[s,"Chemical.Name"],face="bold")) # For long chemical name, size=9
  
  sumfigfolder <- "figures/Summary-Figures"
  ggsave(file.path(sumfigfolder,paste("Index",s,".pdf")),width=4,height=8)
  
}
