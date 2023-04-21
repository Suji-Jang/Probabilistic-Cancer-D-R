library(ggplot2)
library(ggh4x)
library(viridis)
library(scales)
library(ggpubr)
library(lsr)
library(data.table)

# Load data and functions
functionfolder <- "functions"
source(file.path(functionfolder,"Pop_incidence_functions.R"))
source(file.path(functionfolder,"BMD_functions.R"))
source(file.path(functionfolder,"HDMI_functions.R"))
source(file.path(functionfolder,"Extra_risk_functions.R"))
source(file.path(functionfolder,"RSD_functions.R"))
bmdfolder <- "BMD-data"
dose.max <- fread(file.path(bmdfolder,"BMD_data.csv"))
total.data <- fread(file.path(bmdfolder,"Sample_parameters.csv"))
bw.a.df <- fread(file.path(bmdfolder,"BW_data.csv"))
resultsfolder <- "results"
prob.RSD6 <- fread(file.path(resultsfolder,"RSD06_quants.csv"))
hdmi.df <- fread(file.path(resultsfolder,"Supple Table 5 BBMD MA HDMI results.csv"))
figuresfolder <- "figures"

hd.1e6.01.df <- subset(hdmi.df,(hdmi.df$M_indivrisk==1e-6&hdmi.df$I_Incidence==0.01))
hd.1e5.1e3.df <- subset(hdmi.df,(hdmi.df$M_indivrisk==1e-5&hdmi.df$I_Incidence==0.001))
hd.1e4.1e4.df <- subset(hdmi.df,(hdmi.df$M_indivrisk==1e-4&hdmi.df$I_Incidence==1e-4))

dose.max <- aggregate(dose.max$Dose, by = list(dose.max$Study_Index), max)
colnames(dose.max) <- c("Index","Dose.max")
dose.max <- dose.max$Dose.max

bw.a.df <- bw.a.df[,c("Index","bw.a")]
colnames(bw.a.df) <- c("Dataset","bw.a")
bw.h <- 70

# Population incidence from HDMI (For M = 1E-6, I = 0.01)
hd.pop.inc.1e6.df <- data.frame()
for (i in 1:length(dose.max)){
  datasetnum <- i
  samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
  bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
  
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a, # Animal Body Weight
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  
  Dose.conv <- hd.1e6.01.df[datasetnum,"HDMI.p05"]/dose.max[datasetnum]
  Dose.conv <- Dose.conv$HDMI.p05
  pop.incidence.df <- get.pop.incidence(Dose.conv,
                                        uncertainty.samples,
                                        n.pop = 1e4,
                                        prob = c(0.05,0.5,0.95))
  pop.incidence.df$Index <- datasetnum
  hd.pop.inc.1e6.df <- rbind(hd.pop.inc.1e6.df,pop.incidence.df)
}

hd.pop.inc.1e6.df.plot <- hd.pop.inc.1e6.df[order(hd.pop.inc.1e6.df$`95%`,hd.pop.inc.1e6.df$Index),]
hd.pop.inc.1e6.df.plot$order.index <- 1:nrow(hd.pop.inc.1e6.df.plot)
hd.pop.inc.1e6.df.plot <- hd.pop.inc.1e6.df.plot[,c("order.index","5%","50%","95%")]
hd.pop.inc.1e6.df.plot$point.50 <- ifelse(hd.pop.inc.1e6.df.plot$`50%`<=2e-7,0,hd.pop.inc.1e6.df.plot$`50%`) # Run to convert small values to 0

# Population incidence from HDMI (For M = 1E-5, I = 1E-3)
hd.pop.inc.1e5.df <- data.frame()
for (i in 1:length(dose.max)){
  datasetnum <- i
  samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
  bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
  
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a, # Animal Body Weight
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  
  Dose.conv <- hd.1e5.1e3.df[datasetnum,"HDMI.p05"]/dose.max[datasetnum]
  Dose.conv <- Dose.conv$HDMI.p05
  pop.incidence.df <- get.pop.incidence(Dose.conv,
                                        uncertainty.samples,
                                        n.pop = 1e4,
                                        prob = c(0.05,0.5,0.95))
  pop.incidence.df$Index <- datasetnum
  hd.pop.inc.1e5.df <- rbind(hd.pop.inc.1e5.df,pop.incidence.df)
}

hd.pop.inc.1e5.df.plot <- hd.pop.inc.1e5.df[order(hd.pop.inc.1e5.df$`95%`,hd.pop.inc.1e5.df$Index),]
hd.pop.inc.1e5.df.plot$order.index <- 1:nrow(hd.pop.inc.1e5.df.plot)
hd.pop.inc.1e5.df.plot <- hd.pop.inc.1e5.df.plot[,c("order.index","5%","50%","95%")]
hd.pop.inc.1e5.df.plot$point.50 <- ifelse(hd.pop.inc.1e5.df.plot$`50%`<=2e-7,0,hd.pop.inc.1e5.df.plot$`50%`) # Run to convert small values to 0

# Population incidence from HDMI (For M = 1E-4, I = 1E-4)
hd.pop.inc.1e4.df <- data.frame()
for (i in 1:length(dose.max)){
  datasetnum <- i
  samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
  bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
  
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a, # Animal Body Weight
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  
  Dose.conv <- hd.1e4.1e4.df[datasetnum,"HDMI.p05"]/dose.max[datasetnum]
  Dose.conv <- Dose.conv$HDMI.p05
  pop.incidence.df <- get.pop.incidence(Dose.conv,
                                        uncertainty.samples,
                                        n.pop = 1e4,
                                        prob = c(0.05,0.5,0.95))
  pop.incidence.df$Index <- datasetnum
  hd.pop.inc.1e4.df <- rbind(hd.pop.inc.1e4.df,pop.incidence.df)
}

hd.pop.inc.1e4.df.plot <- hd.pop.inc.1e4.df[order(hd.pop.inc.1e4.df$`95%`,hd.pop.inc.1e4.df$Index),]
hd.pop.inc.1e4.df.plot$order.index <- 1:nrow(hd.pop.inc.1e4.df.plot)
hd.pop.inc.1e4.df.plot <- hd.pop.inc.1e4.df.plot[,c("order.index","5%","50%","95%")]
hd.pop.inc.1e4.df.plot$point.50 <- ifelse(hd.pop.inc.1e4.df.plot$`50%`<=2e-7,0,hd.pop.inc.1e4.df.plot$`50%`) # Run to convert small values to 0


# Plot
prob.RSD6 <- prob.RSD6[,c("Dataset","X5.")]
colnames(prob.RSD6) <- c("Index","Prob.5")

prob.rfd.rsd <- data.frame(Index=1:255,RSD6=prob.RSD6$Prob.5,HDMI.1E6=hd.1e6.01.df$HDMI.p05,HDMI.1E5=hd.1e5.1e3.df$HDMI.p05,
                           HDMI.1E4=hd.1e4.1e4.df$HDMI.p05)
colnames(prob.rfd.rsd) <- c("Index","RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")
prob.rfd.rsd[,"RSD6"] <- 10^prob.rfd.rsd[,"RSD6"]

prob.rfd.rsd$HDMI.1E6.RSD6 <- prob.rfd.rsd$HDMI.1E6/prob.rfd.rsd$RSD6
prob.rfd.rsd$HDMI.1E5.RSD6 <- prob.rfd.rsd$HDMI.1E5/prob.rfd.rsd$RSD6
prob.rfd.rsd$HDMI.1E4.RSD6 <- prob.rfd.rsd$HDMI.1E4/prob.rfd.rsd$RSD6

# Figure S6A - Errorbar of population incidence at HD_1E-6^0.01
ymax <- max(c(hd.pop.inc.1e4.df$`95%`,hd.pop.inc.1e5.df$`95%`,hd.pop.inc.1e6.df$`95%`))
p1 <- ggplot(hd.pop.inc.1e6.df.plot) +
  geom_point(aes(x=order.index,y=point.50),color="#27AD81FF") +
  geom_linerange(mapping=aes(x=order.index,y=`50%`,ymin=`5%`, ymax = `95%`),color="#27AD81FF",alpha=0.2)+
  scale_y_log10(breaks=10^(-6:1),labels=trans_format("log10",math_format(10^.x))) + 
  coord_cartesian(ylim=c(5e-7,ymax)) +
  ylab(bquote("Population Incidence at "*HD["1E-6"]^"1E-2")) + labs(tag="A") +
  scale_color_viridis(discrete=TRUE,direction=-1) + theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        plot.margin=margin(0.2,0.2,1.22,0.2,"cm"))

# Figure S6D - Errorbar of population incidence at HD_1E-5^0.001
p4 <- ggplot(hd.pop.inc.1e5.df.plot) +
  geom_point(aes(x=order.index,y=point.50),color="#440154FF") +
  geom_linerange(mapping=aes(x=order.index,y=`50%`,ymin=`5%`, ymax = `95%`),color="#440154FF",alpha=0.2)+
  scale_y_log10(breaks=10^(-6:1),labels=trans_format("log10",math_format(10^.x))) + 
  coord_cartesian(ylim=c(5e-7,ymax)) +
  ylab(bquote("Population Incidence at "*HD["1E-5"]^"1E-3")) + labs(tag="D") +
  scale_color_viridis(discrete=TRUE,direction=-1) + theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        plot.margin=margin(0.2,0.2,1.22,0.2,"cm"))

# Figure S6G - Errorbar of population incidence at HD_1E-4^1E-4
p7 <- ggplot(hd.pop.inc.1e4.df.plot) +
  geom_point(aes(x=order.index,y=point.50),color="#365D8DFF") +
  geom_linerange(mapping=aes(x=order.index,y=`50%`,ymin=`5%`, ymax = `95%`),color="#365D8DFF",alpha=0.2)+
  scale_y_log10(breaks=10^(-6:1),labels=trans_format("log10",math_format(10^.x))) + 
  coord_cartesian(ylim=c(5e-7,ymax)) +
  ylab(bquote("Population Incidence at "*HD["1E-4"]^"1E-4")) + labs(tag="G") +
  scale_color_viridis(discrete=TRUE,direction=-1) + theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        plot.margin=margin(0.2,0.2,1.22,0.2,"cm"))

# Figure S6B - Scatterplot of HD_1E-6^0.01 vs. 5th perc. Prob. RSD_6
p2 <- ggplot(prob.rfd.rsd,aes(x=RSD6,y=HDMI.1E6)) + geom_point(color="#27AD81FF") +
  geom_abline() + labs(tag="B") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote("Prob. RSD6 " * (mg/kg-d))) + ylab(bquote(HD["1E-6"]^"1E-2" * (mg/kg-d))) +
  theme(panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))

# Figure S6E - Scatterplot of HD_1E-5^1E-3 vs. 5th perc. Prob. RSD_6
p5 <- ggplot(prob.rfd.rsd,aes(x=RSD6,y=HDMI.1E5)) + geom_point(color="#440154FF") +
  geom_abline() + labs(tag="E") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote("Prob. RSD6 " * (mg/kg-d))) + ylab(bquote(HD["1E-5"]^"1E-3" * (mg/kg-d))) +
  theme(panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))

# Figure S6H - Scatterplot of HD_1E-4^1E-4 vs. 5th perc. Prob. RSD_6
p8 <- ggplot(prob.rfd.rsd,aes(x=RSD6,y=HDMI.1E4)) + geom_point(color="#365D8DFF") +
  geom_abline() + labs(tag="H") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("RSD6","HDMI.1E6","HDMI.1E5","HDMI.1E4")])),
                minor_breaks = 10^seq(-12,0), guide = "axis_minor",
                breaks=10^seq(-12,0,by=3),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote("Prob. RSD6 " * (mg/kg-d))) + ylab(bquote(HD["1E-4"]^"1E-4" * (mg/kg-d))) +
  theme(panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.24, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))

# Figure S6C - Histogram of ratio HD_1E-6^0.01 : 5th perc. Prob. RSD_6
p3 <- ggplot(prob.rfd.rsd,aes(x=HDMI.1E6.RSD6)) + geom_histogram(fill="#27AD81FF") + labs(tag="C") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")])),
                minor_breaks = 10^seq(-4,8), guide = "axis_minor",
                breaks=10^seq(-4,8,by=2),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote(HD["1E-6"]^"1E-2"*" / Prob. RSD6")) + geom_vline(xintercept = 1,linetype="dashed") +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))

# Figure S6F - Histogram of ratio HD_1E-5^1E-3 : 5th perc. Prob. RSD_6
p6 <- ggplot(prob.rfd.rsd,aes(x=HDMI.1E5.RSD6)) + geom_histogram(fill="#440154FF") + labs(tag="F") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")])),
                minor_breaks = 10^seq(-4,8), guide = "axis_minor",
                breaks=10^seq(-4,8,by=2),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote(HD["1E-5"]^"1E-3"*" / Prob. RSD6")) + geom_vline(xintercept = 1,linetype="dashed") +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))

# Figure S6I - Histogram of ratio HD_1E-4^1E-4 : 5th perc. Prob. RSD_6
p9 <- ggplot(prob.rfd.rsd,aes(x=HDMI.1E4.RSD6)) + geom_histogram(fill="#365D8DFF") + labs(tag="I") +
  scale_x_log10(limits=c(min(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")]),
                         max(prob.rfd.rsd[complete.cases(prob.rfd.rsd),][,c("HDMI.1E6.RSD6","HDMI.1E5.RSD6","HDMI.1E4.RSD6")])),
                minor_breaks = 10^seq(-4,8), guide = "axis_minor",
                breaks=10^seq(-4,8,by=2),labels=trans_format("log10",math_format(10^.x))) +
  theme_classic() + xlab(bquote(HD["1E-4"]^"1E-4"*" / Prob. RSD6")) + geom_vline(xintercept = 1,linetype="dashed") +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95),
        ggh4x.axis.ticks.length.minor = rel(1))


ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9) 
ggsave(file.path(figuresfolder,"Supple Figure 6 - HDMI.pdf"),width=10,height=10)
