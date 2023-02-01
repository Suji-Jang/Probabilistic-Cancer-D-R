library(reshape2)
library(ggpubr)
library(viridis)
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
bmds.df <- fread(file.path(bmdfolder,"BMDS_results.csv"))
resultsfolder <- "results"
total.bmds.df <- fread(file.path(resultsfolder,"Supple Table 2 BMDS RSD results.csv"))
ma.df <- fread(file.path(resultsfolder,"MA_RSD_results.csv"))
hw.df <- fread(file.path(resultsfolder,"Supple Table 3 BBMD HW RSD results.csv"))
rsd6.df <- fread(file.path(resultsfolder,"RSD06_quants.csv"))
figuresfolder <- "figures"

rsd6.df <- rsd6.df[,c("Dataset","X5.","X50.","X95.")]
rsd6.df[,c("X5.","X50.","X95.")] <- 10^rsd6.df[,c("X5.","X50.","X95.")]
colnames(rsd6.df) <- c("Index","Prob.5","Prob.50","Prob.95")
rsd6.df$BMDS.RSD6 <- total.bmds.df$RSD.06
rsd6.df$MA.RSD6 <- ma.df$RSD.06.L

rsd6.df <- rsd6.df[order(rsd6.df$Prob.50,rsd6.df$Index),]
rsd6.df$order.index <- 1:nrow(rsd6.df)

# Plot
rsd6.points.df <- rsd6.df[,c("order.index","MA.RSD6","BMDS.RSD6","Prob.50")]
colnames(rsd6.points.df) <- c("order.index","BBMD-MA + Linear Extrapolation",
                              "BMDS + Linear Extrapolation","BBMD-MA + Probabilistic Extapolation")
rsd6.points.df <- melt(rsd6.points.df,id.vars="order.index")

# Figure 3A - Error bar of Prob. RSD vs. BMDS RSD vs. BBMD-Linear extrapolated RSD
p1 <- ggplot(rsd6.points.df) +
  geom_linerange(data=rsd6.df,mapping=aes(x=order.index,y=Prob.50,
                                          ymin = Prob.95, ymax = Prob.5),alpha=0.5,colour="#6DCD59FF")+
  geom_point(aes(x=order.index,y=value,colour=variable,shape=variable),alpha=0.7) +
  scale_y_log10(limits=c(min(rsd6.df[,c("MA.RSD6","BMDS.RSD6","Prob.5")])
                         ,max(rsd6.df[,c("MA.RSD6","BMDS.RSD6","Prob.95")])),
                breaks=10^c(-11:2),labels=trans_format("log10",math_format(10^.x))) +
  annotation_logticks(sides="l") +
  ylab("RSD6") + labs(tag="A") +
  scale_color_viridis(discrete=TRUE,end=0.8) +
  theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
        legend.title=element_blank(),legend.position = c(0.5, 0.1),
        panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.1, 0.95))

hist.RSD6.df <- data.frame(Index=1:nrow(rsd6.df),Fit.RSD6=rsd6.df[,"Prob.5"],MA.RSD6=rsd6.df[,"MA.RSD6"],
                           BMDS.RSD6=rsd6.df[,"BMDS.RSD6"])
colnames(hist.RSD6.df) <- c("Index","Fit.RSD6","MA.RSD6","BMDS.RSD6")
hist.RSD6.df$MA.Prob <- hist.RSD6.df$MA.RSD6/hist.RSD6.df$Fit.RSD6
hist.RSD6.df$BMDS.Prob <- hist.RSD6.df$BMDS.RSD6/hist.RSD6.df$Fit.RSD6

# Figure 3B - Scatterplot of BBMD-Linear extrapolated RSD vs. 5th perc. Prob. RSD
p2 <- ggplot(hist.RSD6.df,aes(x=Fit.RSD6,y=MA.RSD6)) + geom_point() +
  geom_abline() + labs(tag="B") +
  scale_x_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])
                         ,max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])),
                breaks=10^c(-11:-1),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])
                         ,max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])),
                breaks=10^c(-11:-1),labels=trans_format("log10",math_format(10^.x))) + 
  xlab("RSD6 [BBMD-MA + Prob. Extr. (5th Perc.)]") + ylab("RSD6 [BBMD-MA + Linear Extr.]") + theme_classic() +
  theme(panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.2, 0.95))

# Figure 3C - Scatterplot of BMDS RSD vs. 5th perc. Prob. RSD
p3 <- ggplot(hist.RSD6.df,aes(x=Fit.RSD6,y=BMDS.RSD6)) + geom_point() +
  geom_abline() + labs(tag="C") +
  scale_x_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])
                         ,max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])),
                breaks=10^c(-11:-1),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])
                         ,max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("Fit.RSD6","MA.RSD6","BMDS.RSD6")])),
                breaks=10^c(-11:-1),labels=trans_format("log10",math_format(10^.x))) + 
  xlab("RSD6 [BBMD-MA + Prob. Extr. (5th Perc.)]") + ylab("RSD6 [BMDS + Linear Extr.]") + theme_classic() +
  theme(panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.2, 0.95))

# Figure 3D - Histogram of ratio BBMD-Linear extrapolated RSD : 5th perc. Prob. RSD
p4 <- ggplot(hist.RSD6.df,aes(x=MA.Prob)) + geom_histogram() + labs(tag="D") +
  scale_x_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("MA.Prob","BMDS.Prob")]),
                         max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("MA.Prob","BMDS.Prob")])),
                breaks=10^c(-3:5),labels=trans_format("log10",math_format(10^.x))) + 
  xlab(expression(frac("RSD6 [BBMD-MA + Linear Extr.]","RSD6 [BBMD-MA + Prob. Extr. (5th Perc.)]"))) +
  theme_classic() + geom_vline(xintercept=1,linetype="dashed") +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95))

# Figure 3E - Histogram of ratio BMDS RSD : 5th perc. Prob. RSD
p5 <- ggplot(hist.RSD6.df,aes(x=BMDS.Prob)) + geom_histogram() + labs(tag="E") +
  scale_x_log10(limits=c(min(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("MA.Prob","BMDS.Prob")])
                         ,max(hist.RSD6.df[complete.cases(hist.RSD6.df),][,c("MA.Prob","BMDS.Prob")])),
                breaks=10^c(-3:5),labels=trans_format("log10",math_format(10^.x))) +
  xlab(expression(frac("RSD6 [BMDS + Linear Extr.]","RSD6 [BBMD-MA + Prob. Extr. (5th Perc.)]"))) +
  theme_classic() + geom_vline(xintercept=1,linetype="dashed") +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95))

plot.2 <- ggarrange(p2,p3)
plot.3 <- ggarrange(p4,p5)
sum.plot <- ggarrange(p1,plot.2,plot.3,heights=c(6,4,4),widths=c(10,10),ncol=1)
ggsave(file.path(figuresfolder,"Figure 3 - RSD.pdf"),height=12,width=8)
