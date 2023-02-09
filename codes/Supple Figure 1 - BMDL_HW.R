library(reshape2)
library(ggpubr)
library(viridis)
library(scales)
library(data.table)

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
figuresfolder <- "figures"

errorbar.df <- data.frame(Index=1:255,MA.BMDL=ma.df[,"BMDL.10"],MA.BMD=ma.df[,"BMD.10"],MA.BMDU=ma.df[,"BMDU.10"],
                          HW.BMDL=hw.df[,"BMDL.10"],BMDS.BMDL=bmds.df[,"BMDL"])
colnames(errorbar.df) <- c("Index","MA.BMDL","MA.BMD","MA.BMDU","HW.BMDL","BMDS.BMDL")

errorbar.df <- errorbar.df[order(errorbar.df$MA.BMD,errorbar.df$Index),]
errorbar.df$order.index <- 1:255

# Plot
errorbar.hw.plot.df <- errorbar.df[order(errorbar.df$MA.BMD,errorbar.df$Index),]
errorbar.hw.plot.df$order.index <- 1:255
errorbar.hw.plot.df <- errorbar.hw.plot.df[,c("order.index","MA.BMD","HW.BMDL")]
colnames(errorbar.hw.plot.df) <- c("order.index","BMD [MA]","BMDL [HW]")

errorbar.hw.plot.df <- melt(errorbar.hw.plot.df,id.vars="order.index")

p1 <- ggplot(errorbar.hw.plot.df) +
  geom_linerange(data=errorbar.df,mapping=aes(x=order.index,y=MA.BMD,ymin = MA.BMDL, ymax = MA.BMDU),
                 alpha=0.4,colour="#440154FF")+
  geom_point(aes(x=order.index,y=value,colour=variable,shape=variable),alpha=0.7) +
  scale_y_log10(limits=c(min(errorbar.df[,c("HW.BMDL","MA.BMDL")])
                         ,max(errorbar.df[,c("HW.BMDL","MA.BMDU")])),
                breaks=10^c(-6:4),labels=trans_format("log10",math_format(10^.x))) +
  annotation_logticks(sides="l") + labs(tag="A") +
  ylab(bquote(BMDL[10])) +
  scale_color_viridis(discrete=TRUE,end=0.7) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
        legend.title=element_blank(),legend.position = c(0.5, 0.1),
        panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.1, 0.95))

hist.df <- data.frame(Index=1:255,MA.BMDL=ma.df[,"BMDL.10"],HW.BMDL=hw.df[,"BMDL.10"],
                      BMDS.BMDL=bmds.df[,"BMDL"])
colnames(hist.df) <- c("Index","MA.BMDL","HW.BMDL","BMDS.BMDL")
hist.df$HW.MA <- hist.df$HW.BMDL/hist.df$MA.BMDL
hist.df$BMDS.MA <- hist.df$BMDS.BMDL/hist.df$MA.BMDL

p2 <- ggplot(hist.df,aes(x=MA.BMDL,y=HW.BMDL)) + geom_point() +
  geom_abline() + labs(tag="B") +
  scale_x_log10(limits=c(min(errorbar.df[complete.cases(errorbar.df),][,c("HW.BMDL","BMDS.BMDL")])
                         ,max(errorbar.df[complete.cases(errorbar.df),][,c("HW.BMDL","BMDS.BMDL")])),
                breaks=10^(-6:3),labels=trans_format("log10",math_format(10^.x))) +
  scale_y_log10(limits=c(min(errorbar.df[complete.cases(errorbar.df),][,c("HW.BMDL","BMDS.BMDL")])
                         ,max(errorbar.df[complete.cases(errorbar.df),][,c("HW.BMDL","BMDS.BMDL")])),
                breaks=10^-(6:3),labels=trans_format("log10",math_format(10^.x))) +
  xlab("BMDL [BBMD]") + ylab("BMDL [HW]") + theme_classic() +
  theme(panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.2, 0.95))

p3 <- ggplot(hist.df,aes(x=HW.MA)) + geom_histogram() + labs(tag="C") +
  scale_x_log10(limits=c(min(hist.df$BMDS.MA),max(hist.df$BMDS.MA)),breaks=10^c(-4,-3,-2,-1,0,1,2,3,4),   
                labels=trans_format("log10",math_format(10^.x))) +
  xlab("BMDL [HW] / BMDL [BBMD]") +
  theme_classic() + geom_vline(xintercept=1,linetype="dashed") +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),
        panel.border = element_rect(colour = "black",fill=NA), 
        plot.tag = element_text(size=20, face="bold"),plot.tag.position = c(0.07, 0.95))

hist.scatter <- ggarrange(p2,p3,heights=c(5,5),widths=c(5,5))
sum.plot <- ggarrange(p1,hist.scatter,heights=c(6,4),widths=c(10,10),ncol=1)
ggsave(file.path(figuresfolder,"Supple Figure 1 - BMDL_HW.pdf"),height=10,width=8)