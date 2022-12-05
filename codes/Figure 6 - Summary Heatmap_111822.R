library(ggplot2)
library(gplots)
library(RColorBrewer)
library(viridis)
library(lsr)
library(data.table)
library(GGally)
library(ggcorrplot)
library(scales)
library(cluster)

bmdfolder <- "BMD-data"
weight.df <- fread(file.path(bmdfolder,"Model_weight_092022.csv"))
bmds.df <- fread(file.path(bmdfolder,"BMDS_final_data_092822.csv"))
ma.df <- fread(file.path(bmdfolder,"MA_RSD_092822.csv"))
resultsfolder <- "results"
rsd6.df <- fread(file.path(resultsfolder,"RSD06.quants.csv"))
figuresfolder <- "figures"

bmdl.df <- data.frame(Index=1:255,MA.BMDL=ma.df[,"BMDL.10"],MA.BMDU=ma.df[,"BMDU.10"],BMDS.BMDL=bmds.df[,"BMDL"])
colnames(bmdl.df) <- c("Index","MA.BMDL","MA.BMDU","BMDS.BMDL")

rsd6.df <- rsd6.df[,c("Dataset","X5.","X50.","X95.")]
rsd6.df[,2:4] <- 10^rsd6.df[,2:4]
rsd6.df$MA.RSD6 <- ma.df$RSD.06.L
colnames(rsd6.df) <- c("Index","RSD6.Prob.05","RSD6.Prob.50","RSD6.Prob.95","RSD6.MA.Linear")

heatmap.df <- weight.df[,2:9]
heatmap.df$row.sum <- rowSums(heatmap.df)
heatmap.df[,1:8] <- weight.df[,2:9]/heatmap.df$row.sum
heatmap.df <- heatmap.df[,1:8]
colnames(heatmap.df) <- c("QL","Log","Pro","Wei","MS2","LLog","LPro","DH")

BBMD.BMDS.Ratio <- log10(bmdl.df$BMDS.BMDL/bmdl.df$MA.BMDL)
heatmap.df$BBMD.BMDS.Ratio <- BBMD.BMDS.Ratio
BMD.CI.Ratio <- log10(bmdl.df$MA.BMDU/bmdl.df$MA.BMDL)
heatmap.df$BMD.CI.Ratio <- BMD.CI.Ratio
Prob.Linear.Ratio <- log10(rsd6.df$RSD6.MA.Linear/rsd6.df$RSD6.Prob.05)
heatmap.df$Prob.Linear.Ratio <- Prob.Linear.Ratio
RSD.CI.Ratio <- log10(rsd6.df$RSD6.Prob.95/rsd6.df$RSD6.Prob.05)
heatmap.df$RSD.CI.Ratio <- RSD.CI.Ratio

# All parameters
scattermat <- ggpairs(heatmap.df)
ggsave(file.path(figuresfolder,"Figure 6 - Scatter.pdf"),scattermat,height=5,width=5,scale=2)
corrmat <- ggcorrplot(round(cor(heatmap.df),2), type = "lower", show.diag = TRUE,
                      lab = TRUE)
ggsave(file.path(figuresfolder,"Figure 6 - Corr.pdf"),corrmat,height=5,width=5,scale=1.5)

# Ratio Only
ratio.heatmap.df <- heatmap.df[,c("BBMD.BMDS.Ratio","BMD.CI.Ratio","Prob.Linear.Ratio","RSD.CI.Ratio")]
scattermat <- ggpairs(ratio.heatmap.df)
ggsave(file.path(figuresfolder,"Figure 6 - Scatter_Ratio.pdf"),scattermat,height=5,width=5,scale=2)
corrmat <- ggcorrplot(round(cor(ratio.heatmap.df),2), type = "lower", show.diag = TRUE,
                      lab = TRUE,)
ggsave(file.path(figuresfolder,"Figure 6 - Corr_Ratio.pdf"),corrmat,height=5,width=5,scale=1.5)

# Rescale 0 to 1 for heatmap
heatmap.df$BBMD.BMDS.Ratio <- (BBMD.BMDS.Ratio-min(BBMD.BMDS.Ratio))/(max(BBMD.BMDS.Ratio)-min(BBMD.BMDS.Ratio))
heatmap.df$BMD.CI.Ratio <- (BMD.CI.Ratio-min(BMD.CI.Ratio))/(max(BMD.CI.Ratio)-min(BMD.CI.Ratio))
heatmap.df$Prob.Linear.Ratio <- (Prob.Linear.Ratio-min(Prob.Linear.Ratio))/(max(Prob.Linear.Ratio)-min(Prob.Linear.Ratio))
heatmap.df$RSD.CI.Ratio <- (RSD.CI.Ratio-min(RSD.CI.Ratio))/(max(RSD.CI.Ratio)-min(RSD.CI.Ratio))

mypalette <- (colorRampPalette(brewer.pal(9, "Reds"))(100))

# Rescale with scale function
heatmap.df[,1] <- scale(heatmap.df[,1])
heatmap.df[,2] <- scale(heatmap.df[,2])
heatmap.df[,3] <- scale(heatmap.df[,3])
heatmap.df[,4] <- scale(heatmap.df[,4])
heatmap.df[,5] <- scale(heatmap.df[,5])
heatmap.df[,6] <- scale(heatmap.df[,6])
heatmap.df[,7] <- scale(heatmap.df[,7])
heatmap.df[,8] <- scale(heatmap.df[,8])
heatmap.df$BBMD.BMDS.Ratio <- scale(BBMD.BMDS.Ratio)
heatmap.df$BMD.CI.Ratio <- scale(BMD.CI.Ratio)
heatmap.df$Prob.Linear.Ratio <- scale(Prob.Linear.Ratio)
heatmap.df$RSD.CI.Ratio <- scale(RSD.CI.Ratio)
a <- clusGap(heatmap.df,FUN=kmeans,K.max=)

mypalette <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))

# Mapping
heatmap.data <- t(heatmap.df)
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - All_scale_func.pdf"),height=6,width=20)
heatmap.2(heatmap.data,col=mypalette, trace="none",density.info = "density",
          Rowv=FALSE,dendrogram="column",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,10))
dev.off()

heatmap.data.BMDL <- t(heatmap.df)[c(1:8,9,10),]
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - BMDL_scale_func.pdf"),height=6,width=20)
heatmap.2(heatmap.data.BMDL,col=mypalette, trace="none",density.info = "density",
          Rowv=FALSE,dendrogram="column",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,10))
dev.off()

heatmap.data.RSD <- t(heatmap.df)[c(1:8,11,12),]
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - RSD_scale_func.pdf"),height=6,width=20)
heatmap.2(heatmap.data.RSD,col=mypalette, trace="none",density.info = "density",
          Rowv=FALSE,dendrogram="column",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,10))
dev.off()

heatmap.data.Ratio <- t(heatmap.df)[c(9:12),]
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - Ratio_scale_func.pdf"),height=4,width=20)
heatmap.2(heatmap.data.Ratio,col=mypalette, trace="none",density.info = "density",
          Rowv=FALSE,dendrogram="column",lwid=c(0.8,8),lhei=c(0.4,0.6),margins=c(2,10),cexRow=1)
dev.off()

# ANOVA
bmdl.summary = lm(BMDL.Ratio ~ QL + Log + Pro + Wei + MS2 + LLog + LPro + DH - 1, data = heatmap.df)
rsd.summary = lm(RSD.Ratio ~ QL + Log + Pro + Wei + MS2 + LLog + LPro + DH - 1, data = heatmap.df)

aovbmdl <- aov(bmdl.summary)
print(summary(bmdl.summary))
print(summary(aovbmdl)) 
print(etaSquared(aovbmdl))

aovrsd <- aov(rsd.summary)
print(summary(rsd.summary))
print(summary(aovrsd)) 
print(etaSquared(aovrsd))
