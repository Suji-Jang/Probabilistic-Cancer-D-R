library(gplots)
library(RColorBrewer)
library(ggplot2)
library(lsr)
library(data.table)
library(GGally)
library(ggcorrplot)

bmdfolder <- "BMD-data"
weight.df <- fread(file.path(bmdfolder,"Model_weight_092022.csv"))
bmds.df <- fread(file.path(bmdfolder,"BMDS_final_data_092822.csv"))
ma.df <- fread(file.path(bmdfolder,"MA_RSD_092822.csv"))
resultsfolder <- "results"
rsd6.df <- fread(file.path(resultsfolder,"RSD06.quants.csv"))
figuresfolder <- "figures"

bmdl.df <- data.frame(Index=1:255,MA.BMDL=ma.df[,"BMDL.10"],BMDS.BMDL=bmds.df[,"BMDL"])
colnames(bmdl.df) <- c("Index","MA.BMDL","BMDS.BMDL")

rsd6.df <- rsd6.df[,c("Dataset","X5.","X50.","X95.")]
rsd6.df[,2:4] <- 10^rsd6.df[,2:4]
rsd6.df$MA.RSD6 <- ma.df$RSD.06.L
colnames(rsd6.df) <- c("Index","RSD6.Prob.05","RSD6.Prob.50","RSD6.Prob.95","RSD6.MA.Linear")

heatmap.df <- weight.df[,2:9]
heatmap.df$row.sum <- rowSums(heatmap.df)
heatmap.df[,1:8] <- weight.df[,2:9]/heatmap.df$row.sum
heatmap.df <- heatmap.df[,1:8]
colnames(heatmap.df) <- c("QL","Log","Pro","Wei","MS2","LLog","LPro","DH")

BMDL.Ratio <- log10(bmdl.df$BMDS.BMDL/bmdl.df$MA.BMDL)
heatmap.df$BMDL.Ratio <- BMDL.Ratio
RSD.Ratio <- log10(rsd6.df$RSD6.MA.Linear/rsd6.df$RSD6.Prob.05)
heatmap.df$RSD.Ratio <- RSD.Ratio

scattermat <- ggpairs(heatmap.df)
ggsave(file.path(figuresfolder,"Figure 6 - Scatter.pdf"),scattermat,height=5,width=5,scale=2)
corrmat <- ggcorrplot(round(cor(heatmap.df),2), type = "lower", show.diag = TRUE,
                      lab = TRUE,)
ggsave(file.path(figuresfolder,"Figure 6 - Corr.pdf"),corrmat,height=5,width=5,scale=1.5)

# Rescale 0 to 1 for heatmap
heatmap.df$BMDL.Ratio <- (BMDL.Ratio-min(BMDL.Ratio))/(max(BMDL.Ratio)-min(BMDL.Ratio))
heatmap.df$RSD.Ratio <- (RSD.Ratio-min(RSD.Ratio))/(max(RSD.Ratio)-min(RSD.Ratio))

mypalette <- (colorRampPalette(brewer.pal(9, "Reds"))(100))

heatmap.data <- t(heatmap.df)
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - All.pdf"),height=6,width=20)
heatmap.2(heatmap.data,col=mypalette, trace="none",density.info = "density",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,7))
dev.off()

heatmap.data.BMDL <- t(heatmap.df)[c(1:8,9),]
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - BMDL.pdf"),height=6,width=20)
heatmap.2(heatmap.data.BMDL,col=mypalette, trace="none",density.info = "density",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,7))
dev.off()

heatmap.data.RSD <- t(heatmap.df)[c(1:8,10),]
pdf(file.path(figuresfolder,"Figure 6 - Heatmap - RSD.pdf"),height=6,width=20)
heatmap.2(heatmap.data.RSD,col=mypalette, trace="none",density.info = "density",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,7))
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
