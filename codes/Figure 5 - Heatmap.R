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

# Load data and functions
bmdfolder <- "BMD-data"
weight.df <- fread(file.path(bmdfolder,"Model_weight.csv"))
bmds.df <- fread(file.path(bmdfolder,"BMDS_results.csv"))
resultsfolder <- "results"
ma.df <- fread(file.path(resultsfolder,"MA_RSD_results.csv"))
rsd6.df <- fread(file.path(resultsfolder,"RSD06_quants.csv"))
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

# Rescale using scale function
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

# Create color palette - Blue to Red
mypalette <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))

# Plot
heatmap.data <- t(heatmap.df)
pdf(file.path(figuresfolder,"Figure 5 - Heatmap.pdf"),height=6,width=20)
heatmap.2(heatmap.data,col=mypalette, trace="none",density.info = "density",
          Rowv=FALSE,dendrogram="column",lwid=c(0.8,8),lhei=c(1.5,5),margins=c(2,10))
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
