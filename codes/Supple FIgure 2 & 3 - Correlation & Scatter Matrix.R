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

# Plot
corrmat <- ggcorrplot(round(cor(heatmap.df),2), type = "lower", show.diag = TRUE,
                      lab = TRUE)
ggsave(file.path(figuresfolder,"Supple Figure 2 - Correlation.pdf"),corrmat,height=5,width=5,scale=1.5)
scattermat <- ggpairs(heatmap.df)
ggsave(file.path(figuresfolder,"Supple Figure 3 - Scatter Matrix.pdf"),scattermat,height=5,width=5,scale=2)
