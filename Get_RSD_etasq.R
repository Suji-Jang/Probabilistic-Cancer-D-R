library(data.table)
library(ggplot2)
library(ggpubr)
library(lsr) # Needed for proportion of variation
functionfolder <- "functions"
source(file.path(functionfolder,"Pop_incidence_functions.R"))
source(file.path(functionfolder,"BMD_functions.R"))
source(file.path(functionfolder,"HDMI_functions.R"))
source(file.path(functionfolder,"Extra_risk_functions.R"))
source(file.path(functionfolder,"RSD_functions.R"))
bmdfolder <- "BMD-data"
total.data <- fread(file.path(bmdfolder,"Total_data_121821.csv"))
output.bbmd <- fread(file.path(bmdfolder,"Summary_Output_data.csv")) 
bmddata <- fread(file.path(bmdfolder,"BMD_data_121821.csv"))
bmd.info <- fread(file.path(bmdfolder,"BMD_w_info_112921.csv"))
bw.a.df <- fread(file.path(bmdfolder,"BW_data_final_121821.csv"))
resultsfolder <- "results"

dosemax.df <- aggregate(Dose ~ Study_Index,max,data=bmddata)
bw.a.df <- bw.a.df[,c("Index","bw.a")]
colnames(bw.a.df) <- c("Dataset","bw.a")
bw.h <- 70

for (k in 1:3) {
  if (k==1) {
    risklevel <- 1e-4
    RSD.label <- "RSD04"
  } else if (k==2) {
    risklevel <- 1e-5
    RSD.label <- "RSD05"
  } else if (k==3) {
    risklevel <- 1e-6
    RSD.label <- "RSD06"
  }
  RSD.samps.df <- fread(file.path(resultsfolder,paste0(RSD.label,"samples.csv")))
  n.samp <- ncol(RSD.samps.df)-1
  RSD.etasq.df <- data.frame()
  for (datasetnum in 1:nrow(dosemax.df)) {
    cat(datasetnum,"... ")
    ## Get BBMD parameters
    samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
    bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
    ## Uncertainty samples
    set.seed(314159)
    uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                   bw.a = bw.a, # Animal Body Weight
                                                   bw.h = bw.h,
                                                   n.samp=n.samp)
    RSD.samples <- as.numeric(RSD.samps.df[datasetnum,-1])
    ## Get contribution to variance
    RSD.etasq <- calc.RSD.etasq(RSD.samples,
                                uncertainty.samples)
    RSD.etasq.df <- rbind(RSD.etasq.df,
                          cbind(data.frame(Dataset=datasetnum),RSD.etasq))
  }
  fwrite(RSD.etasq.df,file.path(resultsfolder,paste0(RSD.label,".etasq.csv")))
}
