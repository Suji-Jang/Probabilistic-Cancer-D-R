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
  n.samp <- 2000
  n.pop <- 1e4

  RSD.df <- data.frame()
  # RSD.etasq.df <- data.frame()
  RSD.samps.df <- data.frame()
  for (datasetnum in 1:nrow(dosemax.df)) {
    cat(datasetnum,"... ")
    ## Get BBMD parameters
    samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
    bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
    ## Uncertainty samples
    set.seed(3.14159)
    uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                   bw.a = bw.a, # Animal Body Weight
                                                   bw.h = bw.h,
                                                   n.samp=n.samp)
    ## Get RSD samples
    t1 <- Sys.time()
    RSD.samples <- calc.RSD.samples(pop.incidence = risklevel, # For RSD at 1e-6
                                    dose.max = dosemax.df$Dose[datasetnum], # For dose scaling 
                                    bmdu10 = output.bbmd$`MA-ext10U`[datasetnum],
                                    uncertainty.samples, # One uncertainty sample
                                    n.pop=n.pop)
    print(Sys.time()-t1)
    RSD.samps.df <- rbind(RSD.samps.df,
                            cbind(data.frame(Dataset=datasetnum,as.data.frame(t(RSD.samples)))))
    # Get quantiles
    RSD.quants <- log10(quantile(RSD.samples,prob=seq(0.01,0.99,0.01),na.rm=TRUE))
    RSD.df <- rbind(RSD.df,
                    cbind(data.frame(Dataset=datasetnum,as.data.frame(t(RSD.quants)))))
  }
  fwrite(RSD.samps.df,file.path(resultsfolder,paste0(RSD.label,"samples.csv")))
  fwrite(RSD.df,file.path(resultsfolder,paste0(RSD.label,".quants.csv")))
}
