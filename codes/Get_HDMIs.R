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
total.data <- fread(file.path(bmdfolder,"Sample_parameters.csv"))
bmddata <- fread(file.path(bmdfolder,"BMD_data.csv"))
bw.a.df <- fread(file.path(bmdfolder,"BW_data.csv"))
resultsfolder <- "results"

dosemax.df <- aggregate(Dose ~ Study_Index,max,data=bmddata)
bw.a.df <- bw.a.df[,c("Index","bw.a")]
colnames(bw.a.df) <- c("Dataset","bw.a")
bw.h <- 70

# Incidence from 1e-6 to 0.99
I.vec <- c(10^seq(-6,-2,0.1),seq(0.02,0.99,0.01))

# BMR = individual risk from 1e-6 to 0.1
bmr.vec <- 10^seq(-6,-1)

allhdmi.df <- data.frame()
for (i in 1:nrow(bw.a.df)){
  cat(i,"...")
  datasetnum <- i
  samp.parms <- subset(total.data,Dataset==datasetnum)[,-1] # Remove first column which is dataset #
  bw.a <- subset(bw.a.df,Dataset==datasetnum)$bw.a[1]
  ## Uncertainty samples
  set.seed(3.14159) # use same random seed and n.samp as RSD
  uncertainty.samples <- get.uncertainty.samples(samp.parms,
                                                 bw.a = bw.a, # Animal Body Weight
                                                 bw.h = bw.h,
                                                 n.samp=2000)
  
  # Individual risk of 1e-6, 1e-5... 0.1
  for (indivrisk in bmr.vec) {
    hdmi.df <- get.indivrisk.incidence(indivrisk,
                                       I.vec,
                                       uncertainty.samples,
                                       prob = c(0.05,0.5,0.95),
                                       dose.max = dosemax.df$Dose[i])
    hdmi.df$Index <- datasetnum
    allhdmi.df <- rbind(allhdmi.df,hdmi.df)
  }
}

allhdmi.df <- allhdmi.df[,c("Index",1:9)]
colnames(allhdmi.df) <- c("Index","M_indivrisk","I_Incidence","HDMI.p05","HDMI.p50","HDMI.p95",
                          "HDMI.etasq.bmd","HDMI.etasq.sigmaH","HDMI.etasq.DAF","HDMI.etasq.AHU")
fwrite(allhdmi.df,file.path(resultsfolder,"Supple Table 5 BBMD MA HDMI results.csv"))
