# Get_HDMI.R
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

# Incidence from 1e-6 to 0.99
I.vec <- c(10^seq(-6,-2,0.1),seq(0.02,0.99,0.01))

# BMR = individual risk from 1e-6 to 0.1
bmr.vec <- 10^seq(-6,-1)

allhdmi.df <- data.frame()
for (i in 1:271){
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

write.csv(allhdmi.df,file.path(resultsfolder,"All HDMI - Index 1-271.csv"),row.names = FALSE)
