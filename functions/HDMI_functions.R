# HDMI_functions.R

### Calculate bmd samples
calc.bmd.samples <- function(
  bmr,                # bmr (single value) at which HDMI is calculated
  uncertainty.samples # Uncertainty samples
) {
  # Number of uncertainty samples
  n.samp<-length(uncertainty.samples$sigmaH);
  # bmd for each sample 
  bmd.samples <- numeric(n.samp)
  for (j in 1:n.samp) {
    bmd.samples[j] <- bmd.func(bmr,
                               as.numeric(uncertainty.samples$parms[j,]),
                               uncertainty.samples$modelname[j])
  }      
  bmd.samples
}

### Get incidence of individual risk

get.indivrisk.incidence <- function(
  indivrisk = 0.01, # individual risk value (extra risk = bmr)
  I.vec = 10^seq(-6,0,0.1), # incidence values
  uncertainty.samples, # samples from before
  prob=c(0.05,0.5,0.95), # Quantiles to report
  dose.max = 1
) {
  # cat("Calculating HD50 samples...\n");
  bmd.samples <- dose.max * calc.bmd.samples(indivrisk,uncertainty.samples)
  HD50.samples <- bmd.samples * 
    uncertainty.samples$DAF/(uncertainty.samples$AHU*uncertainty.samples$OU)
  hdmi.vec <- numeric(0)
  for (k in 1:length(I.vec)) {
    z <- qnorm(I.vec[k]);
    hdmi.samples<-HD50.samples*exp(z*uncertainty.samples$sigmaH);
    # Quantiles
    hdmi.quantiles<-as.data.frame(t(quantile(hdmi.samples,prob=prob,na.rm=TRUE)))
    # Etasq=Percent contribution to variance
    tmp <- data.frame(hdmi=log(hdmi.samples),
                      bmd=log(bmd.samples),
                      sigmaH=log(uncertainty.samples$sigmaH),
                      DAF=log(uncertainty.samples$DAF),
                      AHU=log(uncertainty.samples$AHU),
                      OU=log(uncertainty.samples$OU))
    hdmi.aov<-aov(hdmi~bmd+sigmaH+DAF+AHU+OU,data=tmp)
    hdmi.etasq<-as.data.frame(t(etaSquared(hdmi.aov)[,"eta.sq"]))
    # Combined
    hdmi.vec<-rbind(hdmi.vec,cbind(hdmi.quantiles,hdmi.etasq))
  }
  # Return table of quantiles
  cbind(data.frame(indivrisk=indivrisk,Incidence=I.vec,as.data.frame(hdmi.vec)))
}
