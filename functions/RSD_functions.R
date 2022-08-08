### Calculate RSD using root finding
calc.RSD.samples <- function(
  pop.incidence = 1e-6, # For RSD at 1e-6
  dose.max = 1, # For dose scaling 
  bmdu10 = 1, # For range finding
  uncertainty.samples, # Uncertainty samples
  n.samp = length(uncertainty.samples$sigmaH),
  n.pop=1e4,	     # Number of individuals in the population sample
  doplot = FALSE
) {
  # Number of uncertainty samples
  
  RSD.samples <- numeric(n.samp)
  RSD.samples[1:n.samp] <- NA
  for (j in 1:n.samp) {
    func <- function(log10dose) { # 0 when pop.incidence.onesample = pop.incidence
      set.seed(314159)
      dose <- 10^log10dose
      dose.median <- dose * uncertainty.samples$AHU[j]*
        uncertainty.samples$OU[j]/uncertainty.samples$DAF[j];
      dose.distrib <- dose.median*exp(rnorm(n.pop)*uncertainty.samples$sigmaH[j]);
      pop.incidence.onesample<- mean(er.func(dose.distrib,
                                             as.numeric(uncertainty.samples$parms[j,]),
                                             uncertainty.samples$modelname[j]))
      delta <- pop.incidence.onesample/pop.incidence - 1
    }
    if (func(Inf) > 0 & func(-Inf) < 0) {
      try({
        low <- log10(bmdu10) - 20
        high <- log10(bmdu10) + 5
        soln <- uniroot(f = func,interval = c(low,high));
        RSD.samples[j] <- dose.max * 10^soln$root;
        if (abs(soln$f.root) > 0.01) print(paste(soln$f.root,func(soln$root)))
        if (doplot) {
          xx <- log10(bmdu10)+seq(-20,5,0.1)
          yy <- xx*0
          for (k in 1:length(xx)) yy[k] <- func(xx[k])
          plot(dose.max*10^xx,1e-6*(yy+1),log="xy",type="l",main=paste("Sample",j))
          points(RSD.samples[j],pop.incidence)
          abline(h=1e-6)
        }
      })
    }
  }
  return(RSD.samples)
}
### Calculate contribution to variance eta-sq 
calc.RSD.etasq <- function(RSD.samples,
                           uncertainty.samples,
                           printresult=FALSE) {
  tmp <- data.frame(RSD=log10(RSD.samples),
                    Model=uncertainty.samples$modelname, # Model shape
                    sigmaH=log10(uncertainty.samples$sigmaH),
                    DAF=log10(uncertainty.samples$DAF),
                    AHU=log10(uncertainty.samples$AHU))
  if (length(unique(tmp$Model))==1) { # Only one model
    RSD.aov<-aov(RSD~sigmaH+DAF+AHU,data=tmp)
    if (printresult) print(summary(RSD.aov))
    etasq <- etaSquared(RSD.aov,anova = TRUE)
    etasq <- rbind(rep(0,ncol(etasq)),etasq)
    rownames(etasq)[1]<-"Model"
  } else {
    RSD.aov<-aov(RSD~Model+sigmaH+DAF+AHU,data=tmp)
    if (printresult) print(summary(RSD.aov))
    etasq <- etaSquared(RSD.aov,anova = TRUE)
  }
  if (printresult) print(etasq)
  RSD.etasq<-cbind(data.frame(RSD.log10sd=sd(log10(RSD.samples)),
                              RSD.log10mad=mad(log10(RSD.samples)),
                              RSD.log10iqr=IQR(log10(RSD.samples)),
                              RSD.log10ci90=diff(quantile(log10(RSD.samples),prob=c(0.05,0.95))),
                              RSD.log10ci95=diff(quantile(log10(RSD.samples),prob=c(0.025,0.975)))
                              ),
                   as.data.frame(t(etasq[,"eta.sq"])))
  return(RSD.etasq)
}

