# Pop_incidence_functions.R
## Revised - includes proportion of variance
### Get uncertainty samples
get.uncertainty.samples <- function(
  samp.parms,    # data frame first column is model name, remaining columns are parameters to be resampled
  bw.a=0.3,		    # animal body weight (default = rat)
  bw.h=70,	    # human body weight (default 70kg)
  alpha.m=0.7,	    # allometric scaling exponent median
  alpha.sd=0.0243,   # allometric scaling exponent standard deviation
  AHU.gm=1,	    # animal-human uncertainty geometric mean
  AHU.gsd=1.95,	    # animal-human uncertainty standard deviation
  OU.gm=1,	    # other uncertainties geometric mean
  OU.gsd=1,	    # other uncertainties geometric standard dev.
  sigmaH.gm=0.746,   # human variability sigmaH geometric mean
  sigmaH.gsd=1.5935, # other uncertainties geometric standard dev.
  n.samp=1e7	    # number of samples
) {
  # parms - resample model parameters
  if (n.samp > nrow(samp.parms)) {
    # with replacement if number of requested samples > number of parm samples
    resamp.parm.indx<-sample.int(nrow(samp.parms),n.samp,replace=TRUE);
  } else {
    # without replacement
    resamp.parm.indx<-sample.int(nrow(samp.parms),n.samp,replace=FALSE);
  }
  parms.samples <- samp.parms[resamp.parm.indx,];
  # DAF samples
  DAF.gm <- (bw.a/bw.h)^(1-alpha.m);
  DAF.gsd <- (bw.a/bw.h)^-alpha.sd;
  DAF.samples <- exp(rnorm(n.samp,m=log(DAF.gm),sd=log(DAF.gsd)));
  # AHU samples
  AHU.samples <- exp(rnorm(n.samp,m=log(AHU.gm),sd=log(AHU.gsd)));
  # OU samples
  OU.samples <- exp(rnorm(n.samp,m=log(OU.gm),sd=log(OU.gsd)));
  # sigmaH samples
  sigmaH.samples<-exp(rnorm(n.samp,mean=log(sigmaH.gm),sd=log(sigmaH.gsd)));
  # return list of samples
  list(parms=parms.samples[,-1],
       modelname=parms.samples[,1],
       DAF=DAF.samples,
       AHU=AHU.samples,
       OU=OU.samples,
       sigmaH=sigmaH.samples)
}
### Calculate overall population incidence
calc.pop.incidence.samples <- function(
  dose,                # Dose at which population incidence is calculated
  uncertainty.samples, # Uncertainty samples
  n.pop=1e4	     # Number of individuals in the population sample
) {
  # Number of uncertainty samples
  n.samp<-length(uncertainty.samples$sigmaH);
  # Generate population incidence for each sample
  pop.incidence.samples<-numeric(n.samp)
  for (j in 1:n.samp) {
    set.seed(314159) # use same random population each time
    dose.median <- dose * uncertainty.samples$AHU[j]*
      uncertainty.samples$OU[j]/uncertainty.samples$DAF[j];
    dose.distrib <- dose.median*exp(rnorm(n.pop)*uncertainty.samples$sigmaH[j]);
    pop.incidence.samples[j]<- mean(er.func(dose.distrib,
                                            as.numeric(uncertainty.samples$parms[j,]),
                                            uncertainty.samples$modelname[j]))
  }
  # Return samples of the population incidence at this dose
  pop.incidence.samples
}
#
### Calculate quantiles of the population incidence at multiple doses
get.pop.incidence <- function(
  D.vec,                # Doses at which to calculate
  uncertainty.samples,  # Uncertainty samples
  n.pop=1e4, 	      # Number of individual in the population sample
  prob=c(0.05,0.5,0.95) # Quantiles to report
) {
  pop.incidence.vec<-numeric(0);
  for (k in 1:length(D.vec)) {
    cat(k,"...\n");
    dose <- D.vec[k];
    pop.incidence.samples<-calc.pop.incidence.samples(dose,uncertainty.samples,
                                                      n.pop=n.pop);
    # Quantiles
    pop.incidence.quantiles<-as.data.frame(t(quantile(pop.incidence.samples,prob=prob,na.rm=TRUE)))
    # Etasq=Percent contribution to variance
    ## Take double log to regularize
    tmp <- data.frame(pop.incidence=log(-log(pop.incidence.samples)),
                      modelname=uncertainty.samples$modelname,
                      sigmaH=log(uncertainty.samples$sigmaH),
                      DAF=log(uncertainty.samples$DAF),
                      AHU=log(uncertainty.samples$AHU),
                      OU=log(uncertainty.samples$OU))
    tmp <- tmp[!is.infinite(tmp$pop.incidence),] # remove infinities
    if (length(unique(tmp$modelname))>1) {
      pop.incidence.aov<-aov(pop.incidence~modelname+sigmaH+DAF+AHU+OU,data=tmp)
      pop.incidence.etasq<-as.data.frame(t(etaSquared(pop.incidence.aov)[,"eta.sq"]))
    } else {
      pop.incidence.aov<-aov(pop.incidence~sigmaH+DAF+AHU+OU,data=tmp)
      pop.incidence.etasq<-cbind(data.frame(modelname=0),as.data.frame(t(etaSquared(pop.incidence.aov)[,"eta.sq"])))
    }
    # Combined
    pop.incidence.vec<-rbind(pop.incidence.vec,
                             cbind(pop.incidence.quantiles,pop.incidence.etasq))
  }
  # Return table of quantiles and proportion of variance
  cbind(data.frame(Dose=D.vec),as.data.frame(pop.incidence.vec))
}

