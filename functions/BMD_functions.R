################# 2. BMD Modeling Functions ################
#1. Quantal-linear model
get.resp.ql<-function(a,b,dose){
  resp<-a+(1-a)*(1-exp(-b*dose))
  return(resp)
}

get.bmd.ql.extra<-function(a,b,bmr){
  bmd<-log(1-bmr)/(-b)
  return(bmd)
}

#2. Logistic model
get.resp.log<-function(a,b,dose){
  resp<-1/(1+exp(-a-b*dose))
  return(resp)
}

get.bmd.log.extra<-function(a,b,bmr){
  bmd<-log((1-bmr)/(1+bmr*exp(-a)))/(-b)
  return(bmd)
}

#3. Probit model
get.resp.pro<-function(a,b,dose){
  resp<-pnorm(a+b*dose)
  return(resp)
}

get.bmd.pro.extra<-function(a,b,bmr){
  bmd<-(qnorm(bmr-bmr*pnorm(a)+pnorm(a))-a)/b
  return(bmd)
}

#4. Weibull model
get.resp.wei<-function(a,b,c,dose){
  resp<-a+(1-a)*(1-exp(-c*dose^b))
  return(resp)
}

get.bmd.wei.extra<-function(a,b,c,bmr){
  bmd<-exp(log(log((1-bmr*(1-a)-a)/(1-a))/(-c))/b)
  return(bmd)
}

#5. Multistage 2 model
get.resp.ms2<-function(a,b,c,dose){
  resp<-a+(1-a)*(1-exp(-b*dose-c*dose^2))
  return(resp)
}

get.bmd.ms2.extra<-function(a,b,c,bmr){
  if (c==0){
    bmd<-log(1-bmr)/(-b)
  } else {
    bmd<-(-b+sqrt(b^2-4*c*log(1-bmr)))/(2*c)   
  }
  return(bmd)
}

#6. Loglogistic model
get.resp.llog<-function(a,b,c,dose){
  resp<-a+(1-a)/(1+exp(-c-b*log(dose)))
  return(resp)
}

get.bmd.llog.extra<-function(a,b,c,bmr){
  r0<-get.resp.llog(a,b,c,0)
  bmd<-exp((log((1-a)/(bmr*(1-r0)+r0-a)-1)+c)/(-b))
  return(bmd)
}

#7. LogProbit model
get.resp.lpro<-function(a,b,c,dose){
  resp<-a+(1-a)*pnorm(c+b*log(dose))
  return(resp)
}

get.bmd.lpro.extra<-function(a,b,c,bmr){
  r0<-get.resp.lpro(a,b,c,0)
  bmd<-exp((qnorm((bmr*(1-r0)+r0-a)/(1-a))-c)/b)
  return(bmd)
}

#8. DichHill model
get.resp.dh<-function(a,b,c,g,dose){
  resp<-a*g+(a-a*g)/(1+exp(-c-b*log(dose)))
  return(resp)
}

get.bmd.dh.extra<-function(a,b,c,g,bmr){
  r0<-get.resp.dh(a,b,c,g,0)
  bmd<-exp((log((a-a*g)/(bmr*(1-r0)+r0-a*g)-1)+c)/(-b))
  return(bmd)
}

#Combined

bmd.func <- function(bmr,parms,modelname) {
  if (modelname == "Quantal-Linear") {
    bmd <- get.bmd.ql.extra(parms[1],parms[2],bmr)
  } else if (modelname == "Logistic") {
    bmd <- get.bmd.log.extra(parms[1],parms[2],bmr)
  } else if (modelname == "Probit") {
    bmd <- get.bmd.pro.extra(parms[1],parms[2],bmr)
  } else if (modelname == "Weibull") {
    bmd <- get.bmd.wei.extra(parms[1],parms[2],parms[3],bmr)
  } else if (modelname == "Multistage2") {
    bmd <- get.bmd.ms2.extra(parms[1],parms[2],parms[3],bmr)
  } else if (modelname == "LogLogistic") {
    bmd <- get.bmd.llog.extra(parms[1],parms[2],parms[3],bmr)
  } else if (modelname == "LogProbit") {
    bmd <- get.bmd.lpro.extra(parms[1],parms[2],parms[3],bmr)
  } else if (modelname == "DichHill") {
    bmd <- get.bmd.dh.extra(parms[1],parms[2],parms[3],parms[4],bmr)
  } else {
    bmd <- rep(NA,length(bmr))
  }
  return(bmd)
}

