# Extra risk functions

#1. Quantal-linear model
get.er.ql <- function(a,b,dose) {
  er <- 1-exp(-b*dose)
  return(er)
}

#2. Logistic model
get.er.log<-function(a,b,dose){
  er<-1-exp(-b*dose)
  return(er)
}

#3. Probit model
get.er.pro<-function(a,b,dose){
  er<-(pnorm(a+b*dose) - pnorm(a))/(1 - pnorm(a))
  return(er)
}

#4. Weibull model
get.er.wei<-function(a,b,c,dose){
  er<-(1-exp(-c*dose^b))
  return(er)
}

#5. Multistage 2 model
get.er.ms2<-function(a,b,c,dose){
  er<-(1-exp(-b*dose-c*dose^2))
  return(er)
}

#6. Loglogistic model
get.er.llog<-function(a,b,c,dose){
  er<-1/(1+exp(-c-b*log(dose)))
  return(er)
}

#7. LogProbit model
get.er.lpro<-function(a,b,c,dose){
  er<-pnorm(c+b*log(dose))
  return(er)
}

#8. DichHill model
get.er.dh<-function(a,b,c,g,dose){
  resp<-a*g+(a-a*g)/(1+exp(-c-b*log(dose)))
  resp0<-a*g
  er<-(resp - resp0)/(1 - resp0)
  return(er)
}

#Combined

er.func <- function(dose,parms,modelname) {
  if (modelname == "Quantal-Linear") {
    er <- get.er.ql(parms[1],parms[2],dose)
  } else if (modelname == "Logistic") {
    er <- get.er.log(parms[1],parms[2],dose)
  } else if (modelname == "Probit") {
    er <- get.er.pro(parms[1],parms[2],dose)
  } else if (modelname == "Weibull") {
    er <- get.er.wei(parms[1],parms[2],parms[3],dose)
  } else if (modelname == "Multistage2") {
    er <- get.er.ms2(parms[1],parms[2],parms[3],dose)
  } else if (modelname == "LogLogistic") {
    er <- get.er.llog(parms[1],parms[2],parms[3],dose)
  } else if (modelname == "LogProbit") {
    er <- get.er.lpro(parms[1],parms[2],parms[3],dose)
  } else if (modelname == "DichHill") {
    er <- get.er.dh(parms[1],parms[2],parms[3],parms[4],dose)
  } else {
    er <- rep(NA,length(dose))
  }
  return(er)
}

