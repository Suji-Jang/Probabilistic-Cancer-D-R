library(utils)
library(stats4)
library(rstan)
library(MASS)
library(plyr)
library(dplyr)
library(reshape2)
library(viridis)

############ 1. Stan Models ############

log_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
 }
 parameters {
  real log_a;
  real<lower=0> log_b;
 }
 model {
  log_a ~ uniform (-50,50);
  log_b ~ uniform (0,100);
  for (i in 1:len)
    y[i] ~ binomial(n[i],1/(1+exp(-log_a-log_b*d[i])));
 }
"

pro_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len]; 
  int<lower=0> n[len];
  real<lower=0> d[len];
 }
 parameters {
   real pro_a;
   real<lower=0> pro_b;
 }
 model {
  pro_a ~ uniform (-50,50);
  pro_b ~ uniform (0,100);
  for (i in 1:len)
    y[i] ~ binomial(n[i],normal_cdf(pro_a+pro_b*d[i],0,1));
 }
"

ql_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
 }
 parameters {
  real<lower=0,upper=1> ql_a;
  real<lower=0> ql_b;
 }
 model {
  ql_a ~ uniform (0,1);
  ql_b ~ uniform (0,100);
  for (i in 1:len)
    y[i] ~ binomial(n[i],ql_a+(1-ql_a)*(1-exp(-ql_b*d[i])));
 }
"

ms2_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
 }
 parameters {
  real<lower=0,upper=1> ms2_a;
  real<lower=0> ms2_b;
  real<lower=0> ms2_c;
 }
 model {
  ms2_a ~ uniform (0,1);
  ms2_b ~ uniform (0,100);
  ms2_c ~ uniform (0,100);
 for (i in 1:len)
  y[i] ~ binomial(n[i],ms2_a+(1-ms2_a)*(1-exp(-ms2_b*d[i]-ms2_c*(d[i]^2))));
 }
"

wei_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
  int rstrct;
 }
 parameters {
  real <lower=0,upper=1> wei_a;
  real <lower=rstrct> wei_b;
  real <lower=0> wei_c;
 }
 model {
  wei_a ~ uniform (0,1);
  wei_b ~ uniform (0,15);
  wei_c ~ uniform (0,50);
 for (i in 1:len)
   y[i] ~ binomial(n[i],wei_a+(1-wei_a)*(1-exp(-wei_c*(d[i])^wei_b)));
 }
"

llog_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
  int rstrct;
 }
 parameters {
  real <lower=0, upper=1> llog_a;
  real <lower=rstrct> llog_b;
  real llog_c;
 }
 model {
  llog_a ~ uniform (0,1);
  llog_b ~ uniform (0,15);
  llog_c ~ uniform (-5,15);
  for (i in 1:len)
    y[i] ~ binomial(n[i],llog_a+(1-llog_a)/(1+exp(-llog_c-llog_b*log(d[i]))));
 }
"

lpro_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
  int rstrct;
 }
 parameters {
  real <lower=0, upper=1> lpro_a;
  real <lower=rstrct> lpro_b;
  real lpro_c;
 }
 model {
  lpro_a ~ uniform (0,1);
  lpro_b ~ uniform (0,15);
  lpro_c ~ uniform (-5,15);
 for (i in 1:len)
   y[i] ~ binomial(n[i],lpro_a+(1-lpro_a)*normal_cdf(lpro_c+lpro_b*log(d[i]),0,1));
 }
"

dh_model <- "
 data {
  int<lower=0> len;
  int<lower=0> y[len];
  int<lower=0> n[len];
  real<lower=0> d[len];
  int rstrct; 
 }
 parameters {
  real <lower=0, upper=1> dh_a;
  real <lower=rstrct> dh_b;
  real dh_c;
  real <lower=0, upper=1> dh_g;
 }
 model {
  dh_a ~ uniform (0,1);
  dh_b ~ uniform (0,15);
  dh_c ~ uniform (-5,15);
  dh_g ~ uniform (0,1);
  for (i in 1:len)
    y[i] ~ binomial(n[i],dh_a*dh_g+(dh_a-dh_a*dh_g)/(1+exp(-dh_c-dh_b*log(d[i]))));
 }
"
################# End of Stan Models ###############

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

############### 3. Other functions and constants ###############
zeroish<-1e-8  # Globle Variable

get.likelihood.dich<-function(pp,yy,nn){
  pp[pp==0]<-zeroish
  pp[pp==1]<-1-zeroish
  out<-sum(yy*log(pp)+(nn-yy)*log(1-pp))
  return(out)
}

get.para.stats<-function(x){
  out<-rep(NA,5)
  out[1]<-mean(x,na.rm = TRUE)
  out[2]<-sd(x,na.rm = TRUE)
  out[3]<-quantile(x,0.025,na.rm = TRUE)
  out[4]<-quantile(x,0.5,na.rm = TRUE)
  out[5]<-quantile(x,0.975,na.rm = TRUE)
  return(out)
}

get.BMD.stats<-function(x){
  out<-rep(NA,3)
  out[1]<-quantile(x,0.5,na.rm = TRUE)    #BMD
  out[2]<-quantile(x,0.05,na.rm = TRUE)   #BMDL
  out[3]<-quantile(x,0.95,na.rm = TRUE)   #BMDU
  return(out)
}
############# End of other functions and constants ##############

############# 4. Data loading, output preparation, and parameter settings ###############
# Load data sets
setwd("C:/Users/lover/OneDrive - Texas A&M University/Projects/CTV")
datasets_dich<-read.csv("BMD_data_092022.csv",header=TRUE)
study.index<-unique(datasets_dich[,1])
num.study<-length(study.index)
study.index.all<-datasets_dich[,1]
dose.all<-datasets_dich[,2]
subjnum.all<-datasets_dich[,3]
casenum.all<-datasets_dich[,4]

# Prepare Output Matrix #
output.mtx<-matrix(NA,num.study,70)

#output.mtx.ql<-matrix(NA,num.study,35)
#output.mtx.log<-matrix(NA,num.study,35)
#output.mtx.pro<-matrix(NA,num.study,35)
#output.mtx.wei<-matrix(NA,num.study,43)
#output.mtx.ms2<-matrix(NA,num.study,43)
#output.mtx.llog<-matrix(NA,num.study,43)
#output.mtx.lpro<-matrix(NA,num.study,43)
#output.mtx.dh<-matrix(NA,num.study,51)

length.mcmc.sample<-5000  # This value determines how many posterior samples to be printed in output files
# This value cannot exceed mcmc_iter_unit - mcmc_output_len, i.e., 15000 in this default setting

output.par.ql<-matrix(NA,length.mcmc.sample,num.study*2)
output.par.log<-matrix(NA,length.mcmc.sample,num.study*2)
output.par.pro<-matrix(NA,length.mcmc.sample,num.study*2)
output.par.wei<-matrix(NA,length.mcmc.sample,num.study*3)
output.par.ms2<-matrix(NA,length.mcmc.sample,num.study*3)
output.par.llog<-matrix(NA,length.mcmc.sample,num.study*3)
output.par.lpro<-matrix(NA,length.mcmc.sample,num.study*3)
output.par.dh<-matrix(NA,length.mcmc.sample,num.study*4)

# Set some parameters
set.seed(1)
seed.vec<-round(runif(num.study,0,100000))

mcmc_iter_unit = 30000
mcmc_output_len = 15000
mcmc_chain_num = 1
rstrct = 1  #restrict value
############# End of Data loading, output preparation, and parameter settings ###############


for (s in 1:num.study){
  
  mcmc_seed=seed.vec[s]
  start.time<-Sys.time()
  dose.temp<-dose.all[study.index.all==study.index[s]]
  nsub<-subjnum.all[study.index.all==study.index[s]]
  ncase<-casenum.all[study.index.all==study.index[s]]
  dose.max<-max(dose.temp)
  dose<-dose.temp/dose.max
  dose.no0<-dose
  dose.no0[dose.no0==0]<-zeroish/dose.max
  group.num<-length(dose)
  
  # Input data for Logistic model 
  data_log <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose)
  
  # Input data for Probit Model
  data_pro <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose)
  
  # Input data for ql model
  data_ql <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose)
  
  # Input data for ms2 model 
  data_ms2 <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose)
  
  # Input data for Weibull model 
  data_wei <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose,
    rstrct=rstrct)
  
  # Input data for Loglogistic model 
  data_llog <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose.no0,
    rstrct=rstrct)
  
  # Input data for Logprobit model 
  data_lpro = list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose.no0,
    rstrct=rstrct)
  
  # Input data for Dichotomous Hill model 
  data_dh <- list(
    len=group.num,
    y=ncase,
    n=nsub,
    d=dose.no0,
    rstrct=rstrct)
  
  
  mcmc_iter.ql =mcmc_iter_unit
  mcmc_warmup.ql = mcmc_iter.ql - mcmc_output_len
  fit_ql = stan(model_code=ql_model,data=data_ql,iter=mcmc_iter.ql,warmup=mcmc_warmup.ql,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.log =mcmc_iter_unit
  mcmc_warmup.log = mcmc_iter.log - mcmc_output_len
  fit_log = stan(model_code=log_model,data=data_log,iter=mcmc_iter.log,warmup=mcmc_warmup.log,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.pro =mcmc_iter_unit
  mcmc_warmup.pro = mcmc_iter.pro - mcmc_output_len
  fit_pro = stan(model_code=pro_model,data=data_pro,iter=mcmc_iter.pro,warmup=mcmc_warmup.pro,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.wei =mcmc_iter_unit
  mcmc_warmup.wei = mcmc_iter.wei - mcmc_output_len
  fit_wei = stan(model_code=wei_model,data=data_wei,iter=mcmc_iter.wei,warmup=mcmc_warmup.wei,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.ms2 =mcmc_iter_unit
  mcmc_warmup.ms2 = mcmc_iter.ms2 - mcmc_output_len  
  fit_ms2 = stan(model_code=ms2_model,data=data_ms2,iter=mcmc_iter.ms2,warmup=mcmc_warmup.ms2,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.llog =mcmc_iter_unit
  mcmc_warmup.llog = mcmc_iter.llog - mcmc_output_len      
  fit_llog = stan(model_code=llog_model,data=data_llog,iter=mcmc_iter.llog,warmup=mcmc_warmup.llog,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.lpro =mcmc_iter_unit
  mcmc_warmup.lpro = mcmc_iter.lpro - mcmc_output_len  
  fit_lpro = stan(model_code=lpro_model,data=data_lpro,iter=mcmc_iter.lpro,warmup=mcmc_warmup.lpro,chains=mcmc_chain_num,seed=mcmc_seed)
  
  mcmc_iter.dh =mcmc_iter_unit
  mcmc_warmup.dh = mcmc_iter.dh - mcmc_output_len   
  fit_dh = stan(model_code=dh_model,data=data_dh,iter=mcmc_iter.dh,warmup=mcmc_warmup.dh,chains=mcmc_chain_num,seed=mcmc_seed)
  
  ql_a<-as.array(fit_ql)[1:mcmc_output_len]
  ql_b<-as.array(fit_ql)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  
  log_a<-as.array(fit_log)[1:mcmc_output_len]
  log_b<-as.array(fit_log)[(mcmc_output_len+1):(mcmc_output_len*2)]
  
  pro_a<-as.array(fit_pro)[1:mcmc_output_len]
  pro_b<-as.array(fit_pro)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  
  wei_a<-as.array(fit_wei)[1:mcmc_output_len]
  wei_b<-as.array(fit_wei)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  wei_c<-as.array(fit_wei)[(mcmc_output_len*2+1):(mcmc_output_len*3)] 
  
  ms2_a<-as.array(fit_ms2)[1:mcmc_output_len]
  ms2_b<-as.array(fit_ms2)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  ms2_c<-as.array(fit_ms2)[(mcmc_output_len*2+1):(mcmc_output_len*3)] 
  
  llog_a<-as.array(fit_llog)[1:mcmc_output_len]
  llog_b<-as.array(fit_llog)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  llog_c<-as.array(fit_llog)[(mcmc_output_len*2+1):(mcmc_output_len*3)] 
  
  lpro_a<-as.array(fit_lpro)[1:mcmc_output_len]
  lpro_b<-as.array(fit_lpro)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  lpro_c<-as.array(fit_lpro)[(mcmc_output_len*2+1):(mcmc_output_len*3)] 
  
  dh_a<-as.array(fit_dh)[1:mcmc_output_len]
  dh_b<-as.array(fit_dh)[(mcmc_output_len+1):(mcmc_output_len*2)] 
  dh_c<-as.array(fit_dh)[(mcmc_output_len*2+1):(mcmc_output_len*3)] 
  dh_g<-as.array(fit_dh)[(mcmc_output_len*3+1):(mcmc_output_len*4)] 
  
  ## put the MCMC sample in output matrix ##
  mcmcsample.out.index<-(mcmc_output_len-length.mcmc.sample+1):mcmc_output_len
  output.par.ql[,(s-1)*2+1]<-ql_a[mcmcsample.out.index]
  output.par.ql[,(s-1)*2+2]<-ql_b[mcmcsample.out.index]
  output.par.log[,(s-1)*2+1]<-log_a[mcmcsample.out.index]
  output.par.log[,(s-1)*2+2]<-log_b[mcmcsample.out.index]
  output.par.pro[,(s-1)*2+1]<-pro_a[mcmcsample.out.index]
  output.par.pro[,(s-1)*2+2]<-pro_b[mcmcsample.out.index]
  output.par.wei[,(s-1)*3+1]<-wei_a[mcmcsample.out.index]
  output.par.wei[,(s-1)*3+2]<-wei_b[mcmcsample.out.index]
  output.par.wei[,(s-1)*3+3]<-wei_c[mcmcsample.out.index]
  output.par.ms2[,(s-1)*3+1]<-ms2_a[mcmcsample.out.index]
  output.par.ms2[,(s-1)*3+2]<-ms2_b[mcmcsample.out.index]
  output.par.ms2[,(s-1)*3+3]<-ms2_c[mcmcsample.out.index]
  output.par.llog[,(s-1)*3+1]<-llog_a[mcmcsample.out.index]
  output.par.llog[,(s-1)*3+2]<-llog_b[mcmcsample.out.index]
  output.par.llog[,(s-1)*3+3]<-llog_c[mcmcsample.out.index]
  output.par.lpro[,(s-1)*3+1]<-lpro_a[mcmcsample.out.index]
  output.par.lpro[,(s-1)*3+2]<-lpro_b[mcmcsample.out.index]
  output.par.lpro[,(s-1)*3+3]<-lpro_c[mcmcsample.out.index]
  output.par.dh[,(s-1)*4+1]<-dh_a[mcmcsample.out.index]
  output.par.dh[,(s-1)*4+2]<-dh_b[mcmcsample.out.index]
  output.par.dh[,(s-1)*4+3]<-dh_c[mcmcsample.out.index]
  output.par.dh[,(s-1)*4+4]<-dh_g[mcmcsample.out.index]
  
  test_stat_pred_log = rep(NA,mcmc_output_len)
  test_stat_obs_log = rep(NA,mcmc_output_len)
  
  test_stat_pred_pro = rep(NA,mcmc_output_len)
  test_stat_obs_pro = rep(NA,mcmc_output_len) 
  
  test_stat_pred_ql = rep(NA,mcmc_output_len)
  test_stat_obs_ql = rep(NA,mcmc_output_len) 
  
  test_stat_pred_ms2 = rep(NA,mcmc_output_len)
  test_stat_obs_ms2 = rep(NA,mcmc_output_len) 
  
  test_stat_pred_wei = rep(NA,mcmc_output_len)
  test_stat_obs_wei = rep(NA,mcmc_output_len) 
  
  test_stat_pred_llog = rep(NA,mcmc_output_len)
  test_stat_obs_llog = rep(NA,mcmc_output_len) 
  
  test_stat_pred_lpro = rep(NA,mcmc_output_len)
  test_stat_obs_lpro = rep(NA,mcmc_output_len) 
  
  test_stat_pred_dh = rep(NA,mcmc_output_len)
  test_stat_obs_dh = rep(NA,mcmc_output_len)
  
  weight_matrix = matrix(NA,mcmc_output_len,8)  #nine models, now eight
  
  
  for (i in 1:mcmc_output_len){
    
    p_pred_ql<-get.resp.ql(ql_a[i],ql_b[i],dose)
    p_pred_ql[p_pred_ql==0]=zeroish 
    p_pred_ql[p_pred_ql==1]=1-zeroish
    y_pred_ql = rbinom(group.num,nsub,p_pred_ql)
    test_stat_pred_ql[i]=-2*get.likelihood.dich(p_pred_ql,y_pred_ql,nsub)
    test_stat_obs_ql[i]=-2*get.likelihood.dich(p_pred_ql,ncase,nsub)
    
    p_pred_log<-get.resp.log(log_a[i],log_b[i],dose)
    p_pred_log[p_pred_log==0]=zeroish 
    p_pred_log[p_pred_log==1]=1-zeroish
    y_pred_log = rbinom(group.num,nsub,p_pred_log)
    test_stat_pred_log[i]=-2*get.likelihood.dich(p_pred_log,y_pred_log,nsub)
    test_stat_obs_log[i]=-2*get.likelihood.dich(p_pred_log,ncase,nsub)
    
    p_pred_pro<-get.resp.pro(pro_a[i],pro_b[i],dose)
    p_pred_pro[p_pred_pro==0]=zeroish 
    p_pred_pro[p_pred_pro==1]=1-zeroish
    y_pred_pro = rbinom(group.num,nsub,p_pred_pro)
    test_stat_pred_pro[i]=-2*get.likelihood.dich(p_pred_pro,y_pred_pro,nsub)
    test_stat_obs_pro[i]=-2*get.likelihood.dich(p_pred_pro,ncase,nsub)
    
    p_pred_wei<-get.resp.wei(wei_a[i],wei_b[i],wei_c[i],dose)
    p_pred_wei[p_pred_wei==0]=zeroish 
    p_pred_wei[p_pred_wei==1]=1-zeroish
    y_pred_wei = rbinom(group.num,nsub,p_pred_wei)
    test_stat_pred_wei[i]=-2*get.likelihood.dich(p_pred_wei,y_pred_wei,nsub)
    test_stat_obs_wei[i]=-2*get.likelihood.dich(p_pred_wei,ncase,nsub) 
    
    p_pred_ms2<-get.resp.ms2(ms2_a[i],ms2_b[i],ms2_c[i],dose)
    p_pred_ms2[p_pred_ms2==0]=zeroish 
    p_pred_ms2[p_pred_ms2==1]=1-zeroish
    y_pred_ms2 = rbinom(group.num,nsub,p_pred_ms2)
    test_stat_pred_ms2[i]=-2*get.likelihood.dich(p_pred_ms2,y_pred_ms2,nsub)
    test_stat_obs_ms2[i]=-2*get.likelihood.dich(p_pred_ms2,ncase,nsub)
    
    p_pred_llog<-get.resp.llog(llog_a[i],llog_b[i],llog_c[i],dose)
    p_pred_llog[p_pred_llog==0]=zeroish 
    p_pred_llog[p_pred_llog==1]=1-zeroish
    y_pred_llog = rbinom(group.num,nsub,p_pred_llog)
    test_stat_pred_llog[i]=-2*get.likelihood.dich(p_pred_llog,y_pred_llog,nsub)
    test_stat_obs_llog[i]=-2*get.likelihood.dich(p_pred_llog,ncase,nsub)
    
    p_pred_lpro<-get.resp.lpro(lpro_a[i],lpro_b[i],lpro_c[i],dose)
    p_pred_lpro[p_pred_lpro==0]=zeroish 
    p_pred_lpro[p_pred_lpro==1]=1-zeroish
    y_pred_lpro = rbinom(group.num,nsub,p_pred_lpro)
    test_stat_pred_lpro[i]=-2*get.likelihood.dich(p_pred_lpro,y_pred_lpro,nsub)
    test_stat_obs_lpro[i]=-2*get.likelihood.dich(p_pred_lpro,ncase,nsub)
    
    p_pred_dh<-get.resp.dh(dh_a[i],dh_b[i],dh_c[i],dh_g[i],dose)
    p_pred_dh[p_pred_dh==0]=zeroish 
    p_pred_dh[p_pred_dh==1]=1-zeroish
    y_pred_dh = rbinom(group.num,nsub,p_pred_dh)
    test_stat_pred_dh[i]=-2*get.likelihood.dich(p_pred_dh,y_pred_dh,nsub)
    test_stat_obs_dh[i]=-2*get.likelihood.dich(p_pred_dh,ncase,nsub)
    
    m_vector<-c(
      -0.5*test_stat_obs_ql[i]-2*log(group.num)/2, -0.5*test_stat_obs_log[i]-2*log(group.num)/2,-0.5*test_stat_obs_pro[i]-2*log(group.num)/2,
      -0.5*test_stat_obs_wei[i]-3*log(group.num)/2,-0.5*test_stat_obs_ms2[i]-3*log(group.num)/2,-0.5*test_stat_obs_llog[i]-3*log(group.num)/2,
      -0.5*test_stat_obs_lpro[i]-3*log(group.num)/2,-0.5*test_stat_obs_dh[i]-4*log(group.num)/2)
    m_vector<-exp(m_vector-min(m_vector))
    
    weight_matrix[i,]<-m_vector/sum(m_vector)
    
  }
  
  test_stat_ql = test_stat_pred_ql - test_stat_obs_ql
  test_stat_ql[test_stat_ql<0] = NA
  ppp_ql = sum(is.na(test_stat_ql))/length(test_stat_ql)
  
  test_stat_log = test_stat_pred_log - test_stat_obs_log
  test_stat_log[test_stat_log<0] = NA
  ppp_log = sum(is.na(test_stat_log))/length(test_stat_log)
  
  test_stat_pro = test_stat_pred_pro - test_stat_obs_pro
  test_stat_pro[test_stat_pro<0] = NA
  ppp_pro = sum(is.na(test_stat_pro))/length(test_stat_pro)
  
  test_stat_wei = test_stat_pred_wei - test_stat_obs_wei
  test_stat_wei[test_stat_wei<0] = NA
  ppp_wei = sum(is.na(test_stat_wei))/length(test_stat_wei)
  
  test_stat_ms2 = test_stat_pred_ms2 - test_stat_obs_ms2
  test_stat_ms2[test_stat_ms2<0] = NA
  ppp_ms2 = sum(is.na(test_stat_ms2))/length(test_stat_ms2)
  
  test_stat_llog = test_stat_pred_llog - test_stat_obs_llog
  test_stat_llog[test_stat_llog<0] = NA
  ppp_llog = sum(is.na(test_stat_llog))/length(test_stat_llog)
  
  test_stat_lpro = test_stat_pred_lpro - test_stat_obs_lpro
  test_stat_lpro[test_stat_lpro<0] = NA
  ppp_lpro = sum(is.na(test_stat_lpro))/length(test_stat_lpro)
  
  test_stat_dh = test_stat_pred_dh - test_stat_obs_dh
  test_stat_dh[test_stat_dh<0] = NA
  ppp_dh = sum(is.na(test_stat_dh))/length(test_stat_dh)
  
  weight_mean<-colMeans(weight_matrix,na.rm=TRUE)
  
  #prepare output matrix#
  # Quantile-Linear
  bmd.extra.ql.10<-get.bmd.ql.extra(ql_a,ql_b,0.1)*dose.max
  bmd.extra.ql.1<-get.bmd.ql.extra(ql_a,ql_b,0.01)*dose.max
  output.mtx[s,1]<-ppp_ql
  output.mtx[s,2]<-weight_mean[1]
  output.mtx[s,3:5]<-get.BMD.stats(bmd.extra.ql.10)
  output.mtx[s,6:8]<-get.BMD.stats(bmd.extra.ql.1)
  
  # Logistic
  bmd.extra.log.10<-get.bmd.log.extra(log_a,log_b,0.1)*dose.max
  bmd.extra.log.1<-get.bmd.log.extra(log_a,log_b,0.01)*dose.max
  output.mtx[s,9]<-ppp_log
  output.mtx[s,10]<-weight_mean[2]
  output.mtx[s,11:13]<-get.BMD.stats(bmd.extra.log.10)
  output.mtx[s,14:16]<-get.BMD.stats(bmd.extra.log.1)
  
  
  # Probit
  bmd.extra.pro.10<-get.bmd.pro.extra(pro_a,pro_b,0.1)*dose.max
  bmd.extra.pro.1<-get.bmd.pro.extra(pro_a,pro_b,0.01)*dose.max
  output.mtx[s,17]<-ppp_pro
  output.mtx[s,18]<-weight_mean[3]
  output.mtx[s,19:21]<-get.BMD.stats(bmd.extra.pro.10)
  output.mtx[s,22:24]<-get.BMD.stats(bmd.extra.pro.1)
  
  
  # Weibull
  bmd.extra.wei.10<-get.bmd.wei.extra(wei_a,wei_b,wei_c,0.1)*dose.max
  bmd.extra.wei.1<-get.bmd.wei.extra(wei_a,wei_b,wei_c,0.01)*dose.max
  output.mtx[s,25]<-ppp_wei
  output.mtx[s,26]<-weight_mean[4]
  output.mtx[s,27:29]<-get.BMD.stats(bmd.extra.wei.10)
  output.mtx[s,30:32]<-get.BMD.stats(bmd.extra.wei.1)
  
  
  # Multistage 2
  bmd.extra.ms2.10<-get.bmd.ms2.extra(ms2_a,ms2_b,ms2_c,0.1)*dose.max
  bmd.extra.ms2.1<-get.bmd.ms2.extra(ms2_a,ms2_b,ms2_c,0.01)*dose.max
  output.mtx[s,33]<-ppp_ms2
  output.mtx[s,34]<-weight_mean[5]
  output.mtx[s,35:37]<-get.BMD.stats(bmd.extra.ms2.10)
  output.mtx[s,38:40]<-get.BMD.stats(bmd.extra.ms2.1)
  
  
  # Loglogistic
  bmd.extra.llog.10<-get.bmd.llog.extra(llog_a,llog_b,llog_c,0.1)*dose.max
  bmd.extra.llog.1<-get.bmd.llog.extra(llog_a,llog_b,llog_c,0.01)*dose.max
  output.mtx[s,41]<-ppp_llog
  output.mtx[s,42]<-weight_mean[6]
  output.mtx[s,43:45]<-get.BMD.stats(bmd.extra.llog.10)
  output.mtx[s,46:48]<-get.BMD.stats(bmd.extra.llog.1)
  
  
  # LogProbit
  bmd.extra.lpro.10<-get.bmd.lpro.extra(lpro_a,lpro_b,lpro_c,0.1)*dose.max
  bmd.extra.lpro.1<-get.bmd.lpro.extra(lpro_a,lpro_b,lpro_c,0.01)*dose.max
  output.mtx[s,49]<-ppp_lpro
  output.mtx[s,50]<-weight_mean[7]
  output.mtx[s,51:53]<-get.BMD.stats(bmd.extra.lpro.10)
  output.mtx[s,54:56]<-get.BMD.stats(bmd.extra.lpro.1)
  
  
  # Dichotomous Hill
  bmd.extra.dh.10<-get.bmd.dh.extra(dh_a,dh_b,dh_c,dh_g,0.1)*dose.max
  bmd.extra.dh.1<-get.bmd.dh.extra(dh_a,dh_b,dh_c,dh_g,0.01)*dose.max
  output.mtx[s,57]<-ppp_dh
  output.mtx[s,58]<-weight_mean[8]
  output.mtx[s,59:61]<-get.BMD.stats(bmd.extra.dh.10)
  output.mtx[s,62:64]<-get.BMD.stats(bmd.extra.dh.1)
  
  
  #Calculate model averaged BMD
  samplesize.vector<-round(weight_mean*mcmc_output_len)
  
  if (samplesize.vector[1]!=0){
    bmd.ext10.ql<-sample(bmd.extra.ql.10,samplesize.vector[1])
    bmd.ext1.ql<-sample(bmd.extra.ql.1,samplesize.vector[1])
  } else {
    bmd.ext10.ql<-vector(mode="numeric", length=0)
    bmd.ext1.ql<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[2]!=0){
    bmd.ext10.log<-sample(bmd.extra.log.10,samplesize.vector[2])
    bmd.ext1.log<-sample(bmd.extra.log.1,samplesize.vector[2])
  } else {
    bmd.ext10.log<-vector(mode="numeric", length=0)
    bmd.ext1.log<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[3]!=0){
    bmd.ext10.pro<-sample(bmd.extra.pro.10,samplesize.vector[3])
    bmd.ext1.pro<-sample(bmd.extra.pro.1,samplesize.vector[3])
  } else {
    bmd.ext10.pro<-vector(mode="numeric", length=0)
    bmd.ext1.pro<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[4]!=0){
    bmd.ext10.wei<-sample(bmd.extra.wei.10,samplesize.vector[4])
    bmd.ext1.wei<-sample(bmd.extra.wei.1,samplesize.vector[4])
  } else {
    bmd.ext10.wei<-vector(mode="numeric", length=0)
    bmd.ext1.wei<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[5]!=0){
    bmd.ext10.ms2<-sample(bmd.extra.ms2.10,samplesize.vector[5])
    bmd.ext1.ms2<-sample(bmd.extra.ms2.1,samplesize.vector[5])
  } else {
    bmd.ext10.ms2<-vector(mode="numeric", length=0)
    bmd.ext1.ms2<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[6]!=0){
    bmd.ext10.llog<-sample(bmd.extra.llog.10,samplesize.vector[6])
    bmd.ext1.llog<-sample(bmd.extra.llog.1,samplesize.vector[6])
  } else {
    bmd.ext10.llog<-vector(mode="numeric", length=0)
    bmd.ext1.llog<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[7]!=0){
    bmd.ext10.lpro<-sample(bmd.extra.lpro.10,samplesize.vector[7])
    bmd.ext1.lpro<-sample(bmd.extra.lpro.1,samplesize.vector[7])
  } else {
    bmd.ext10.lpro<-vector(mode="numeric", length=0)
    bmd.ext1.lpro<-vector(mode="numeric", length=0)
  }
  
  if (samplesize.vector[8]!=0){
    bmd.ext10.dh<-sample(bmd.extra.dh.10,samplesize.vector[8])
    bmd.ext1.dh<-sample(bmd.extra.dh.1,samplesize.vector[8])
  } else {
    bmd.ext10.dh<-vector(mode="numeric", length=0)
    bmd.ext1.dh<-vector(mode="numeric", length=0)
  }
  
  bmd.extra.10.ma<-c(bmd.ext10.ql,bmd.ext10.log,bmd.ext10.pro,bmd.ext10.wei,bmd.ext10.ms2,bmd.ext10.llog,bmd.ext10.lpro,bmd.ext10.dh)
  bmd.extra.1.ma<-c(bmd.ext1.ql,bmd.ext1.log,bmd.ext1.pro,bmd.ext1.wei,bmd.ext1.ms2,bmd.ext1.llog,bmd.ext1.lpro,bmd.ext1.dh)
  
  output.mtx[s,65:67]<-get.BMD.stats(bmd.extra.10.ma)
  output.mtx[s,68:70]<-get.BMD.stats(bmd.extra.1.ma)
  
  
  end.time<-Sys.time()
  
  print(s)
}  

rname<-study.index[1:num.study]
cname1<-c("ql-ppp","ql-weight","ql-ext10","ql-ext10L","ql-ext10U","ql-ext1","ql-ext1L","ql-ext1U",
          "log-ppp","log-weight","log-ext10","log-ext10L","log-ext10U","log-ext1","log-ext1L","log-ext1U",
          "pro-ppp","pro-weight","pro-ext10","pro-ext10L","pro-ext10U","pro-ext1","pro-ext1L","pro-ext1U",
          "wei-ppp","wei-weight","wei-ext10","wei-ext10L","wei-ext10U","wei-ext1","wei-ext1L","wei-ext1U", 
          "ms2-ppp","ms2-weight","ms2-ext10","ms2-ext10L","ms2-ext10U","ms2-ext1","ms2-ext1L","ms2-ext1U",
          "llog-ppp","llog-weight","llog-ext10","llog-ext10L","llog-ext10U","llog-ext1","llog-ext1L","llog-ext1U", 
          "lpro-ppp","lpro-weight","lpro-ext10","lpro-ext10L","lpro-ext10U","lpro-ext1","lpro-ext1L","lpro-ext1U",
          "dh-ppp","dh-weight","dh-ext10","dh-ext10L","dh-ext10U","dh-ext1","dh-ext1L","dh-ext1U",
          "MA-ext10","MA-ext10L","MA-ext10U","MA-ext1","MA-ext1L","MA-ext1U")


rownames(output.mtx)<-rname
colnames(output.mtx)<-cname1

write.csv(output.mtx,"Summary_Output_data_092022.csv")
###################################################################################### From here!

write.csv(output.par.ql,"MCMC_Sample_ql_data_092022.csv")
write.csv(output.par.log,"MCMC_Sample_log_data_092022.csv")
write.csv(output.par.pro,"MCMC_Sample_pro_data_092022.csv")
write.csv(output.par.wei,"MCMC_Sample_wei_data_092022.csv")
write.csv(output.par.ms2,"MCMC_Sample_ms2_data_092022.csv")
write.csv(output.par.llog,"MCMC_Sample_llog_data_092022.csv")
write.csv(output.par.lpro,"MCMC_Sample_lpro_data_092022.csv")
write.csv(output.par.dh,"MCMC_Sample_dh_data_092022.csv")

################# Parameter Samples ######################################

################# When running some samples ##########################
setwd("C:/Users/lover/Documents/Projects/CTV")

output.par.ql <- read.csv("MCMC_Sample_ql_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.log <- read.csv("MCMC_Sample_log_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.pro <- read.csv("MCMC_Sample_pro_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.wei <- read.csv("MCMC_Sample_wei_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.ms2 <- read.csv("MCMC_Sample_ms2_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.llog <- read.csv("MCMC_Sample_llog_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.lpro <- read.csv("MCMC_Sample_lpro_data_092022.csv",as.is=TRUE,row.names = 1)
output.par.dh <- read.csv("MCMC_Sample_dh_data_092022.csv",as.is=TRUE,row.names = 1)

mod.summary <- read.csv("Summary_Output_data_092022.csv",as.is=TRUE)

mod.summary <- mod.summary[,c("ql.weight","log.weight","pro.weight","wei.weight","ms2.weight",
                              "llog.weight","lpro.weight","dh.weight")]

#############################################################################
mod.summary <- round(mod.summary * 5000, digits=0)
write.csv(mod.summary,"Model_weight_092022.csv") 

ql.list <- mod.summary$ql.weight
log.list <- mod.summary$log.weight
pro.list <- mod.summary$pro.weight
wei.list <- mod.summary$wei.weight
ms2.list <- mod.summary$ms2.weight
llog.list <- mod.summary$llog.weight
lpro.list <- mod.summary$lpro.weight
dh.list <- mod.summary$dh.weight

ql.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(ql.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(ql.list)) {
  if (ql.list[i]!=0){
    rand.par <- output.par.ql[sample(nrow(output.par.ql),size=ql.list[i],replace=FALSE),(2*i-1):(2*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2")
    rand.par$Dataset <- i
    rand.par$Model <- "Quantal-Linear"
    rand.par$Par3 <- NA
    rand.par$Par4 <- NA
    ql.parms <- rbind(ql.parms,rand.par)
  }
  else next
}
sum(mod.summary$ql.weight)

log.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(log.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(log.list)) {
  if (log.list[i]!=0){
    rand.par <- output.par.log[sample(nrow(output.par.log),size=log.list[i],replace=FALSE),(2*i-1):(2*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2")
    rand.par$Dataset <- i
    rand.par$Model <- "Logistic"
    rand.par$Par3 <- NA
    rand.par$Par4 <- NA
    log.parms <- rbind(log.parms,rand.par)
  }
  else next
}
sum(mod.summary$log.weight)

pro.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(pro.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(pro.list)) {
  if (pro.list[i]!=0){
    rand.par <- output.par.pro[sample(nrow(output.par.pro),size=pro.list[i],replace=FALSE),(2*i-1):(2*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2")
    rand.par$Dataset <- i
    rand.par$Model <- "Probit"
    rand.par$Par3 <- NA
    rand.par$Par4 <- NA
    pro.parms <- rbind(pro.parms,rand.par)
  }
  else next
}
sum(mod.summary$pro.weight)

wei.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(wei.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(wei.list)) {
  if (wei.list[i]!=0){
    rand.par <- output.par.wei[sample(nrow(output.par.wei),size=wei.list[i],replace=FALSE),(3*i-2):(3*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2","Par3")
    rand.par$Dataset <- i
    rand.par$Model <- "Weibull"
    rand.par$Par4 <- NA
    wei.parms <- rbind(wei.parms,rand.par)
  }
  else next
}
sum(mod.summary$wei.weight)

ms2.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(ms2.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(ms2.list)) {
  if (ms2.list[i]!=0){
    rand.par <- output.par.ms2[sample(nrow(output.par.ms2),size=ms2.list[i],replace=FALSE),(3*i-2):(3*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2","Par3")
    rand.par$Dataset <- i
    rand.par$Model <- "Multistage2"
    rand.par$Par4 <- NA
    ms2.parms <- rbind(ms2.parms,rand.par)
  }
  else next
}
sum(mod.summary$ms2.weight)

llog.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(llog.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(llog.list)) {
  if (llog.list[i]!=0){
    rand.par <- output.par.llog[sample(nrow(output.par.llog),size=llog.list[i],replace=FALSE),(3*i-2):(3*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2","Par3")
    rand.par$Dataset <- i
    rand.par$Model <- "LogLogistic"
    rand.par$Par4 <- NA
    llog.parms <- rbind(llog.parms,rand.par)
  }
  else next
}
sum(mod.summary$llog.weight)

lpro.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(lpro.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(lpro.list)) {
  if (lpro.list[i]!=0){
    rand.par <- output.par.lpro[sample(nrow(output.par.lpro),size=lpro.list[i],replace=FALSE),(3*i-2):(3*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2","Par3")
    rand.par$Dataset <- i
    rand.par$Model <- "LogProbit"
    rand.par$Par4 <- NA
    lpro.parms <- rbind(lpro.parms,rand.par)
  }
  else next
}
sum(mod.summary$lpro.weight)

dh.parms <- as.data.frame(matrix(NA,nrow=0,ncol=6))
colnames(dh.parms) <- c("Dataset","Model","Par1","Par2","Par3","Par4")
for (i in 1:length(dh.list)) {
  if (dh.list[i]!=0){
    rand.par <- output.par.dh[sample(nrow(output.par.dh),size=dh.list[i],replace=FALSE),(4*i-3):(4*i)]
    rand.par <- as.data.frame(rand.par)
    colnames(rand.par) <- c("Par1","Par2","Par3","Par4")
    rand.par$Dataset <- i
    rand.par$Model <- "DichHill"
    dh.parms <- rbind(dh.parms,rand.par)
  }
  else next
}
sum(mod.summary$dh.weight)

colSums(mod.summary)
rowSums(mod.summary)
total.data <- rbind(ql.parms,log.parms,pro.parms,wei.parms,ms2.parms,llog.parms,lpro.parms,dh.parms)
total.data <- total.data[c("Dataset","Model","Par1","Par2","Par3","Par4")]
write.csv(total.data,"Total_data_092022.csv", row.names=FALSE)

################# End of BMD Modeling Functions ################

