library(utils)
library(stats4)
library(rstan)
library(MASS)
library(plyr)
library(dplyr)
library(reshape2)

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
