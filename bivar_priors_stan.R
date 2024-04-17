## Last modified: 04/26:
## Normal prior now set to variance = 300 and fixed a typo in FC prior code. 

library(rstan)
#set_cppo("fast")
library(ggplot2)
library(plyr)
library(reshape2)

########### Horseshoe Plus ########
stan.hsplus.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta_step; 
vector<lower=0>[J] lambda; 
vector<lower=0>[J] eta;
real<lower=0> tau;
}
transformed parameters {
vector[J] theta; 
theta <- ((theta_step .* lambda) .* eta) * tau;
}
model {
tau ~ cauchy(0, 1);
eta ~ cauchy(0, 1);
lambda ~ cauchy(0, 1);
theta_step ~ normal(0, 1);
for (k in 1:K) {
y[k] ~ normal(theta, 1);
}
}  
"
######### Horseshoe Prior ########
stan.hs.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta_step; 
vector<lower=0>[J] lambda; 
real<lower=0> tau;
}
transformed parameters {
vector[J] theta; 
theta <- (theta_step .* lambda) * tau;
}
model {
tau ~ cauchy(0, 1);
lambda ~ cauchy(0, 1);
theta_step ~ normal(0, 1);
for (k in 1:K) {
y[k] ~ normal(theta, 1);
}
}  
"
########## Normal vague, N(0, sqrt(300)) prior ####
stan.norm.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta; 
}
model {
theta ~ normal(0, sqrt(300));
for (k in 1:K) {
y[k] ~ normal(theta, 1);
}
}  
"
##### Laplace Shrinkage Prior #####
stan.lap.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
real<lower=0> psi;
real<lower=0> d;
}
parameters {
vector[J] theta_step; 
vector<lower=0>[J] lambda_step;
real<lower=0> tau;
}
transformed parameters {
vector[J] theta; 
vector<lower =0>[J] lambda;
for (j in 1:J)
lambda[j] <- sqrt(lambda_step[j]);
theta <- (theta_step .* lambda);
}
model {
tau ~ inv_gamma(psi/2,psi*d^2/2);
lambda_step ~ exponential(1/(2*tau^2));
theta_step ~ normal(0, 1);
for (k in 1:K) {
y[k] ~ normal(theta, 1);
}
}  
"
##### Global Shinkage Code ###
stan.global.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta; 
real<lower=0> tau;
}
model {
tau ~ cauchy(0, 1);
theta ~ normal(0, tau);
for (k in 1:K) {
  y[k] ~ normal(theta, 1);
}
}  
"
##### Pure Local Shrinkage ####
stan.local.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta_step; 
vector<lower=0>[J] lambda; 
}
transformed parameters {
vector[J] theta; 
theta <- (theta_step .* lambda);
}
model {
lambda ~ cauchy(0, 1);
theta_step ~ normal(0, 1);
for (k in 1:K) {
  y[k] ~ normal(theta, 1);
}
}  
"
######## Reference Prior for Product Means ######
stan.ref.code = "
data {
int<lower=0> K; 
int<lower=0> J; 
vector[J] y[K]; 
}
parameters {
vector[J] theta; 
}
model {
increment_log_prob(0.5*log(theta'*theta));
for (k in 1:K) {
y[k] ~ normal(theta, 1);
}
}  
"
######## Reference prior for Fieller-Creasy ###

stan.fc.code = "
data {
int<lower=0> K; 
vector[2] y[K]; 
}
parameters {
vector[2] theta; 
}
transformed parameters {
real nu;
real lambda;
nu <- theta[1]/theta[2];
lambda <- theta[2];
}
model {
##increment_log_prob(0.5*log(theta'*theta));
for (k in 1:K) {
increment_log_prob(-0.5*log(1+nu^2));
y[k] ~ normal(theta, 1);
}
}
"
#stan.fc.fit = stan_model(model_code=stan.fc.code, model_name="creasy")