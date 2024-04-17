

## This R code calls the Stan code for bivariate shrinkage priors 
## for the Efron's product means problem, i.e. where the parameter
## of interest is psi = theta_1*theta_2. 
## Please see the 2016 Biometrika paper: "Default Bayes" by Bhadra et al. 


rm(list=ls())

setwd("~/GitHub/BivariateHSstan")

source("bivar_priors_stan_0324.R")

set.seed(657)
hist.ci.pct = 0.95
A=10
J=2 
K = 100
theta.true = c(0,0)
## Gaussian Data
# test.data = list('J'=J, 'K'=K,
#                  'y'=t(replicate(K, rnorm(J,theta.true,1))),'psi' = 1, 'd' = 1)
## Non-Gaussian data 
test.data = list('J'=J, 'K'=K,
                 'y'=t(replicate(K, rt(J,df = 3))),'psi' = 1, 'd' = 1)
cat(colMeans(y),mean(y[,1]*y[,2]),mean(y[,1]/y[,2]),"\n",
    apply(y,2,median),median(y[,1]*y[,2]),median(y[,1]/y[,2]))

stan.hsplus.fit = stan_model(model_code=stan.hsplus.code, model_name="hs+ cauchy")
stan.hs.fit = stan_model(model_code=stan.hs.code, model_name="hs cauchy")
stan.norm.fit = stan_model(model_code=stan.norm.code, model_name="normal")
stan.lap.fit = stan_model(model_code=stan.lap.code, model_name="laplace")
stan.local.fit = stan_model(model_code=stan.local.code, model_name="local")
stan.global.fit = stan_model(model_code=stan.global.code, model_name="global")
stan.ref.fit = stan_model(model_code=stan.ref.code, model_name="reference")
stan.fc.fit = stan_model(model_code=stan.fc.code, model_name="creasy")
##########
n.iters = 1000
n.chains = 4
rng.seed = 768

## Horseshoe+ Results
smpls.hsplus.res = sampling(stan.hsplus.fit, 
                            data = test.data, 
                            iter = n.iters,
                            #init = 0,
                            seed = rng.seed, 
                            thin = 2, 
                            chains = n.chains)
theta.smpls.hsplus = extract(smpls.hsplus.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.hsplus) = c("theta1", "theta2")

hist.hsplus.ci.pct = 0.95

# hsplus.2d.density<- ggplot(as.data.frame(theta.smpls.hsplus), aes(theta1, theta2)) + 
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("HS+")
# print(hsplus.2d.density)

prior.prod.theta.hsplus = replicate(nrow(theta.smpls.hsplus), 
{
  tau = rcauchy(1, 0, 1)
  eta = rcauchy(2, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau * eta))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]*theta[2]
})
post.prod.theta.hsplus = apply(theta.smpls.hsplus, 1, function(x) x[1]*x[2])
prod.dist.hsplus = rbind(data.frame(type="prior", value=prior.prod.theta.hsplus),
                         data.frame(type="posterior", value=post.prod.theta.hsplus))

hist.quants.hsplus = quantile(prod.dist.hsplus$value, probs=c(hist.hsplus.ci.pct, 1.0-hist.hsplus.ci.pct))
hist.data.hsplus = subset(prod.dist.hsplus, min(hist.quants.hsplus) < value & value < max(hist.quants.hsplus))
hist.breaks.hsplus = hist(hist.data.hsplus$value, plot=FALSE, breaks="Scott")$breaks
# hsplus.hist<- ggplot(hist.data.hsplus, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.hsplus) + 
#   xlab(expression(theta[1]*theta[2]))
# print(hsplus.hist)
summary(prior.prod.theta.hsplus)
summary(post.prod.theta.hsplus)

### Horseshoe Results
smpls.hs.res = sampling(stan.hs.fit, 
                        data = test.data, 
                        iter = n.iters,
                        #init = 0,
                        seed = rng.seed, 
                        thin = 2, 
                        chains = n.chains)
theta.smpls.hs = extract(smpls.hs.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.hs) = c("theta1", "theta2")

hist.hs.ci.pct = 0.95

# hs.2d.density<- ggplot(as.data.frame(theta.smpls.hs), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("HS")
# print(hs.2d.density)

prior.prod.theta.hs = replicate(nrow(theta.smpls.hs), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]*theta[2]
})
post.prod.theta.hs = apply(theta.smpls.hs, 1, function(x) x[1]*x[2])
prod.dist.hs = rbind(data.frame(type="prior", value=prior.prod.theta.hs),
                     data.frame(type="posterior", value=post.prod.theta.hs))

hist.quants.hs = quantile(prod.dist.hs$value, probs=c(hist.hs.ci.pct, 1.0-hist.hs.ci.pct))
hist.data.hs = subset(prod.dist.hs, min(hist.quants.hs) < value & value < max(hist.quants.hs))
hist.breaks.hs = hist(hist.data.hs$value, plot=FALSE, breaks="Scott")$breaks
# hs.hist<- ggplot(hist.data.hs, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.hs) + 
#   xlab(expression(theta[1]*theta[2]))
# print(hs.hist)

summary(prior.prod.theta.hs)
summary(post.prod.theta.hs)

#----------------#
# Normal Results #
smpls.norm.res = sampling(stan.norm.fit, 
                          data = test.data, 
                          iter = n.iters,
                          #init = 0,
                          seed = rng.seed, 
                          thin = 2,
                          chains = n.chains)
theta.smpls.norm = extract(smpls.norm.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.norm) = c("theta1", "theta2")

hist.norm.ci.pct = 0.95

# norm.2d.density<-ggplot(as.data.frame(theta.smpls.norm), aes(theta1, theta2)) + 
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("Normal")
# print(norm.2d.density)

prior.prod.theta.norm = replicate(nrow(theta.smpls.norm), 
{
  theta = rnorm(2, 0, 1) 
  theta[1]*theta[2]
})
post.prod.theta.norm = apply(theta.smpls.norm, 1, function(x) x[1]*x[2])
prod.dist.norm = rbind(data.frame(type="prior", value=prior.prod.theta.norm),
                       data.frame(type="posterior", value=post.prod.theta.norm))

hist.quants.norm = quantile(prod.dist.norm$value, probs=c(hist.norm.ci.pct, 1.0-hist.norm.ci.pct))
hist.data.norm = subset(prod.dist.norm, min(hist.quants.norm) < value & value < max(hist.quants.norm))
hist.breaks.norm = hist(hist.data.norm$value, plot=FALSE, breaks="Scott")$breaks
# norm.hist<-ggplot(hist.data.norm, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.norm) + 
#   xlab(expression(theta[1]*theta[2]))
# print(norm.hist)
summary(prior.prod.theta.norm)
summary(post.prod.theta.norm)

#### Laplace 

smpls.lap.res = sampling(stan.lap.fit, 
                        data = test.data, 
                        iter = n.iters,
                        #init = 0,
                        seed = rng.seed, 
                        thin = 2,
                        chains = n.chains)
theta.smpls.lap = extract(smpls.lap.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.lap) = c("theta1", "theta2")

hist.lap.ci.pct = 0.95

# lap.2d.density<- ggplot(as.data.frame(theta.smpls.lap), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("lap")
# print(lap.2d.density)

prior.prod.theta.lap = replicate(nrow(theta.smpls.lap), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]*theta[2]
})
post.prod.theta.lap = apply(theta.smpls.lap, 1, function(x) x[1]*x[2])
prod.dist.lap = rbind(data.frame(type="prior", value=prior.prod.theta.lap),
                     data.frame(type="posterior", value=post.prod.theta.lap))

hist.quants.lap = quantile(prod.dist.lap$value, probs=c(hist.lap.ci.pct, 1.0-hist.lap.ci.pct))
hist.data.lap = subset(prod.dist.lap, min(hist.quants.lap) < value & value < max(hist.quants.lap))
hist.breaks.lap = hist(hist.data.lap$value, plot=FALSE, breaks="Scott")$breaks
# lap.hist<- ggplot(hist.data.lap, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.lap) + 
#   xlab(expression(theta[1]*theta[2]))
# print(lap.hist)

summary(prior.prod.theta.lap)
summary(post.prod.theta.lap)

#########################
##------- LOCAL -------##

smpls.local.res = sampling(stan.local.fit, 
                        data = test.data, 
                        iter = n.iters,
                        #init = 0,
                        seed = rng.seed, 
                        thin = 2,
                        chains = n.chains)
theta.smpls.local = extract(smpls.local.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.local) = c("theta1", "theta2")

hist.local.ci.pct = 0.95

# local.2d.density<- ggplot(as.data.frame(theta.smpls.local), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("local")
# print(local.2d.density)

prior.prod.theta.local = replicate(nrow(theta.smpls.local), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]*theta[2]
})
post.prod.theta.local = apply(theta.smpls.local, 1, function(x) x[1]*x[2])
prod.dist.local = rbind(data.frame(type="prior", value=prior.prod.theta.local),
                     data.frame(type="posterior", value=post.prod.theta.local))

hist.quants.local = quantile(prod.dist.local$value, probs=c(hist.local.ci.pct, 1.0-hist.local.ci.pct))
hist.data.local = subset(prod.dist.local, min(hist.quants.local) < value & value < max(hist.quants.local))
hist.breaks.local = hist(hist.data.local$value, plot=FALSE, breaks="Scott")$breaks
# local.hist<- ggplot(hist.data.local, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.local) + 
#   xlab(expression(theta[1]*theta[2]))
# print(local.hist)

summary(prior.prod.theta.local)
summary(post.prod.theta.local)

############################
##------- Global -------###

smpls.global.res = sampling(stan.global.fit, 
                        data = test.data, 
                        iter = n.iters,
                        #init = 0,
                        seed = rng.seed, 
                        thin = 2,
                        chains = n.chains)
theta.smpls.global = extract(smpls.global.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.global) = c("theta1", "theta2")

hist.global.ci.pct = 0.95

# global.2d.density<- ggplot(as.data.frame(theta.smpls.global), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("global")
# print(global.2d.density)

prior.prod.theta.global = replicate(nrow(theta.smpls.global), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]*theta[2]
})
post.prod.theta.global = apply(theta.smpls.global, 1, function(x) x[1]*x[2])
prod.dist.global = rbind(data.frame(type="prior", value=prior.prod.theta.global),
                     data.frame(type="posterior", value=post.prod.theta.global))

hist.quants.global = quantile(prod.dist.global$value, probs=c(hist.global.ci.pct, 1.0-hist.global.ci.pct))
hist.data.global = subset(prod.dist.global, min(hist.quants.global) < value & value < max(hist.quants.global))
hist.breaks.global = hist(hist.data.global$value, plot=FALSE, breaks="Scott")$breaks
# global.hist<- ggplot(hist.data.global, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.global) + 
#   xlab(expression(theta[1]*theta[2]))
# print(global.hist)

summary(prior.prod.theta.global)
summary(post.prod.theta.global)

############ Reference ########
smpls.ref.res = sampling(stan.ref.fit, 
                            data = test.data, 
                            iter = n.iters,
                            #init = 0,
                            seed = rng.seed, 
                            thin = 2,
                            chains = n.chains)
theta.smpls.ref = extract(smpls.ref.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.ref) = c("theta1", "theta2")

hist.ref.ci.pct = 0.95

# ref.2d.density<- ggplot(as.data.frame(theta.smpls.ref), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("ref")
# print(ref.2d.density)

prior.prod.theta.ref = replicate(nrow(theta.smpls.ref), 
                                    {
                                      tau = rcauchy(1, 0, 1)
                                      lambda = rcauchy(2, c(0,0), abs(tau))
                                      theta = rnorm(2, lambda, c(1,1)) 
                                      theta[1]*theta[2]
                                    })
post.prod.theta.ref = apply(theta.smpls.ref, 1, function(x) x[1]*x[2])
prod.dist.ref = rbind(data.frame(type="prior", value=prior.prod.theta.ref),
                         data.frame(type="posterior", value=post.prod.theta.ref))

hist.quants.ref = quantile(prod.dist.ref$value, probs=c(hist.ref.ci.pct, 1.0-hist.ref.ci.pct))
hist.data.ref = subset(prod.dist.ref, min(hist.quants.ref) < value & value < max(hist.quants.ref))
hist.breaks.ref = hist(hist.data.ref$value, plot=FALSE, breaks="Scott")$breaks
# ref.hist<- ggplot(hist.data.ref, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.ref) + 
#   xlab(expression(theta[1]*theta[2]))
# print(ref.hist)

summary(prior.prod.theta.ref)
summary(post.prod.theta.ref)

###############
smpls.fc.res = sampling(stan.fc.fit, 
                            data = test.data, 
                            iter = n.iters,
                            #init = 0,
                            seed = rng.seed, 
                            thin = 2, 
                            chains = n.chains)

theta.smpls.fc = extract(smpls.fc.res, pars=c("theta"), permuted=TRUE)[[1]]
colnames(theta.smpls.fc) = c("theta1", "theta2")

hist.fc.ci.pct = 0.95

# fc.2d.density<- ggplot(as.data.frame(theta.smpls.fc), aes(theta1, theta2)) +
#   geom_density2d()+stat_density2d(aes(fill = ..level..), geom="polygon")+
#   xlab(expression(theta[1]))+ylab(expression(theta[2]))+ggtitle("fc")
# print(fc.2d.density)

prior.prod.theta.fc = replicate(nrow(theta.smpls.fc), 
                                 {
                                   tau = rcauchy(1, 0, 1)
                                   lambda = rcauchy(2, c(0,0), abs(tau))
                                   theta = rnorm(2, lambda, c(1,1)) 
                                   theta[1]*theta[2]
                                 })
post.prod.theta.fc = apply(theta.smpls.fc, 1, function(x) x[1]*x[2])
prod.dist.fc = rbind(data.frame(type="prior", value=prior.prod.theta.fc),
                      data.frame(type="posterior", value=post.prod.theta.fc))

hist.quants.fc = quantile(prod.dist.fc$value, probs=c(hist.fc.ci.pct, 1.0-hist.fc.ci.pct))
hist.data.fc = subset(prod.dist.fc, min(hist.quants.fc) < value & value < max(hist.quants.fc))
hist.breaks.fc = hist(hist.data.fc$value, plot=FALSE, breaks="Scott")$breaks
# fc.hist<- ggplot(hist.data.fc, 
#        aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.fc) + 
#   xlab(expression(theta[1]*theta[2]))
# print(fc.hist)

summary(prior.prod.theta.fc)
summary(post.prod.theta.fc)

#-------------------------#
# library(gridExtra)
# setwd("C:/Users/Jyotishka/OneDrive/Documents/latex files/hsplus-applied-newer/art")
# # cairo_ps(file='all_2d_fc.eps', width=12, height=7)
# # grid.arrange(hsplus.2d.density,hs.2d.density,norm.2d.density,ncol=3)
# # dev.off()
# cairo_ps(file='all_hist_pm.eps', width=12, height=7)
# grid.arrange(hsplus.hist,hs.hist,norm.hist,ncol=3)
# dev.off()

##########

pm.data = rbind(data.frame(method="HSPlus",theta1=theta.smpls.hsplus[,1],theta2=theta.smpls.hsplus[,2]),
                data.frame(method="HS",theta1=theta.smpls.hs[,1],theta2=theta.smpls.hs[,2]),
                data.frame(method="Laplace",theta1=theta.smpls.lap[,1],theta2=theta.smpls.lap[,2]),
                data.frame(method="Normal",theta1=theta.smpls.norm[,1],theta2=theta.smpls.norm[,2]),
                data.frame(method="Pure-Local",theta1=theta.smpls.local[,1],theta2=theta.smpls.local[,2]),
                data.frame(method="Pure-Global",theta1=theta.smpls.global[,1],theta2=theta.smpls.global[,2]),
                data.frame(method="Product Mean",theta1=theta.smpls.ref[,1],theta2=theta.smpls.ref[,2]),
                data.frame(method="Fieller-Creasy",theta1=theta.smpls.fc[,1],theta2=theta.smpls.fc[,2]))

## 04/26: removed fill =..level.. from stat_density2d args and added colour = "black"
setwd("C:\\Users\\Jyotishka\\Dropbox\\biomet-hs\\biometrika_revised_draft\\art")
pm.data<-read.table(paste("0426_FCPM.csv",sep=""),sep=",") ##Normal plots
## For t-distributed
pm.data<-read.table(paste("0426_t_FCPM.csv",sep=""),sep=",") ##t plots

library(plyr)
pm.data$method = revalue(pm.data$method, c("HSPlus"="HS+"))

var.plot<-ggplot(pm.data, aes(theta1, theta2))+#geom_point()+
  geom_density2d()+stat_density2d(aes(alpha = ..level..),geom="polygon",colour="black")+
    scale_alpha(range = c(0.0, 0.5)) + 
    theme(legend.position = "none", axis.title = element_blank(), text = element_text(size = 12))+
  xlim(-0.5,0.5)+theme_bw()+#ylim(-0.1,0.4)+
  geom_vline(xintercept=0)+geom_hline(yintercept=0)+scale_x_continuous(breaks=c(-0.25,0,0.25))+
  xlab(expression(theta[1]))+ylab(expression(theta[2]))+facet_wrap(~method,ncol=4)

var.plot <- var.plot + theme(axis.title.y = element_text(size = rel(2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(2)))
var.plot<- var.plot+ theme(axis.text = element_text(size = rel(2)))+
  theme(legend.title=element_text(size=15, face="bold"),legend.text=element_text(size=18))
var.plot <- var.plot+theme(strip.text.x = element_text(size=18, face="bold"),strip.text.y = element_text(size=18, face="bold"))
print(var.plot)

setwd("C:\\Users\\Jyotishka\\Dropbox\\biomet-hs\\biometrika_revised_draft\\art\\revised_060716")
# setwd("C:/Users/Jyotishka/OneDrive/Documents/latex files/hsplus-applied-newer/art")
# ggsave(var.plot,file='7-contour-plots-pm.eps', width=7, height=6)
# 
# setwd("C:/Users/Jyotishka/OneDrive/Documents/latex files/hsplus-applied-newer/art")
cairo_ps(file='t-8-contour-plots-grayscale-legend.eps', width=10, height=6)
var.plot
dev.off()

#################
## Histogram for Product Means Problem ##
## I am not sure if this adds any insight ##
#################

hist.data.hs$method="HS"
#hist.data.hs$breaks = hist.breaks.hs
hist.data.hsplus$method="HSPlus"
#hist.data.hsplus$breaks = hist.breaks.hsplus
hist.data.norm$method="Normal"
#hist.data.norm$breaks = hist.breaks.norm

hist.data = rbind(hist.data.hsplus,hist.data.hs,hist.data.norm)

hist.plot<-ggplot(hist.data,aes(x=value, group=type))+
  geom_histogram(aes(fill=type, y=..density..)) + 
  xlab(expression(theta[1]*theta[2]))+facet_wrap(~method,scales="free")
print(hist.plot)

###############################################################
###############################################################
################# Fieller _ Creasy Problem ###################
###############################################################
###############################################################

## Only calculate Posterior Dist of ratios ##
## Horseshoe+ Prior ##
prior.ratio.theta.hsplus = replicate(nrow(theta.smpls.hsplus), 
{
  tau = rcauchy(1, 0, 1)
  eta = rcauchy(2, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau * eta))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]/theta[2]
})
post.ratio.theta.hsplus = apply(theta.smpls.hsplus, 1, function(x) x[1]/x[2])
ratio.dist.hsplus = rbind(data.frame(type="prior", value=prior.ratio.theta.hsplus),
                          data.frame(type="posterior", value=post.ratio.theta.hsplus))

hist.quants.hsplus = quantile(ratio.dist.hsplus$value, probs=c(hist.hsplus.ci.pct, 1.0-hist.hsplus.ci.pct))
hist.data.hsplus = subset(ratio.dist.hsplus, min(hist.quants.hsplus) < value & value < max(hist.quants.hsplus))
hist.breaks.hsplus = hist(hist.data.hsplus$value, plot=FALSE, breaks="Scott")$breaks
hsplus.hist<- ggplot(hist.data.hsplus, 
                     aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.hsplus) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe+")
print(hsplus.hist)

summary(prior.ratio.theta.hsplus)
summary(post.ratio.theta.hsplus)

##########################
###----- Horseshoe Prior ----###

prior.ratio.theta.hs = replicate(nrow(theta.smpls.hs), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]/theta[2]
})
post.ratio.theta.hs = apply(theta.smpls.hs, 1, function(x) x[1]/x[2])
ratio.dist.hs = rbind(data.frame(type="prior", value=prior.ratio.theta.hs),
                      data.frame(type="posterior", value=post.ratio.theta.hs))

hist.quants.hs = quantile(ratio.dist.hs$value, probs=c(hist.hs.ci.pct, 1.0-hist.hs.ci.pct))
hist.data.hs = subset(ratio.dist.hs, min(hist.quants.hs) < value & value < max(hist.quants.hs))
hist.breaks.hs = hist(hist.data.hs$value, plot=FALSE, breaks="Scott")$breaks
hs.hist<- ggplot(hist.data.hs, 
                 aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.hs) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(hs.hist)

summary(prior.ratio.theta.hs)
summary(post.ratio.theta.hs)

##########################################
####---------- Normal ---------###########

prior.ratio.theta.norm = replicate(nrow(theta.smpls.norm), 
{
  theta = rnorm(2, 0, 1) 
  theta[1]/theta[2]
})
post.ratio.theta.norm = apply(theta.smpls.norm, 1, function(x) x[1]/x[2])
ratio.dist.norm = rbind(data.frame(type="prior", value=prior.ratio.theta.norm),
                        data.frame(type="posterior", value=post.ratio.theta.norm))

hist.quants.norm = quantile(ratio.dist.norm$value, probs=c(hist.norm.ci.pct, 1.0-hist.norm.ci.pct))
hist.data.norm = subset(ratio.dist.norm, min(hist.quants.norm) < value & value < max(hist.quants.norm))
hist.breaks.norm = hist(hist.data.norm$value, plot=FALSE, breaks="Scott")$breaks
norm.hist<-ggplot(hist.data.norm, 
                  aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.norm) + 
  xlab(expression(theta[1]/theta[2])) + ggtitle("Normal")
print(norm.hist)

## Summary
summary(prior.ratio.theta.norm)
summary(post.ratio.theta.norm)

#################################################
########-------- Laplace--------#################
prior.ratio.theta.lap = replicate(nrow(theta.smpls.lap), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]/theta[2]
})
post.ratio.theta.lap = apply(theta.smpls.lap, 1, function(x) x[1]/x[2])
ratio.dist.lap = rbind(data.frame(type="prior", value=prior.ratio.theta.lap),
                      data.frame(type="posterior", value=post.ratio.theta.lap))

hist.quants.lap = quantile(ratio.dist.lap$value, probs=c(hist.lap.ci.pct, 1.0-hist.lap.ci.pct))
hist.data.lap = subset(ratio.dist.lap, min(hist.quants.lap) < value & value < max(hist.quants.lap))
hist.breaks.lap = hist(hist.data.lap$value, plot=FALSE, breaks="Scott")$breaks
lap.hist<- ggplot(hist.data.lap, 
                 aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.lap) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(lap.hist)

summary(prior.ratio.theta.lap)
summary(post.ratio.theta.lap)

##########################################
############-------- LOCAL-------#########

prior.ratio.theta.local = replicate(nrow(theta.smpls.local), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]/theta[2]
})
post.ratio.theta.local = apply(theta.smpls.local, 1, function(x) x[1]/x[2])
ratio.dist.local = rbind(data.frame(type="prior", value=prior.ratio.theta.local),
                       data.frame(type="posterior", value=post.ratio.theta.local))

hist.quants.local = quantile(ratio.dist.local$value, probs=c(hist.local.ci.pct, 1.0-hist.local.ci.pct))
hist.data.local = subset(ratio.dist.local, min(hist.quants.local) < value & value < max(hist.quants.local))
hist.breaks.local = hist(hist.data.local$value, plot=FALSE, breaks="Scott")$breaks
local.hist<- ggplot(hist.data.local, 
                  aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.local) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(local.hist)

summary(prior.ratio.theta.local)
summary(post.ratio.theta.local)

############################################
######---- Global-------#########

prior.ratio.theta.global = replicate(nrow(theta.smpls.global), 
{
  tau = rcauchy(1, 0, 1)
  lambda = rcauchy(2, c(0,0), abs(tau))
  theta = rnorm(2, lambda, c(1,1)) 
  theta[1]/theta[2]
})
post.ratio.theta.global = apply(theta.smpls.global, 1, function(x) x[1]/x[2])
ratio.dist.global = rbind(data.frame(type="prior", value=prior.ratio.theta.global),
                         data.frame(type="posterior", value=post.ratio.theta.global))

hist.quants.global = quantile(ratio.dist.global$value, probs=c(hist.global.ci.pct, 1.0-hist.global.ci.pct))
hist.data.global = subset(ratio.dist.global, min(hist.quants.global) < value & value < max(hist.quants.global))
hist.breaks.global = hist(hist.data.global$value, plot=FALSE, breaks="Scott")$breaks
global.hist<- ggplot(hist.data.global, 
                    aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.global) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(global.hist)

summary(prior.ratio.theta.global)
summary(post.ratio.theta.global)

######---- Product-------#########

prior.ratio.theta.ref = replicate(nrow(theta.smpls.ref), 
                                     {
                                       tau = rcauchy(1, 0, 1)
                                       lambda = rcauchy(2, c(0,0), abs(tau))
                                       theta = rnorm(2, lambda, c(1,1)) 
                                       theta[1]/theta[2]
                                     })
post.ratio.theta.ref = apply(theta.smpls.ref, 1, function(x) x[1]/x[2])
ratio.dist.ref = rbind(data.frame(type="prior", value=prior.ratio.theta.ref),
                          data.frame(type="posterior", value=post.ratio.theta.ref))

hist.quants.ref = quantile(ratio.dist.ref$value, probs=c(hist.ref.ci.pct, 1.0-hist.ref.ci.pct))
hist.data.ref = subset(ratio.dist.ref, min(hist.quants.ref) < value & value < max(hist.quants.ref))
hist.breaks.ref = hist(hist.data.ref$value, plot=FALSE, breaks="Scott")$breaks
ref.hist<- ggplot(hist.data.ref, 
                     aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.ref) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(ref.hist)

summary(prior.ratio.theta.ref)
summary(post.ratio.theta.ref)

############## Fieller Creasy #########

prior.ratio.theta.fc = replicate(nrow(theta.smpls.fc), 
                                  {
                                    tau = rcauchy(1, 0, 1)
                                    lambda = rcauchy(2, c(0,0), abs(tau))
                                    theta = rnorm(2, lambda, c(1,1)) 
                                    theta[1]/theta[2]
                                  })
post.ratio.theta.fc = apply(theta.smpls.fc, 1, function(x) x[1]/x[2])
ratio.dist.fc = rbind(data.frame(type="prior", value=prior.ratio.theta.fc),
                       data.frame(type="posterior", value=post.ratio.theta.fc))

hist.quants.fc = quantile(ratio.dist.fc$value, probs=c(hist.fc.ci.pct, 1.0-hist.fc.ci.pct))
hist.data.fc = subset(ratio.dist.fc, min(hist.quants.fc) < value & value < max(hist.quants.fc))
hist.breaks.fc = hist(hist.data.fc$value, plot=FALSE, breaks="Scott")$breaks
fc.hist<- ggplot(hist.data.fc, 
                  aes(x=value, group=type)) + geom_histogram(aes(fill=type, y=..density..), breaks=hist.breaks.fc) + 
  xlab(expression(theta[1]/theta[2]))+ggtitle("Horseshoe")
print(fc.hist)

summary(prior.ratio.theta.fc)
summary(post.ratio.theta.fc)

############## Table ###########
rbind(c(summary(post.ratio.theta.hsplus),sd(post.ratio.theta.hsplus)),
      c(summary(post.ratio.theta.hs),sd(post.ratio.theta.hs)),
      c(summary(post.ratio.theta.lap),sd(post.ratio.theta.lap)),
      c(summary(post.ratio.theta.norm),sd(post.ratio.theta.norm)),
      c(summary(post.ratio.theta.local),sd(post.ratio.theta.local)),
      c(summary(post.ratio.theta.global),sd(post.ratio.theta.global)),
      c(summary(post.ratio.theta.fc),sd(post.ratio.theta.fc)))

rbind(c(summary(post.prod.theta.hsplus),sd(post.prod.theta.hsplus)),
      c(summary(post.prod.theta.hs),sd(post.prod.theta.hs)),
      c(summary(post.prod.theta.lap),sd(post.prod.theta.lap)),
      c(summary(post.prod.theta.norm),sd(post.prod.theta.norm)),
      c(summary(post.prod.theta.local),sd(post.prod.theta.local)),
      c(summary(post.prod.theta.global),sd(post.prod.theta.global)),
      c(summary(post.prod.theta.ref),sd(post.prod.theta.ref)))

setwd("C:\\Users\\Jyotishka\\Dropbox\\biomet-hs\\biometrika_revised_draft\\art")
write.table(pm.data, paste("0426_t_FCPM.csv",sep=""),sep=",")
write.table(test.data, paste("0426_t_testdata_FCPM.csv",sep=""),sep=",")

pm.data<-read.table(paste("0426_t_FCPM.csv",sep=""),sep=",") ##t plots
pm.data<-read.table(paste("0426_FCPM.csv",sep=""),sep=",") ##Normal plots