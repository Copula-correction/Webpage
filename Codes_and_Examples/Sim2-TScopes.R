library(dplyr)
library(Formula)
library(data.table)
library(scales)
library(car)
library(MASS)
library(ivreg)
library(np)
library(mvtnorm)
library(truncnorm)


#####################################################################
## Use recent R version 4.3.0 or later
#####################################################################
if (getRversion() <"4.3.0") stop ("Use R version later than 4.3.0")

rm(list = ls())

## Load Rcope package
library(Rcope)

################################################################################################
######################################## Simulation Studies #####################################
##################################################################################################

nsim = 20
n = 500
ols.res = tscope.res =tscopenp.res = NULL


for (i in 1:nsim){
  
  print(i)
  temp <- rmvnorm(n, mean=c(0,0,0), sigma=cbind(c(1,0.5,0.5),
                                                c(0.5,1,0),
                                                c(0.5,0,1)))
  
  x= qgamma(pnorm(temp[,1]), shape=1, rate =1, lower.tail = TRUE,log.p = FALSE)
  w= qnorm(pnorm(temp[,2]))
  y <- 1 + x - w + temp[,3]
  
  ## Estimation:
  ##1. OLS estimation
  ols.reg<- lm(y~x+w)
  ols.res= rbind(ols.res, ols.reg$coef)
  
  ##2. 2sCOPE:
  data = as.data.frame(cbind(y,x,w))
  formula = 'y~x+w|x|w'
  res <- tscope(formula=formula, data=data,nboot=1)
  tscope.res = rbind(tscope.res, res[1:3,])
  
  ##3. 2sCOPE-np
  pccdf<-npcdist(x~w)
  s1star<- qnorm(pccdf$condist)
  s2reg<-lm(y~x+w+s1star)
  tscopenp.res <-rbind(tscopenp.res, coef(s2reg)[1:3])
  
  print(apply(ols.res[1:i,,drop=F],2,mean)) #OLS
  print(apply(tscope.res[1:i,,drop=F],2,mean)) #2sCOPE
  print(apply(tscopenp.res[1:i,, drop=F],2,mean)) #2scope-np
  
}



## Examine the estimates

apply(ols.res,2,mean)
apply(tscope.res,2,mean)
apply(tscopenp.res,2,mean)