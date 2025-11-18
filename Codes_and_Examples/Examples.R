library(dplyr)
library(Formula)
library(boot)
library(data.table)
library(mvtnorm)
library(truncnorm)
library(scales)
library(car)
library(spatstat)
library(MASS)

rm(list = ls())

library(Rcope) #Load Rcope package

### Data sets for example 1 and 2 are perturbed data from the original data sets
### Thus they generate very close results but not exactly the same results as those presented in the paper
### Despite numerical differences, the main conclusions using these perturbed data remain the same as using the original data set.
############################################################################################
############################## Example 1: Demonstrate the use of copula correction  ####################################################
###########################################################################################

## Example 1
### Add into the code below the name of the folder where Example1.csv is downloaded and save to.
data1<-read.csv("C:/Users/anton/Desktop/school/package/tscope scratch files/webpage/Codes_and_Examples/Example1.csv")

store1.ols<-lm(logVol ~ logPrice + Fshare + week + Q2 + Q3 + Q4, data=data1)
summary(store1.ols)


## Check Step 2 and conclude insufficient nonnormality of endogenouos regressor
ks.test(scale(data1$logPrice), y="pnorm")

## Check Step 3 to find at least one exogenous regressor with sufficient nonnormality and association with P.

ccf_data1 <- CCF(logVol ~ logPrice + Fshare  + week + Q2 + Q3 + Q4|logPrice|Fshare  + week + Q2 + Q3 + Q4, data=data1)
relevance_test(ccf_data1)
ks.test(scale(data1$Fshare),y="pnorm")

## Run 2sCOPE
ymodel1 = logVol ~ logPrice + Fshare + week + Q2 + Q3 + Q4 |logPrice|Fshare + week + Q2 + Q3 + Q4
store1.2scope = tscope(ymodel1, data=data1, nboot = 1)
store1.2scope


######################################################################################
## Example 2: Demonstrate the empirical differences between copula correction
##            with and without copula interaction term.
#####################################################################################

### Add into the code below the name of the folder where Example2.csv is downloaded and save to.
data2 = read.csv("C:/Users/anton/Desktop/school/package/tscope scratch files/webpage/Codes_and_Examples/Example2.csv")

## OLS regression
store2.ols<-lm(logVol ~ logPrice + Fshare+ I(logPrice*Fshare) + week + Q2 + Q3 + Q4, data=data2)
summary(store2.ols)


## Check Step 2 and conclude sufficient nonnormality of endogenous regressor.
ks.test(scale(data2$logPrice), y="pnorm")

## Check Step 3 if endogenous regressor is correlated with exogenous regressor
## conclude correlation exist
ccf_data2 <- CCF(logVol ~ logPrice+ Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4|logPrice|
                   Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4, data=data2)

pstar <- ccf_data2$pstar
cor.test(pstar, data2$Fshare);
cor.test(pstar, data2$week);
cor.test(pstar, data2$Q2);
cor.test(pstar, data2$Q3);
cor.test(pstar, data2$Q4);

## run 2sCOPE regression
## 2sCOPE without copula interaction term
example2.2scope=tscope(logVol ~ logPrice+ Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4|logPrice|
                         Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4, data=data2, nboot = 1)
example2.2scope

## 2sCOPE with copula interaction term
example2.2scopewint=tscope(logVol ~ logPrice+ Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4|logPrice +I(logPrice*Fshare)|
                             Fshare + week + Q2 + Q3 + Q4,
                           data=data2, nboot = 1)
example2.2scopewint

