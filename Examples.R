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

##############################################################
### R functions to run 2sCOPE
################################################################
tscope<- function(formula, data, nboot=500, first.stage=F){
  ## Formula gives model: Y~ P+W |P, where Y is response, P is a vector of endogenous regressors and W is a vector of exogenous regressors
  ## data: names of data set
  ## nboot: number of bootstrap samples for estimating standard errors, defaults to 500
  ## first.stage: =F to obtain 2sCOPE result; =T to examine which W has first stage (Type III) F stat >10
  ## It implments the 2sCOPE method for Yang, Qian and Xie (2024) JMR paper.
  ## If no control variables, then 2sCOPE reduces to the P&G method.



  if (first.stage) {fun.tscope(formula, data, first.stage=T); return(cat(""))} else
  {
  # print output
  cat("========================================================",fill=TRUE)
  cat(paste("2sCOPE estimates", sep=" "),fill=TRUE)
  est<-fun.tscope(formula, data)[1,]
  # print output
  cat("========================================================",fill=TRUE)
  cat("Bootstrap standard errors",fill=TRUE)


  if (nboot==1)  res<-est else {
    se <- matrix(0, nboot, length(est))
    for (i in 1:nboot) {
      bootdata=data[sample(1:nrow(data), size=nrow(data), replace=T),]
      se[i,]<-fun.tscope(formula, bootdata)[1,]
    }
    res<-matrix(0, ncol=2, nrow=length(est))
    res[,1]<-est
    res[,2]<-apply(se,2,sd)
    rownames(res)<-names(est); colnames(res)<-c("Est", "boot.SE")
  }
  res
  }

}


fun.tscope <- function(formula, data, first.stage=F){
  F.formula <- as.Formula(formula)
  ymodel = formula(F.formula, lhs = 1, rhs = 1)
  x <- as.matrix(model.matrix(F.formula, data = data, rhs = 1)) #all regressors including intercept
  y <- unlist(model.part(F.formula, data = data, lhs = 1)[, 1]) #dependent variable
  n = nrow(x)
  if (ncol(x) <= 1)  print("No predictor variables specified for the outcome")

  #obtain endogenous regressors and exogenous regressors
  endox <- as.matrix(model.matrix(F.formula, data = data, rhs = length(F.formula)[2]))
  endox <- endox[,colnames(endox)!= "(Intercept)", drop=F] ## remove the intercept to obtain endogenous regressors only.
  nendox <- ncol(endox)
  endoxstar <- matrix(0,n,nendox)
  w <- x[ , -which(colnames(x) %in% c("(Intercept)",colnames(endox))), drop=F]
  wnames<-colnames(w)
  pnames<-colnames(endox)
  nw = ncol(w)
  wstar <- matrix(0,n,nw)

  #calculate xstar
  for (i in 1:nendox){
    endoxstartemp = ecdf(endox[,i])(endox[,i])
    endoxstartemp[endoxstartemp==1]=n/(n+1)
    endoxstartemp= qnorm(endoxstartemp)
    endoxstar[,i] = copy(endoxstartemp)
  }

  if ((nendox >= ncol(x)-1)) {
    ##print("No exogenous regressors specified for the outcome, and use copula_origin directly")
    ##print("using P&G")
    res <- lm(y ~ -1 + x + endoxstar)
  }else{
    #calculate wstar
    for (i in 1:nw){
      wstartemp = ecdf(w[,i])(w[,i])
      wstartemp[wstartemp==1]=n/(n+1)
      wstartemp= qnorm(wstartemp)
      wstar[,i] = copy(wstartemp)
    }
    stage1_resid = matrix(0,n,nendox)
    colnames(wstar)=wnames
    for (j in 1:nendox){
      dataj<-data.frame(pj=endoxstar[,j], wstar)
      regj<-lm(pj~., data=dataj)
      stage1_resid[,j]<- regj$resid
      if (first.stage) {
        stage1F<- Anova(regj, type=3)
        cat(paste("The following variable(s) have F-stat>10 in the first stage regression for:", pnames[j]))
        ##print(pnames[j], quote=F)
        cat("",fill=TRUE)
        cat("",fill=TRUE)
        print(as.data.frame(stage1F[stage1F[,3]>10 & !is.na(stage1F[,3]),]))
      }
    }
    res <- lm(y ~ -1 + x + stage1_resid)
  }

  resid2=y-c(x%*%res$coef[1:ncol(x)])
  corr_xerror_tscope = p_corr=rep(0,nendox)
  for (j in 1:nendox){
    corr_xerror_tscope[j]<- cor(endoxstar[,j],resid2)
    p_corr[j] <- cor.test(endoxstar[,j],resid2, method="pearson")$p.value
  }

  res.tab = matrix(0,2,length(res$coef)+nendox+1)
  res.tab[1,] = c(t(res$coef),t(corr_xerror_tscope),sd(resid2))
  res.tab[2, ] = c(summary(res)$coef[,2], p_corr, 0)
  colnames(res.tab)<-c(names(res$coef), rep("cor", nendox), "sigma")

  return(res.tab)
}




### Data sets for example 1 and 2 are perturbed data from the original data sets
### Thus they generate very close results but not exactly the same results as those presented in the paper
### Despite numerical differences, the main conclusions using these perturbed data remain the same as using the original data set.
############################################################################################
############################## Example 1: Demonstrate the use of copula correction  ####################################################
###########################################################################################

## Example 1
### Add into the code below the name of the folder where Example1.csv is downloaded and save to.
data1<-read.csv("Example1.csv")

store1.ols<-lm(logVol ~ logPrice + Fshare + week + Q2 + Q3 + Q4, data=data1)
summary(store1.ols)


## Check Step 2 and conclude insufficient nonnormality of endogenouos regressor
ks.test(scale(data1$logPrice), y="pnorm")

## Check Step 3 to find at least one exogenous regressor with sufficient nonnormality and association with P.

tscope(logVol ~ logPrice + Fshare  + week + Q2 + Q3 + Q4|logPrice, data=data1, first.stage=T)
ks.test(scale(data1$Fshare),y="pnorm")

## Run 2sCOPE
ymodel1 = logVol ~ logPrice + Fshare + week + Q2 + Q3 + Q4 |logPrice
store1.2scope = tscope(ymodel1, data=data1)
store1.2scope


######################################################################################
## Example 2: Demonstrate the empirical differences between copula correction
##            with and without copual interaction term.
#####################################################################################

### Add into the code below the name of the folder where Example2.csv is downloaded and save to.
data2 = read.csv("Example2.csv")

## OLS regression
store2.ols<-lm(logVol ~ logPrice + Fshare+ I(logPrice*Fshare) + week + Q2 + Q3 + Q4, data=data2)
summary(store2.ols)




## Check Step 2 and conclude sufficient nonnormality of endogenous regressor.
ks.test(scale(data2$logPrice), y="pnorm")

## Check Step 3 if endogenous regressor is correlated with exogenous regressor
## conclude correlation exist
pstar<-ecdf(data2$logPrice)(data2$logPrice); pstar[pstar==1]=length(pstar)/(length(pstar)+1); pstar<-qnorm(pstar)
cor.test(pstar, data2$Fshare);
cor.test(pstar, data2$week);
cor.test(pstar, data2$Q2);
cor.test(pstar, data2$Q3);
cor.test(pstar, data2$Q4);

## run 2sCOPE regression
## 2sCOPE without copula interaction term
example2.2scope=tscope(logVol ~ logPrice+ Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4|logPrice, data=data2)
example2.2scope

## 2sCOPE with copula interaction term
example2.2scopewint=tscope(logVol ~ logPrice+ Fshare + I(logPrice*Fshare) + week + Q2 + Q3 + Q4|logPrice +I(logPrice*Fshare),
                           data=data2)
example2.2scopewint

