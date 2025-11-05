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


##############################################################
### R functions to run 2sCOPE
################################################################
tscope<- function(formula, data, nboot=500, first.stage=F){
  ## Formula gives model: Y~ P+W |P, where Y is response, P is a vector of endogenous regressors and W is a vector of exogenous regressors
  ## data: names of data set
  ## nboot: number of bootstrap samples for estimating standard errors, defaults to 500
  ## first.stage: =F to obtain 2sCOPE result; =T to examine which W has first stage (Type III) F stat >10
  ## It implements the 2sCOPE method for Yang, Qian and Xie (2024) JMR paper.
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




################################################################################################
######################################## Simulation Studies #####################################
##################################################################################################

nsim = 1000
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
  formula = 'y~x+w|x'
  res <- tscope(formula=formula, data=data,nboot=1)
  tscope.res = rbind(tscope.res, res[1:3])

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



