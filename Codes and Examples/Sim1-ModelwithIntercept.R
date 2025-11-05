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

####################################################################################################
###################################### Simulation Study for Models with intercept ##################
####################################################################################################
## To ensure consistency with prior literature, this simulation study set up is the same as the simulation setting in Becker, Proksch, and Ringle (2022).
## The simulation study considers (1) the ECDF using the fixed percentile value for the maximum rank from Becker, Proksch, and Ringle (2022) in that algorithm in Equation 5 of the paper is
## (2) the ECDF using n/(n+1) for the maximum rank used to construct rank-based copula transformation.
## (3) the integrating kernel density estimate to cont

simulate_copula_intercept<-function(calculationContext){
## To ensure consistency with prior literature, this R function  is the same as the simulation setting and code
## provided in Becker, Proksch, and Ringle (2021).

  # Get variables from execution context
  # Simulation Repetition
  i = calculationContext$i
  # The sample size
  samplesize = calculationContext$samplesize
  # The error correlation rho
  rho = calculationContext$p1
  # The intercept level
  intercept = calculationContext$intercept
  # factorial combination index for seeding
  factor_index <- calculationContext$seed_index

  # Important to set a unique seed for parallel execution
  set.seed(factor_index*1000+i)

  # Create correlation matrix for regressor error correlation
  Sigma <- matrix(c(1,rho,rho,1),2,2)
  # Create random variables from correlation matrix
  dat <- data.frame(mvrnorm(samplesize,rep(0,2),Sigma,empirical = FALSE))

  colnames(dat) <- c("Epsilon","P_star")
  dat$P <- pnorm(dat$P_star)

  # Regression equation with specified intercept
  dat$Y <- intercept -1*dat$P + dat$Epsilon

  # Centering all variables
  centered_dat <- data.frame(scale(dat,center = TRUE, scale = FALSE))

  ################################################################################################
  ####  Results of regression equation without copula and without intercept (centered)
  ################################################################################################
  endog.c.fit       <- summary(lm(Y~0+P,data=centered_dat))

  endog.c.param_p   <- endog.c.fit$coefficients[1,1]
  endog.c.se_p      <- endog.c.fit$coefficients[1,2]

  ################################################################################################
  ####  Results of regression equation without copula and without intercept  (non-centered)
  ################################################################################################
  endog.fit        <- summary(lm(Y~0+P,data=dat))

  endog.param_p    <- endog.fit$coefficients[1,1]
  endog.se_p       <- endog.fit$coefficients[1,2]

  ################################################################################################
  ####   Results of regression equation without copula (normal OLS) with intercept (non-centered)
  ################################################################################################
  endog.i.fit       <- summary(lm(Y~1+P,data=dat))

  endog.i.param_p   <- endog.i.fit$coefficients[2,1]
  endog.i.param_i   <- endog.i.fit$coefficients[1,1]
  endog.i.se_p      <- endog.i.fit$coefficients[2,2]
  endog.i.se_i      <- endog.i.fit$coefficients[1,2]

  ################################################################################
  ####   Results of copula regression using the control function approach
  ####   (Park and Gupta Proposed Method I) generating P_Star via integrating nonparametric density estimation.
  ################################################################################
  dat$Copula1  <- createCopulaNew(dat$P, method="density")
  
  copula.d.fit.lm       <- lm(Y~0+P+Copula1,data=dat)          # without intercept (non-centered)
  copula.id.fit.lm     <- lm(Y~1+P+Copula1,data=dat)          # with intercept (non-centered)
  
  copula.d.fit        <- summary(copula.d.fit.lm)
  copula.id.fit        <- summary(copula.id.fit.lm)
 
  ## without intercept
  copula.d.param_p      <- copula.d.fit$coefficients[1]
  copula.d.param_rho  <- copula.d.fit$coefficients[2]
  
  ## with intercept
  copula.id.param_i      <- copula.id.fit$coefficients[1]
  copula.id.param_p      <- copula.id.fit$coefficients[2]
  copula.id.param_rho  <- copula.id.fit$coefficients[3]
  
  ################################################################################
  ####   Results of copula regression using the control function approach
  ####   (Park and Gupta Proposed Method II) generating P_Star via ecdf() using Becker et al. (2022)
  ################################################################################
  dat$Copula1  <- createCopula(dat$P)
  centered_dat$Copula2  <- createCopula(centered_dat$P)

  copula.fit.lm       <- lm(Y~0+P+Copula1,data=dat)          # without intercept (non-centered)
  copula.i.fit.lm     <- lm(Y~1+P+Copula1,data=dat)          # with intercept (non-centered)
  copula.c.fit.lm     <- lm(Y~0+P+Copula2,data=centered_dat) # without intercept (centered)

  copula.fit          <- summary(copula.fit.lm)
  copula.i.fit        <- summary(copula.i.fit.lm)
  copula.c.fit        <- summary(copula.c.fit.lm)

  # without intercept (non-centered)
  copula.param_p      <- copula.fit$coefficients[1,1]
  copula.param_pstar  <- copula.fit$coefficients[2,1]

  copula.param_p_vif      <- vif(copula.fit.lm)[1]
  copula.param_pstar_vif  <- vif(copula.fit.lm)[2]

  # with intercept (non-centered)
  copula.i.param_i      <- copula.i.fit$coefficients[1,1]
  copula.i.param_p      <- copula.i.fit$coefficients[2,1]
  copula.i.param_pstar  <- copula.i.fit$coefficients[3,1]

  copula.i.param_p_vif      <- vif(copula.i.fit.lm)[1]
  copula.i.param_pstar_vif  <- vif(copula.i.fit.lm)[2]

  # without intercept (centered)
  copula.c.param_p      <- copula.c.fit$coefficients[1,1]
  copula.c.param_pstar  <- copula.c.fit$coefficients[2,1]

  copula.c.param_p_vif      <- vif(copula.c.fit.lm)[1]
  copula.c.param_pstar_vif  <- vif(copula.c.fit.lm)[2]
  
  ################################################################################
  ####   Results of copula regression using the control function approach
  ####   (Park and Gupta Proposed Method II) generating P_Star via ecdf() using Adaptive algorithm. 
  ################################################################################
  dat$Copula1  <- createCopulaNew(dat$P, method="ecdf")
  centered_dat$Copula2  <- createCopulaNew(centered_dat$P, method="ecdf")
  
  copula.fit.lm       <- lm(Y~0+P+Copula1,data=dat)          # without intercept (non-centered)
  copula.i.fit.lm     <- lm(Y~1+P+Copula1,data=dat)          # with intercept (non-centered)
  copula.c.fit.lm     <- lm(Y~0+P+Copula2,data=centered_dat) # without intercept (centered)
  
  copula.fit          <- summary(copula.fit.lm)
  copula.i.fit        <- summary(copula.i.fit.lm)
  copula.c.fit        <- summary(copula.c.fit.lm)
  
  # without intercept (non-centered)
  copulaA.param_p      <- copula.fit$coefficients[1,1]
  copulaA.param_pstar  <- copula.fit$coefficients[2,1]
  
  copulaA.param_p_vif      <- vif(copula.fit.lm)[1]
  copulaA.param_pstar_vif  <- vif(copula.fit.lm)[2]
  
  # with intercept (non-centered)
  copulaA.i.param_i      <- copula.i.fit$coefficients[1,1]
  copulaA.i.param_p      <- copula.i.fit$coefficients[2,1]
  copulaA.i.param_pstar  <- copula.i.fit$coefficients[3,1]
  
  copulaA.i.param_p_vif      <- vif(copula.i.fit.lm)[1]
  copulaA.i.param_pstar_vif  <- vif(copula.i.fit.lm)[2]
  
  # without intercept (centered)
  copulaA.c.param_p      <- copula.c.fit$coefficients[1,1]
  copulaA.c.param_pstar  <- copula.c.fit$coefficients[2,1]
  
  copulaA.c.param_p_vif      <- vif(copula.c.fit.lm)[1]
  copulaA.c.param_pstar_vif  <- vif(copula.c.fit.lm)[2]
  

  result_list <- c(i,samplesize,rho,intercept,
                   # Endogeneous parameters
                   endog.param_p, endog.se_p,
                   endog.i.param_p, endog.i.se_p,
                   endog.i.param_i, endog.i.se_i,
                   endog.c.param_p, endog.c.se_p,
                   # Copula parameters using Becker et al. (2022)
                   copula.param_p,
                   copula.param_pstar, 
                   copula.i.param_p,
                   copula.i.param_pstar, 
                   copula.i.param_i, 
                   copula.c.param_p, 
                   copula.c.param_pstar, 
                   ## Copula parameters using Adaptive algorithm 
                   copulaA.param_p, 
                   copulaA.param_pstar, 
                   copulaA.i.param_p, 
                   copulaA.i.param_pstar, 
                   copulaA.i.param_i, 
                   copulaA.c.param_p, 
                   copulaA.c.param_pstar, #
                   ## Copula parameters with integrating nonparametric density estimates 
                   copula.d.param_p, 
                   copula.d.param_rho,
                   copula.id.param_p,
                   copula.id.param_rho,
                   copula.id.param_i 
                   )
  return(result_list)
}

#Package for parallel calculation
library(snowfall)


# Function to create Gausian Copula
# From Gui, Raluca, Markus Meierer, and Rene Algesheimer (2017), 
# "R Package REndo: Fitting Linear Models with Endogenous Regressors using
# Latent Instrumental Variables (Version 1.3)," https://cran.r-project.org/web/packages/REndo/
createCopula <- function(P){
  H.p <- stats::ecdf(P)
  H.p <- H.p(P)
  H.p <- ifelse(H.p==0,0.0000001,H.p)
  H.p <- ifelse(H.p==1,0.9999999,H.p)
  U.p <- H.p
  p.star <- stats::qnorm(U.p)
  return(p.star)	
}



createCopulaNew<- function(P, method="ecdf"){
  ## Method: ecdf
  if (method=="ecdf") {
    H.p <- stats::ecdf(P)
    H.p <- H.p(P)
    H.p <- ifelse(H.p==1,length(P)/(length(P)+1),H.p)
    U.p <- H.p
    p.star <- stats::qnorm(U.p)} else 
  if (method=="density") {
    density_estimate <-density(P, kernel="epanechnikov" )
    H.p <- CDF(density_estimate)(P)## cdf_from_density(P,density_estimate)
    p.star <- stats::qnorm(H.p)
  } else 
  {
    cat("Method not implmented")
  }
  return(p.star)
}

cdf_from_density <- function(x, density_estimate) {
  # Interpolation function for density
  density_function <- function(x) {
    approx(density_estimate$x, density_estimate$y, xout = x, rule = 2)$y
  }
  
  # Compute CDF using trapezoidal rule
  cdf_value <- sapply(x, function(x_val) {
    # Define range for integration
    x_range <- density_estimate$x[density_estimate$x <= x_val]
    if (length(x_range) == 0) return(0)
    
    # Density values corresponding to the x_range
    y_range <- density_function(x_range)
    
    # Apply trapezoidal rule
    trapz <- sum((diff(x_range) * (head(y_range, -1) + tail(y_range, -1)) / 2))
    return(trapz)
  })
  
  return(cdf_value)
}




# Create the design list for subsequent execution
calculationContextList <- list()
z <- 1
seed_ind <- 0
Simulationrepetitions <- 1000

# Parameters for simulation design
Samplesizes <- c(100, 200, 400, 600, 800, 1000, 2000, 4000, 6000, 8000, 10000, 20000, 40000, 60000)
Intercept <- 0

result <-NULL

for(s in Samplesizes)
{
  print("sample size= "); print(s)
  for (k in Intercept)
  {
    seed_ind <- seed_ind+1

    for(i in 0:(Simulationrepetitions-1))
    {
      if (i %in% c(10,200,500,800)) print(i)
      calculationContextList[[z]]<-list(i=i, samplesize=s, p1=0.5, intercept=k, seed_index=seed_ind)

      result <-rbind(result, simulate_copula_intercept(calculationContextList[[z]]))
      z<-z+1
     ## write.csv(result, "result24.csv")
    }
  }
}



colnames(result) <- c("i","samplesize","p1","intercept",
                      # Endogeneous parameters
                      "endog.param_p", "endog.se_p",
                      "endog.i.param_p", "endog.i.se_p",
                      "endog.i.param_i", "endog.i.se_i",
                      "endog.c.param_p", "endog.c.se_p",
                      # Copula parameters using Becker et al. (2022)
                      "copula.param_p",
                      "copula.param_pstar",
                      "copula.i.param_p", 
                      "copula.i.param_pstar",
                      "copula.i.param_i", 
                      "copula.c.param_p",
                      "copula.c.param_pstar",
                      ## Copula parameters uing Adaptive algorithm
                      "copulaA.param_p",
                      "copulaA.param_pstar",
                      "copulaA.i.param_p", 
                      "copulaA.i.param_pstar",
                      "copulaA.i.param_i", 
                      "copulaA.c.param_p",
                      "copulaA.c.param_pstar",
                      ## Copula parameters with integrating nonparametric density estimation (IKDE)
                      "copula.d.param_p", "copula.d.param_rho","copula.id.param_p", "copula.id.param_rho","copula.id.param_i"
)
bias<-aggregate(result, by=list(ss=result[,2]), FUN=mean)


################### Make Plots #########################################

postscript("biasIKDE.eps")
plot(log(bias$ss), bias$copulaA.i.param_p+1, type="b", xaxt="n", ylim=c(0,2), lty=1, pch=0, 
     xlab="Sample Size (Log Scale)", ylab="Parameter Bias",cex.lab=1.5)
axis(side=1, at=log(bias$ss), labels=bias$ss)
legend(log(100),0.25, legend="Copula regression with intercept", pch=0, bty="n")
text(log(400),0.12, "(this work)")
lines(log(bias$ss), bias$copula.i.param_p+1, type="b", pch=4, lty=2)
legend(log(600),0.5,  "Copula regression with intercept (Becker et al. 2022)",pch=4, bty="n")
lines(log(bias$ss), bias$endog.i.param_p+1, type="b", pch=1)
legend(log(600),1.9, "OLS with intercept",pch=1, bty="n")
lines(log(bias$ss),abs(bias$copula.id.param_p+1), type="b", pch=5)
legend(log(600),1.4, "Copula regression with intercept and IKDE",pch=5, bty="n")
dev.off()

