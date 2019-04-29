# Clear workspace
rm(list=ls())

# Required packages
library(tidyverse)
library(socialmixr) # POLYMOD data
library(eurostat)   # EUROSTAT data

####################
## Initialisation ##
####################

## Age-dependent
age <- sero$AGE
y <- sero$indic
rij <- contact_w
muy = predict(demfit, type = "response")
N <- sum(PS)
D <- 6 / 365
A <- 0.5
Lmax <- 80 
prop <- "loglin"
print <- 0
plots <- plot
prop_fac <- prop


nEstimate <- 1 # no. of parameters to be estimated 
nIter <- 1e4 # number of iterations per chain
propScalar <- nEstimate * 1
pilotVar <- 1 


############################
## Set uniform priors for ##
## input parameters       ##
############################

prior <- c(1e-6, 1e3)
## check this
disp.prop <-  1e-6 * prior[2]
acc <- 0

mcmc.out <- rep(NA, (nEstimate+1)*nIter)
dim(mcmc.out) <- c(nIter,(nEstimate+1))

######################
## Initialise chain ##
######################

param0 <- c(0.5, 0.3) # Two values for age-dependent FoI
lnLike0 <- lnLike(param0, otherParams, countryData)

#########################
## Likelihood function ##
#########################
logLikeHPV <- function(parameters, seroData){

prev <- parameters

nPos <- prevData$A+prevData$C
nNeg <- prevData$B+prevData$D+prevData$E

lnLike <- sum(na.omit(nPos*log(prev) + nNeg*log(1-prev))) # Binomial log-likelihood

return(lnLike)
}


###################
## MCMC function ##
###################
ptm <- proc.time()[3]
for(i in 1:nIter){

param1 <- drop(rmvnorm(1, mean=param0, sigma=disp.prop))
###################
for(j in 1:length(prior[,1])){
param1[j] <- ifelse(param1[j] > prior[j,2], prior[j,2]-((param1[j]-prior[j,2])%%(prior[j,2]-prior[j,1])), ifelse(param1[j] < prior[j,1], prior[j,1]+((prior[j,1]-param1[j])%%(prior[j,2]-prior[j,1])), param1[j]))
}
###########################
lnLike1 <- logLikeHPV(param1, prevData)

if(!is.na(lnLike1)){
logAlpha <- (lnLike1 + dmvnorm(param0, param1, disp.prop, log=T)) - (lnLike0 + dmvnorm(param1, param0, disp.prop, log=T)) 

alpha <- exp(logAlpha)
if(alpha > 1) alpha <- 1
}

else alpha <- 0

u <- runif(1,0,1)
if(u<alpha){

	param0 <- param1
	lnLike0 <- lnLike1
	acc <- acc+1
}

acc.rate <-  (acc*100)/i
cat("\nAcceptance rate at jump", i, "is", round(acc.rate, digits=2) ,"%\n")

mcmc.out[i,1:nEstimate] <- param0
mcmc.out[i,(nEstimate+1)] <- lnLike0
}

cat("\nAlgorithm required", (proc.time()[3]-ptm)/60, "minutes\nwith final acceptance rate", round(acc.rate, digits=2), "%\n")

mcmc.out <- as.mcmc(mcmc.out) 
HPV.res <- window(mcmc.out, start=1001)
#plot(HPV.res)
ESS <- effectiveSize(HPV.res)
sampleOutput <- HPV.res [sample(nrow(HPV.res ), 25, replace=TRUE),1:29]



































