# Clear workspace
#rm(list=ls())

# Required packages
#library(tidyverse)
#library(socialmixr) # POLYMOD data
#library(eurostat)   # EUROSTAT data

####################
## Initialisation ##
####################

seroData <- sero # remove this when Maria has finished editing load_data

## Group fixed parameters 
otherParams <- list(rij = contact_w,     # country-specific contact matrix
					N <- sum(popSize),   # population size
					D <- 6 / 365,        # duration of infection?!
					A <- 0.5,            # duration of maternally-derived immunity
					propFac <- "loglin", # type of proportionality factor
					Lmax <- 80           # life expectancy
					)
 
## Initialise MCMC params
nEstimate <- 1               # no. of parameters to be estimated 
nIter <- 1e4                 # number of iterations per chain
prior <- c(1e-6, 1e3)        # set uniform prior for q params
dispProp <-  1e-6 * prior[2] # dispersal parameter (inc. scaling coefficient)
acc <- 0                     # Initialise acceptance

mcmc.out <- rep(NA, (nEstimate+1)*nIter)  # Initialise output
dim(mcmc.out) <- c(nIter,(nEstimate+1))

######################
## Initialise chain ##
######################

param0 <- c(0.5, 0.3) # initial values for q parameters
lnLike0 <- FoI(param0, otherParams, seroData)$lnLike # Initialise log-likelihood

#####################################
## Estimate force of infection     ##
## and calculate associated R0, R, ##
## prevalence and log-likelihood   ##
#####################################

FoI <- function(fitParams, otherParams, seroData){

# Prepare serological data
age     <- seroData$AGE
obsData <- seroData$indic

## Generate transmission matrix
    if(propFac == "loglin"){ # First option (age-dependent susceptibility)
     
		# Check start parameters fit with prop. factor type
		if(length(param0) != 2){
			stop("length of param0 does not fit choice of proportionality factor")
		}
	  
		# Create transmission matrix with log-linear parameters
		qFunction <- function(x, obsData){exp(fitParams[1] + fitParams[2] * x)}
		qij <- outer(c(1 : Lmax), c(1 : Lmax), qFunction)
		bij <- 365 * qij * rij)[1 : Lmax, 1 : Lmax]
		}
	if(sum(is.na(bij)) > 0){warning("There are missing values in the transmission matrix beta_ij, which may break the tolerance parameter in the optimisation")}

## Initialise for iterative estimation of force of infection
    foiiprev <- rep(0.01, Lmax)
    tol <- 1 # initialise tolerance
    it <- 0 # iterator

## Iterate to estimate force of infection
    while((tol > 1e-10) & (it < 2000)){
	# Social contact FOI
	# lambda(a) is calculated based on beta(a, a') and I(a', t)
		foii <- (N / L) * D * bij %*%
			(as.matrix(foiiprev / (foiiprev + muy)) *
				matrix(c(1 - exp(- (1 - A) * (foiiprev[1] + muy[1])),
                exp(- (1 - A) * (foiiprev[1] + muy[1]) -
					c(0, cumsum(foiiprev[- 1] + muy[- 1])[1 : (Lmax - 2)])) -
                exp(- (1 - A) * (foiiprev[1] + muy[1]) - c(cumsum(foiiprev[- 1] + muy[- 1])))), ncol = 1))
		
		foii <- apply(cbind(0, foii), 1, max) # Bigger than 0
		foii <- apply(cbind(1, foii), 1, min) # Smaller than 1
		tol <- sum((foii - foiiprev) ^ 2) # Tolerance is sum of squared difference
		it <- it + 1
		foiiprev <- foii
    }
	
## Initialise modelled prevalence
	prev <- rep(NA, length(age))
	
## Estimate seroprevalence from fitted FoI model
	for(i in 1 : length(age)){
		# pi(a) = 1 - f(a)
		prev[i] <- (1 - exp(c(- (age[i] - A) * foii[1],
                            - (1 - A) * foii[1] - cumsum(c(0, foii[- 1])) -
                            (foii[- 1])[floor(age[i])] * (age[i] - floor(age[i])))))[floor(age[i]) + 1]
    }
	
## Estimate R0 from next-generation matrix
    R0ij <- (N / L) * D * bij[1 : Lmax, 1 : Lmax]
    Mij <- diag(c(My[1 : Lmax])) # g(a) (see notes)
    NGM <- Mij %*% R0ij
    NGM[is.na(NGM)] <- 0 # Replace missings with zeros
    R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
	R0 <- max(as.double(R0vec)) # spectral radius
	
## Estimate R
    seroTab <- table(floor(age), obsData) # Table with seroprevalence for each single year of age
	suscp <- seroTab[, 1] / rowSums(seroTab)
    # len included to reflect situations where NGM might be smaller than the size
    # of seroTab
    len <- min(dim(NGM)[1], length(suscp))
    Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
	R <- max(as.double(Rvec)) # spectral radius
	
## Calculate log-likelihood
	ll[i] <- obsData[i] * log(prev[i] + 1e-8) + (1 - obsData[i]) * log(1 - prev[i] + 1e-8) # This is currently done piecewise. Could use seroTab for grouping.
    
	return(list(lnLike = ll,     # log-likelihood
			    R0 = R0,         # R0
				R = R,           # R
				prev = prev,     # estimated seroprevalence
				foi = foiiprev)) # estimated force of infection
}

###################
## MCMC function ##
###################
ptm <- proc.time()[3]
for(i in 1:nIter){

	param1 <- drop(rmvnorm(1, mean=param0, sigma=dispProp))

	## Check what this is doing? Bouncing?
	for(j in 1:length(prior[,1])){
		param1[j] <- ifelse(param1[j] > prior[j,2], prior[j,2]-((param1[j]-prior[j,2])%%(prior[j,2]-prior[j,1])), ifelse(param1[j] < prior[j,1], prior[j,1]+((prior[j,1]-param1[j])%%(prior[j,2]-prior[j,1])), param1[j]))
	}

	lnLike1 <- FoI(param1, otherParams, seroData)$ll

	if(!is.na(lnLike1)){
		logAlpha <- (lnLike1 + dmvnorm(param0, param1, dispProp, log=T)) - (lnLike0 + dmvnorm(param1, param0, dispProp, log=T)) 

		alpha <- exp(logAlpha)
		if(alpha > 1) alpha <- 1
	}	

	else alpha <- 0

	u <- runif(1,0,1)
	if(u < alpha){
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



































