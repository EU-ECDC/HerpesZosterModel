library(mvtnorm)
library(coda)
library(purrrlyr) # included in tidyverse?

####################
## Initialisation ##
####################

## Helen, tidy this up!
rij <- contact_w     # country-specific contact matrix
					N <- sum(popSize)   # population size
					D <- 6 / 365        # duration of infection?!
					A <- 0.5            # duration of maternally-derived immunity
					propFac <- "constant" # type of proportionality factor
					Lmax <- 80 
					
muy <- predict(demfit, type = "response")
muy <- muy[1 : Lmax] # Ensure not longer than life expectancy
  My <- exp(- cumsum(muy)) # Type I mortality
  L <- Lmax * mean(My) # Life expectancy
  My <- My[1 : Lmax] # Ensure no longer than life expectancy
 
 
## Group fixed parameters 
otherParams <- list(rij = contact_w,     # country-specific contact matrix
					N = nrow(seroData),  # sum(popSize),   # population size
					D = 6 / 365,        # duration of infection?!
					A = 0.5,            # duration of maternally-derived immunity
					propFac = "constant", # type of proportionality factor
					Lmax = 80           # life expectancy
					)
 
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
    if(propFac == "constant"){
      if(length(param0) != 1){
        stop("length of param0 does not fit choice of constant proportionality factor")
      }
      # Create transmission matrix with constant parameter
      bij <- 365 * fitParams[1] * (rij)[1 : Lmax, 1 : Lmax]
    }
	
    if(propFac == "loglin"){ # First option (age-dependent susceptibility)
     
		# Check start parameters fit with prop. factor type
		if(length(param0) != 2){
			stop("length of param0 does not fit choice of proportionality factor")
		}
	  
		# Create transmission matrix with log-linear parameters
		qFunction <- function(x, obsData){exp(fitParams[1] + (fitParams[2] * x))}
		qij <- outer(c(1 : Lmax), c(1 : Lmax), qFunction)
		bij <- 365 * qij * rij[1 : Lmax, 1 : Lmax]
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
	ll <- rep(NA, length(age)) # initialise log-like
	
## Estimate seroprevalence from fitted FoI model
	for(k in 1 : length(age)){
        prev[k] <- (1 - exp(c(- (age[k] - A) * foii[1],
                      - (1 - A) * foii[1] - cumsum(c(0, foii[- 1]))
					  - (foii[- 1])[floor(age[k])] * (age[k] - floor(age[k])))))[floor(age[k]) + 1]
      # Bernoulli loglikelihood used to estimate FOI
      ll[k] <- obsData[k] * log(prev[k]) + (1 - obsData[k]) * log(1 - prev[k])
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

    ll[!is.finite(ll)] <- 0 
	return(list(lnLike =  2 * sum(ll),     # log-likelihood
			    R0 = R0,         # R0
				R = R,           # R
				prev = prev,     # estimated seroprevalence
				foi = foiiprev)) # estimated force of infection
				

}

###################
## MCMC function ##
###################

## Initialise MCMC params(
nEstimate <- 1                  # no. of parameters to be estimated 
nIter <- 2e3                    # number of iterations per chain
prior <- c(0,0.2) #(c(0,1), c(-1,1)) # set uniform prior for q params
dim(prior) <- c(2,nEstimate)
prior <- t(prior) 
scaleDisp <- 1e-3				# set scaling factor for dispersion           
dispProp <-  diag(nEstimate) * scaleDisp * (prior[,2]-prior[,1]) # dispersal parameter
acc <- 0                        # Initialise acceptance

mcmcOutput <- rep(NA, (nEstimate+1)*nIter)  # Initialise output
dim(mcmcOutput) <- c(nIter,(nEstimate+1))

######################
## Initialise chain ##
######################
param0 <- 0.05
#param0 <- c(0.5, 0.5) # initial values for q parameters
lnLike0 <- FoI(param0, otherParams, seroData)$lnLike # Initialise log-likelihood

ptm <- proc.time()[3] # set clock
for(i in 1:nIter){

	param1 <- drop(rmvnorm(1, mean=param0, sigma=dispProp))

	# Ensure that sampled parameters fall within prior. If not, then bounce off boundary
	for(j in 1:length(prior[,1])){
		param1[j] <- ifelse(param1[j] > prior[j,2], prior[j,2]-((param1[j]-prior[j,2])%%(prior[j,2]-prior[j,1])), ifelse(param1[j] < prior[j,1], prior[j,1]+((prior[j,1]-param1[j])%%(prior[j,2]-prior[j,1])), param1[j]))
	}

	lnLike1 <- FoI(param1, otherParams, seroData)$lnLike

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

	mcmcOutput[i,1:nEstimate] <- param0
	mcmcOutput[i,(nEstimate+1)] <- lnLike0
}

cat("\nAlgorithm required", (proc.time()[3]-ptm)/60, "minutes\nwith final acceptance rate", round(acc.rate, digits=2), "%\n")

mcmcOutput <- as.mcmc(mcmcOutput) 
mcmcResults <- window(mcmcOutput, start=1) # 
plot(mcmcResults)
ESS <- effectiveSize(mcmcResults)

mcmcOutput <- as_tibble(mcmcOutput)
burnIn <- 100
sampleSize <- 250
postSample <- mcmcOutput %>% tail(-burnIn) %>%
				sample_n(sampleSize, replace=TRUE) %>%
				select(V1) %>%
				rename(gamma0 = V1) %>%
				mutate(q1 = exp(gamma0))
	
# Generate results for posterior samples	
sampledResults <- as.vector(postSample$gamma0) %>%
                  map(FoI, otherParams, seroData)

# R and R0 values from sampled parameters
sampledR <- sampledResults %>% {
     				tibble(
						R = map_dbl(., "R"),
					   R0 = map_dbl(., "R0")
					   )
					}

# Re-format force of infection estimates from sampled parameters
sampledFoI <- map_dfc(sampledResults, extract, "foi")
sampledFoI <- as_tibble(sampledFoI) %>% # combine with age data and reduce to unique values (i.e. not list for whole serosurvey sample)
					mutate(age = seq(0.5, Lmax-0.5)) %>% # check this. Should it indeed be mid-points?
					gather(key=foiName, value=foi, -age) %>%
                    mutate(set = sapply(strsplit(foiName, split='i', fixed=TRUE),function(x) (x[2]))) %>% # derive set name
					mutate(set = replace(set, is.na(set), 0)) %>% # rename the first set as 0
					mutate(set = as.integer(set) + 1) %>% # to avoid 0/1 confusion
					select(set, age, foi) # re-order for clarity

width1 <- 0.95

summaryFoI <- sampledFoI %>% 
				group_by(age) %>% 
					summarise(n = n(),
						  midFoi = median(foi),
						  lower = quantile(foi, probs = (1-width1)),
						  upper = quantile(foi, probs = width1)) %>%
				mutate(width = width1) %>%
				mutate(point = "median")

# Plot force of infection
ggplot(summaryFoI, aes(x=age, y=midFoi)) +
		geom_line(size=0.01, alpha=0.8) +
		geom_ribbon(aes(ymin=lower, ymax=upper) ,fill="blue", alpha=0.2) +
		ylim(0,1)
   

# Re-format prevalence estimates from sampled parameters
sampledPrev <- map_dfc(sampledResults, extract, "prev")
sampledPrev <- as_tibble(unique(cbind(seroData$AGE, sampledPrev))) %>% # combine with age data and reduce to unique values (i.e. not list for whole serosurvey sample)
					gather(key=prevName, value=prev, -`seroData$AGE`) %>%
                    mutate(set = sapply(strsplit(prevName, split='v', fixed=TRUE),function(x) (x[2]))) %>% # derive set name
					mutate(set = replace(set, is.na(set), 0)) %>% # rename the first set as 0
					mutate(set = as.integer(set) + 1) %>% # to avoid 0/1 confusion
					rename(age = 'seroData$AGE') %>%
					select(set, age, prev) # re-order for clarity

summaryPrev <- sampledPrev %>% 
				group_by(age) %>% 
					summarise(n = n(),
					midPrev = median(prev),
					lower = quantile(prev, probs = (1-width1)),
					upper = quantile(prev, probs = width1)) %>%
				mutate(width = width1) %>%
				mutate(point = "median")
				
# Plot prevalence
ggplot(summaryPrev, aes(x=age, y=midPrev)) +
		geom_line(size=0.01, alpha=0.8) +
		geom_ribbon(aes(ymin=lower, ymax=upper) ,fill="blue", alpha=0.2) +
		ylim(0,1)
		




