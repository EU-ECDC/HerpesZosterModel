library(mvtnorm)
library(coda)
library(purrrlyr) # included in tidyverse?


myColours <- c("26 107 133", "241 214 118", "168 45 23")
ECDCcol <- sapply(strsplit(myColours, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))  # convert to hexadecimal

colours1 <- colorRampPalette(ECDCcol)(5)



####################
## Initialisation ##
####################

## Helen, tidy this up!
rij <- contact_w     # country-specific contact matrix
N <- sum(popSize)   # population size
D <- 6 / 365        # duration of infection?!
A <- 0.5            # duration of maternally-derived immunity
propFac <- "extloglin" # type of proportionality factor
Lmax <- 70 

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
					Lmax = 70           # life expectancy
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
		if(length(fitParams) != 1){
			stop("length of param0 does not fit choice of constant proportionality factor")
		}
		# Create transmission matrix with constant parameter
		qij <- exp(fitParams[1])
		bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
    }
	
    if(propFac == "loglin"){ # First option (age-dependent susceptibility)
		# Check start parameters fit with prop. factor type
		if(length(fitParams) != 2){
			stop("length of param0 does not fit choice of proportionality factor")
		}
	  
		# Create transmission matrix with log-linear parameters
		qFunction <- function(x, obsData){exp(fitParams[1] + (fitParams[2] * x))}
		qij <- outer(c(1 : Lmax), c(1 : Lmax), qFunction)
		bij <- 365 * qij * rij[1 : Lmax, 1 : Lmax]
	}
	
    if(propFac == "extloglin"){
      if(length(fitParams) != 3){
        stop("length of param0 does not fit choice of proportionality factor")
      }
      # Create transmission matrix with extended log-linear parameters
      q.f <- function(x, y){exp(fitParams[1] + fitParams[2] * x + fitParams[3] * y)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
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
	return(list(lnLike = sum(ll),  # log-likelihood
			    R0 = R0,         # R0
				R = R,           # R
				prev = prev,     # estimated seroprevalence
				foi = foiiprev)) # estimated force of infection
				

}

###################
## MCMC function ##
###################

## Initialise MCMC params
# no. of parameters to be estimated 
nEstimate <- case_when(
				propFac == "constant" ~ 1,
				propFac == "loglin" ~ 2,
				propFac == "extloglin" ~ 3
				) 
				
nIter <- 2e3 # number of iterations per chain

prior <- c(-2, -1, -0.1, 0, -0.5, 0) # set uniform prior for q params
#prior <- c(-5, 5)
dim(prior) <- c(2, nEstimate)
prior <- t(prior)

if(dim(prior)[1] / nEstimate != 1)
  warning("Prior length does not match desired propFac")

init.scale.sd  <- 1e-2	# set scaling factor for proposal distribution  
covmat.proposal <-  diag(nEstimate) * init.scale.sd  * (prior[,2]-prior[,1]) # proposal distribution
adapt.size.start <- 100
adapt.size.cooling <- 0.99
adapt.shape.start <- 75
adapt.shape.stop <- NULL
max.scaling.sd <- 50
   
acc <- 0  # Initialise acceptance
acc.rate <- 1 #Initalise acceptance rate		  

mcmcOutput <- rep(NA, (nEstimate+1)*nIter)  # Initialise output
dim(mcmcOutput) <- c(nIter,(nEstimate+1))


## Function for updating covariance matrix
updateCovmat <- function(covmat, param.mean, param0, i) {

    residual <- as.vector(param0-param.mean)
    covmat <- (covmat*(i-1)+(i-1)/i*residual%*%t(residual))/i
    param.mean <- param.mean + residual/i

    return(list(covmat=covmat,param.mean=param.mean))
}

######################
## Initialise chain ##
######################
#param0 <- -0.5
param0 <- c(-1.5, -0.05, -0.25) # initial values for q parameters
if(length(param0) / nEstimate != 1)
  warning("Parameter input length does not match desired propFac")
lnLike0 <- FoI(param0, otherParams, seroData)$lnLike # Initialise log-likelihood

# Initialise covariance matrix
covmat.proposal.init <- covmat.proposal

adapting.size <- FALSE # will be set to TRUE once we start adapting the size               
adapting.shape <- 0  # will be set to the iteration at which adaptation starts

# scaling factor for covmat size
scaling.sd  <- 1

# scaling multiplier
scaling.multiplier <- 1

# empirical covariance matrix (0 everywhere initially)
covmat.empirical <- covmat.proposal
covmat.empirical[,] <- 0

# empirical mean vector
param.mean <- param0

ptm <- proc.time()[3] # set clock
for(i.iteration in 1:nIter){

        # Adaptation of proposal distribution
        if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
           (is.null(adapt.shape.start) || acc.rate*i.iteration < adapt.shape.start)) {
            if (!adapting.size) {
                message("\n---> Start adapting size of covariance matrix")
                adapting.size <- TRUE
            }
            # adapt size of covmat until we get enough accepted jumps
            scaling.multiplier <- exp(adapt.size.cooling^(i.iteration-adapt.size.start) * (acc.rate - 0.234))
            scaling.sd <- scaling.sd * scaling.multiplier
            scaling.sd <- min(c(scaling.sd, max.scaling.sd))
            # only scale if it doesn't reduce the covariance matrix to 0
            covmat.proposal.new <- scaling.sd^2 * covmat.proposal.init
            if (!(any(diag(covmat.proposal.new) <
                .Machine$double.eps))) {
                covmat.proposal <- covmat.proposal.new
            }

        } else if (!is.null(adapt.shape.start) &&
                   acc.rate*i.iteration >= adapt.shape.start &&
                   (adapting.shape == 0 || is.null(adapt.shape.stop) || 
                    i.iteration < adapting.shape + adapt.shape.stop)) {
            if (!adapting.shape) {
                message("\n---> Start adapting shape of covariance matrix")
                # flush.console()
                adapting.shape <- i.iteration
            }

            ## adapt shape of covariance matrix using optimal scaling factor for multivariate target distributions
            scaling.sd <- 2.38/sqrt(nEstimate)

            covmat.proposal <- scaling.sd^2 * covmat.empirical
        } else if (adapting.shape > 0) {
            message("\n---> Stop adapting shape of covariance matrix")
            adapting.shape <- -1
        }


	param1 <- drop(rmvnorm(1, mean=param0, sigma=covmat.proposal))

	# Ensure that sampled parameters fall within prior. If not, then bounce off boundary
	for(j in 1:length(prior[,1])){
		param1[j] <- ifelse(param1[j] > prior[j,2], prior[j,2]-((param1[j]-prior[j,2])%%(prior[j,2]-prior[j,1])), ifelse(param1[j] < prior[j,1], prior[j,1]+((prior[j,1]-param1[j])%%(prior[j,2]-prior[j,1])), param1[j]))
	}

	lnLike1 <- FoI(param1, otherParams, seroData)$lnLike

	if(!is.na(lnLike1)){
		logAlpha <- (lnLike1 + dmvnorm(param0, param1, covmat.proposal, log=T)) - (lnLike0 + dmvnorm(param1, param0, covmat.proposal, log=T)) 

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

 # update empirical covariance matrix
    if (adapting.shape >= 0) {
        tmp <- updateCovmat(covmat.empirical, param.mean,
                            param0, i.iteration)
        covmat.empirical <- tmp$covmat
        param.mean <- tmp$param.mean
        }
		
	acc.rate <- acc/i.iteration
	cat("\nAcceptance rate at jump", i.iteration, "is", round(acc.rate, digits=3) ,"\n")

	mcmcOutput[i.iteration,1:nEstimate] <- param0
	mcmcOutput[i.iteration,(nEstimate+1)] <- lnLike0
}

cat("\nAlgorithm required", (proc.time()[3]-ptm)/60, "minutes\nwith final acceptance rate", round(acc.rate, digits=3), "\n")

burnIn <- 500

mcmcOutput <- as.mcmc(mcmcOutput) 
mcmcResults <- window(mcmcOutput, start=(1)) # 
plot(mcmcResults)
ESS <- effectiveSize(mcmcResults)

AIC <- (2*nEstimate) + 2 * max(mcmcOutput[,ncol(mcmcOutput)]) # Akaike information criterion

mcmcOutput <- as_tibble(mcmcOutput)
sampleSize <- 100

if(propFac == "constant"){
postSample <- mcmcOutput %>% tail(-burnIn) %>%
				sample_n(sampleSize, replace=TRUE) %>%
				select(V1) %>%
				rename(gamma0 = V1) %>%
				mutate(q1 = exp(gamma0))
	
# Generate results for posterior samples	
sampledResults <- as.vector(postSample$gamma0) %>%
                  map(FoI, otherParams, seroData)
				  
}

if(propFac == "loglin"){
postSample <- mcmcOutput %>% tail(-burnIn) %>%
				sample_n(sampleSize, replace=TRUE) %>%
				select(V1, V2) %>%
				rename(gamma0 = V1) %>%
				rename(gamma1 = V2) %>%
				mutate(q1 = exp(gamma0 + (gamma1 * mean(seq(0.5, Lmax-0.5)))))
	
# Generate results for posterior samples	
sampledResults <- mapply(c, postSample$gamma0, postSample$gamma1, SIMPLIFY = FALSE) %>% # list of pairs (gamma0, gamma1)
                  map(FoI, otherParams, seroData)
				  
}

if(propFac == "extloglin"){
postSample <- mcmcOutput %>% tail(-burnIn) %>%
				sample_n(sampleSize, replace=TRUE) %>%
				select(V1, V2, V3) %>%
				rename(gamma0 = V1) %>%
				rename(gamma1 = V2) %>%
				rename(gamma2 = V3) %>%
				mutate(q1 = exp(gamma0 + (gamma1 * mean(seq(0.5, Lmax-0.5))) + (gamma2 * mean(seq(0.5, Lmax-0.5)))))
	
# Generate results for posterior samples	
sampledResults <- mapply(c, postSample$gamma0, postSample$gamma1, postSample$gamma2, SIMPLIFY = FALSE) %>% # list of pairs (gamma0, gamma1, gamma2)
                  map(FoI, otherParams, seroData)
				  
}

# R and R0 values from sampled parameters
sampledR <- sampledResults %>% {
     				tibble(
						R = map_dbl(., "R"),
					   R0 = map_dbl(., "R0")
					   )
					}
										
sampledR <- bind_cols(postSample, sampledR)
#sampledR %>% filter(R >= 0.95 & R < 1.05) # filter for parameter sets close to endmeic eqm

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

# Plot force of infection and prevalence
htab <- table(
  floor(
    seroData[(seroData$AGE > 0.5) & 
               (seroData$AGE < 80) &
               (!is.na(seroData$AGE)) & 
               !is.na(seroData$indic), ]$AGE[order(seroData[(seroData$AGE > 0.5) & 
                                                              (seroData$AGE < 80) &
                                                              (!is.na(seroData$AGE)) & 
                                                              !is.na(seroData$indic), ]$AGE)]),
  seroData[(seroData$AGE > 0.5) &
             (seroData$AGE < 80) &
             (!is.na(seroData$AGE)) & 
             !is.na(seroData$indic), ]$indic[order(seroData[(seroData$AGE > 0.5) & 
                                                              (seroData$AGE < 80) &
                                                              (!is.na(seroData$AGE)) & 
                                                              !is.na(seroData$indic), ]$AGE)])

(p <- ggplot() +
    geom_line(data = summaryFoI, mapping = aes(x = age, y = midFoi),
              size = 0.01, alpha = 0.8) +
    geom_ribbon(data = summaryFoI, mapping = aes(x = age, ymin = lower, ymax = upper),
                fill = "blue", alpha = 0.2) +
    geom_line(data = summaryPrev, mapping = aes(x = age, y = midPrev),
              size = 0.01, alpha = 0.8) +
    geom_ribbon(data = summaryPrev, mapping = aes(x = age, ymin = lower, ymax = upper),
                fill = "dark green", alpha = 0.3) +
    ylim(0, 1) +
    geom_point(data = as.data.frame(cbind(age = as.numeric(row.names(htab)), 
                                          prop = htab[, 2] / rowSums(htab),
                                          tot = rowSums(htab))),
               mapping = aes(x = age, y = prop, size = tot),
               pch = 1) +
    labs(y = "", size = "number of samples"))
