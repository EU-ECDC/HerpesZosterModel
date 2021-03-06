# Updated code from original by Hens et al 2012 (function called contact)
# available from
# http://www.simid.be/modeling-infectious-disease-parameters-based-on-serological-and-social-contact-data/

FOI <- function(age, y, rij, muy, N, Dur, Lmax, A, startpar,
                prop = "constant", print = 0){
  muy <- muy[1 : Lmax] # Ensure not longer than life expectancy
  My <- exp(- cumsum(muy)) # Type I mortality
  L <- Lmax * mean(My) # Life expectancy
  My <- My[1 : Lmax] # Ensure no longer than life expectancy
  htab <- table(floor(age), y) # Table with seroprevalence for each single year of age
  qproc <- function(age, y, qpar, rij, Lmax, N, Dur, prop_fac = prop){
    # qpar is the parameter we optimise wrt, its starting value is startpar
    if(Lmax > 100){return("Please specify Lmax<100")} 
    # Lmax restriction is inherited from original code
    if(!(prop_fac %in% c("constant", "loglin", "ext_loglin", "logpoly")))
      stop("prop_fac must be one of constant, loglin, ext_loglin, or logpoly")
    if(prop_fac == "constant"){
      if(length(startpar) != 1){
        stop("length of startpar does not fit choice of proportionality factor")
      }
      # Create transmission matrix with constant parameter
      bij <- 365 * qpar * (rij)[1 : Lmax, 1 : Lmax]
    }
    # Below are other options for the social contact hypothesis
    if(prop_fac == "loglin"){
      if(length(startpar) != 2){
        stop("length of startpar does not fit choice of proportionality factor")
      }
      # Create transmission matrix with log-linear parameters
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
    }
    if(prop_fac == "ext_loglin"){
      warning("NB large values may produce an exponential value of Inf")
      if(length(startpar) != 2){
        stop("length of startpar does not fit choice of proportionality factor")
      }
      # Create transmission matrix with extended log-linear parameters
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x + qpar[2] * y)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
    }
    if(prop_fac == "logpoly"){
      if(length(startpar) != 3){
        stop("length of startpar does not fit choice of proportionality factor")
      }
      # Create transmission matrix with log-polynomial parameters
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x + qpar[3] * x * x)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
    }
    if(prop_fac == "ext_logpoly"){
      if(length(startpar) != 3){
        stop("length of startpar does not fit choice of proportionality factor")
      }
      # Create transmission matrix with extended log-polynomial parameters
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x + qpar[2] * y + qpar[3] * y * x)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
    }
    if(sum(is.na(bij)) > 0){warning("There are missing values in the transmission matrix beta_ij, which may break the tolerance parameter in the optimisation")}
    foiiprev <- rep(0.01, Lmax)
    tol <- 1 # tolerance level
    it <- 0 # iterator
    while((tol > 1e-10) & (it < 2000)){
      # Social contact FOI
      # lambda(a) is calculated based on beta(a, a') and I(a', t)
      foii <- (N / L) * Dur * bij %*%
        (as.matrix(foiiprev / (foiiprev + muy)) *
           matrix(c(1 - exp(- (1 - A) * (foiiprev[1] + muy[1])),
                    exp(- (1 - A) * (foiiprev[1] + muy[1]) -
                          c(0, cumsum(foiiprev[- 1] + muy[- 1])[1 : (Lmax - 2)])) -
                      exp(- (1 - A) * (foiiprev[1] + muy[1]) - c(cumsum(foiiprev[- 1] + muy[- 1])))), ncol = 1))
      foii <- apply(cbind(0, foii), 1, max) # Bigger than 0
      foii <- apply(cbind(1, foii), 1, min) # Smaller than 1
      tol <- sum((foii - foiiprev) ^ 2)
      it <- it + 1
      foiiprev <- foii
    }
    prev <- rep(NA, length(age))
    ll <- rep(NA, length(age))
    for(i in 1 : length(age)){
      # Each value is repeated by the number sampled in that age group
      # pi(a) = 1 - f(a)
      prev[i] <- (1 - exp(c(- (age[i] - A) * foii[1],
                            - (1 - A) * foii[1] - cumsum(c(0, foii[- 1])) -
                              (foii[- 1])[floor(age[i])] * (age[i] - floor(age[i])))))[floor(age[i]) + 1]
      # Bernoulli loglikelihood used to estimate FOI
      ll[i] <- y[i] * log(prev[i] + 1e-8) + (1 - y[i]) * log(1 - prev[i] + 1e-8)
    }
    R0ij <- (N / L) * Dur * bij[1 : Lmax, 1 : Lmax]
    Mij <- diag(c(My[1 : Lmax])) # g(a) (see notes)
    NGM <- Mij %*% R0ij
    NGM[is.na(NGM)] <- 0 # Replace missings with zeros
    R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
    # We use the proportion of susceptibles to obtain Rhat
    suscp <- htab[, 1] / rowSums(htab)
    # len included to reflect situations where NGM might be smaller than the size
    # of htab
    len <- min(dim(NGM)[1], length(suscp))
    Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
    return(list(ll = - 2 * sum(ll), eivalues = R0vec, bij = bij,
                foi = foiiprev, sero = prev, Rvalues = Rvec))
  }
  # Minimise the above with respect to the loglikelihood
  q.result <- nlm(function(qpar){return(qproc(age, y, qpar, rij, Lmax, N, Dur)$ll)},
                  startpar, print.level = print)
  # Return the values obtained for the value that minimises the above
  result.global <- qproc(age = age, y = y, q = q.result$estimate, rij = rij,
                         Lmax = Lmax, N = N, Dur = Dur)
  #if(isTRUE(all.equal(q.result$estimate, startpar))){warning("Optimisation did not start")}
  return(list(qhat = q.result$estimate, deviance = q.result$minimum, 
              aic = q.result$minimum + 2,
              bic = q.result$minimum + log(length(y)), bij = result.global$bij,
              R0 = max(as.double(result.global$eivalues)), # Spectral radius
              lambda = result.global$foi,
              R = max(as.double(result.global$Rvalues)), # Spectral radius
              pi = result.global$sero, inputs = list(start = startpar, prop = prop,
              age = age, y = y, rij = rij, muy = muy, N = N, Dur = Dur,
              Lmax = Lmax, A = A, print = print),
              iterations = q.result$iterations))
}

# # Updates include
## - Improved whitespace for readability;
## - Replaced the obsolete as.real with as.double
## - Incorporated contact.fitter.loglinear and contact.fitter through a type argument;
## - Included duration of maternal immunity A as argument
## - Moved subsetting of muy by Lmax into the function
## - Returned additional outputs
## - Used htab to calculate estimates of number of susceptibles
## - Included NGM expression prior to calculation of R0vec allowing us to include Rvec 
## and obtain  Rhat values from susceptibles
## - Replaced use of = as assignment operator with <-
## - Replaced T and F with TRUE and FALSE
## (Compare outputs from T <- FALSE or T <- 0 with TRUE <- FALSE to see why 
## using the shortened versions is undesirable)
## - Removed obsolete EISPACK argument from nlm optimisation
