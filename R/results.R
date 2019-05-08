library(XLConnect)
library(lubridate)

#calc <- function(code){
  
  get_data(code)
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

  prior <- c(-2, -1, -0.1, 0, -0.5, 0) # set uniform prior for q params
  #prior <- c(-5, 5)
  #param0 <- -0.5
  param0 <- c(-1.5, -0.05, -0.25) # initial values for q parameters
  
  source("MCMC.r")
  vals <- list(country = code, propFac = propFac, prior = prior, 
               init.scale.sd = init.scale.sd, nIter = nIter, acc.rate = acc.rate, 
               param0 = param0, mcmcOutput = mcmcOutput, postSample = postSample, 
               sampledResults = sampledResults, ESS = ESS, AIC = AIC, 
               date = today())
#  return(vals)
#}

code <- "BE"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resBE_extloglin_1.txt")
# Currently saving plots manually

code <- "DE"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resDE_extloglin_1.txt")

code <- "FI"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resFI_extloglin_1.txt")

code <- "IE"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resIE_extloglin_1.txt")

code <- "IT"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resIT_extloglin_1.txt")

code <- "LU"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resLU_extloglin_1.txt")

code <- "NL"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resNL_extloglin_1.txt")

code <- "SK"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resSK_extloglin_1.txt")

code <- "UK"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resUK_extloglin_1.txt")

code <- "RS"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resRS_extloglin_1.txt")

code <- "SI"
dput(vals, file = "S:/HelenJohnson/Herpes Zoster/resSI_extloglin_1.txt")
