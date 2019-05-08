# Load data and required packages
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")

library(ggplot2)
library(gridExtra)

####################
## Initialisation ##
####################
init <- function(code, D, A, propFac, Lmax, prior,
                 param0){
  get_data(code)
  otherParams <- list(rij = contact_w,     # country-specific contact matrix
                      N = nrow(seroData),  # sum(popSize),   # population size
                      D = D,        # duration of infection?!
                      A = A,            # duration of maternally-derived immunity
                      propFac = propFac, # type of proportionality factor
                      Lmax = Lmax,           # life expectancy
                      muy = predict(demfit, type = "response")[1 : Lmax],
                      L = Lmax * mean(exp(- cumsum(predict(demfit, type = "response")[1 : Lmax]))),
                      My = exp(- cumsum(predict(demfit, type = "response")[1 : Lmax]))[1 : Lmax])
  out <- list(otherParams = otherParams, prior = prior, param0 = param0,
              seroData = seroData)
  return(out)
}

# Load FoI and MCMC functions
source("MCMC.r")

# Example -- Belgium -- Example -- Belgium -- Example -- Belgium -- Example --
input <- init(code = "BE",
              D = 6 / 365,
              A = 0.5,
              Lmax = 70,
              prior = c(-2, 1, -0.1, 0, -0.5, 0),
              param0 = c(-1.5, -0.05, -0.25),
              propFac = "extloglin")

param0 <- input$param0
prior <- input$prior
otherParams <- input$otherParams
seroData <- input$seroData
propFac <- input$otherParams$propFac

# Run MCMC
res <- MCMC(nIter = 2e3, # number of iterations per chain
     init.scale.sd = 1e-2,	# set scaling factor for proposal distribution  
     adapt.size.start = 100,
     adapt.size.cooling = 0.99,
     adapt.shape.start = 75,
     adapt.shape.stop = NULL,
     max.scaling.sd = 50,
     burnIn = 500,
     sampleSize = 100,
     width1 = 0.95,
     Lmax = input$otherParams$Lmax)

# Calculate htab
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

# Plot results
p <- ggplot() +
  geom_line(data = res$summaryFoI, mapping = aes(x = age, y = midFoi),
            size = 0.01, alpha = 0.8) +
  geom_ribbon(data = res$summaryFoI, mapping = aes(x = age, ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.2) +
  geom_line(data = res$summaryPrev, mapping = aes(x = age, y = midPrev),
            size = 0.01, alpha = 0.8) +
  geom_ribbon(data = res$summaryPrev, mapping = aes(x = age, ymin = lower, ymax = upper),
              fill = "dark green", alpha = 0.3) + 
  labs(y = "", title = "BE") +
  theme(plot.title = element_text(hjust = 1,
                                  margin = margin(t = 10, b = - 20)),
        title = element_text(family = "serif"))
p
# These numbers seem to be too small ! They look very odd

# All countries
## Gather options used for the various countries
opts <- list(country = rep(use, each = 3),
               propFac = rep(c("constant", "loglin", "extloglin"), times = length(use)),
               D = rep(6 / 365, times = 3 * length(use)),
               A = rep(0.5, times = 3 * length(use)),
               Lmax = rep(70, times = 3 * length(use)),
               prior = list(constant = rep(c(-2, 1), times = length(use)),
                            loglin = rep(c(-2, 1, -0.1, 0), times = length(use)),
                            extloglin = rep(c(-2, 1, -0.1, 0, -0.5, 0), times = length(use))),
               param0 = list(constant = rep(c(-1.5), times = length(use)),
                             loglin = rep(c(-1.5, -0.05), times = length(use)),
                             extloglin = rep(c(-1.5, -0.05, -0.25), times = length(use))))
# Or have one list for each propFac?