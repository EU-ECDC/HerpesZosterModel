source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
#source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(ggplot2)
theme_set(theme_classic())
library(gridExtra)

# Plot
get_results <- function(code, ...){
  get_data(code)
  #source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/MCMC.r")
  source("MCMC.r") # Use the local version
}
# TODO optimise such that results are only retrieved once - save in a list perhaps
# Should set a seed also?

library(parallel)

plot_list <- mclapply(use, function(code){get_results(code)
  ggplot() +
    geom_line(data = summaryFoI, mapping = aes(x = age, y = midFoi),
              size = 0.01, alpha = 0.8) +
    geom_ribbon(data = summaryFoI, mapping = aes(x = age, ymin = lower, ymax = upper),
                fill = "blue", alpha = 0.2) +
    geom_line(data = summaryPrev, mapping = aes(x = age, y = midPrev),
              size = 0.01, alpha = 0.8) +
    geom_ribbon(data = summaryPrev, mapping = aes(x = age, ymin = lower, ymax = upper),
                fill = "blue", alpha = 0.2) +
    ylim(0, 1) +
    geom_point(data = as.data.frame(cbind(age = as.numeric(row.names(htab)), 
                                          prop = htab[, 2] / rowSums(htab),
                                          tot = rowSums(htab))),
               mapping = aes(x = age, y = prop, size = tot),
               pch = 1) +
    labs(y = "", size = "number of samples", title = code) +
    scale_size_continuous(limits = c(1, 500),
                          breaks=c(50, 100, 150, 200,
                                   250, 300, 400, 500)) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(family = "serif"))})
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

# Save plot
tiff("S:/HelenJohnson/Herpes Zoster/Figures/overview_all.tif",
     width = 1400, height = 800)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IE, IT, LU, NL, SK, UK, RS, SI)
while(!is.null(dev.list())) dev.off()

# Plot rates
# First combine the values for all the countries
tmp <- sapply(use, function(code){
  get_results(code)
  return(as.data.frame(cbind(sampledR, ID = rep(code, dim(sampledR)[1]))))})
# Rework data into a format that is useful for plotting
tmp <- do.call(rbind, split(tmp, rep(1 : ncol(tmp), each = nrow(tmp))))
library(tidyr)
dat <- unnest(as.data.frame(tmp))
colnames(dat) <- c("R", "R0", "country")

library(tidybayes)
library(ggjoy)
# Save plots of the rates
tiff("S:/HelenJohnson/Herpes Zoster/Figures/rates_R.tif",
     width = 800, height = 1700)
if(FALSE){ # Currently broken
  as.mcmc.list(lapply(dat, mcmc)) %>% # R
    spread_draws(R[country]) %>%
    ggplot(aes(x = country, y = R)) +
    geom_halfeyeh() 
}
ggplot() + 
  geom_joy(data = dat, 
           mapping = aes(x = R, y = country, height = ..density..),
           scale = 0.85) +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(family = "serif"))
while(!is.null(dev.list())) dev.off()

tiff("S:/HelenJohnson/Herpes Zoster/Figures/rates_R0.tif",
     width = 800, height = 1700)
ggplot() + 
  geom_joy(data = dat, 
           mapping = aes(x = R0, y = country, height = ..density..),
           scale = 0.85) +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(family = "serif"))
while(!is.null(dev.list())) dev.off()

# Save values
save_vals <- function(code){
  get_results(code)
  return(c(code,
           propFac,
           prior,
           dispProp,
           nIter,
           acc.rate,
           ESS))
}

library(XLConnect)
sapply(use, function(code){
  mcmc_save <- loadWorkbook("MCMC settings.xlsx")
  appendWorksheet(mcmc_save, t(save_vals(code)), sheet = 1)
  saveWorkbook(mcmc_save)
})