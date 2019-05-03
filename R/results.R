source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(ggplot2)
theme_set(theme_classic() %+replace%
            theme(plot.title = element_text(hjust = 0.5))) # Ensure centred titles
library(gridExtra)

# Plot
plot_results <- function(code, ...){
  get_data(code)
  source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/MCMC.r")
  p + labs(title = code)
}
library(parallel)
plot_list <- mclapply(use, plot_results)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

# Save plot
tiff("S:/HelenJohnson/Herpes Zoster/Figures/overview_all.tif",
     width = 1400, height = 800)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IE, IT, LU, NL, SK, UK, RS, SI)
while(!is.null(dev.list())) dev.off()
