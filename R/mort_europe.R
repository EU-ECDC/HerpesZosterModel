# Prerequisites

library(eurostat)
library(tidylog)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(here)
source("https://raw.githubusercontent.com/EU-ECDC/EcdcSupportSnippets/master/R/eurostat_wrappers.R")
# source("https://raw.githubusercontent.com/EU-ECDC/EcdcSupportSnippets/master/R/save_outputs.R")

# ------------------------------------------------------------------------------

# Load data
countries_of_interest <- c("AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL",
                           "ES", "FI", "FR", "HR", "HU", "IE", "IS", "IT", "LI",
                           "LT", "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO",
                           "SE", "SI", "SK", "UK")

pop_data <- get_eurostat(id = "demo_pjan")
mort_data <- get_eurostat(id = "demo_magec")

use <- countries_of_interest[countries_of_interest %in% intersect(pop_data$geo, mort_data$geo)]
# Remove LI temporarily as it's causing errors
use <- setdiff(use, "LI")

# Fit mortality

fit_mort <- function(country, time = "2003-01-01", ...){
  PS <- pop_data %>% 
    filter(geo == country) %>%
    filter(sex == "T") %>%
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
    filter(time == time)
  AGE <- droplevels(PS$age)
  levels(AGE)[length(levels(AGE))] <- "0.5"
  AGE <- gsub("[A-z]","", AGE)
  AGE <- as.numeric(AGE)
  AGE <- sort(AGE, index.return = TRUE)$`x`
  PS <- PS[sort(AGE, index.return = TRUE)$ix, ]$values
  ND <- mort_data %>% 
    filter(geo == country) %>%
    filter(sex == "T") %>%
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
    filter(time == time)
  ND <- ND[sort(AGE, index.return = TRUE)$ix, ]$values
  # Fit contact data
  demfit <- mgcv::gam(ND ~ s(AGE), offset = log(PS), family = "poisson", link = "log")
  vals <- plot(demfit)
  return(as.data.frame(vals[[1]][c("x", "se", "fit")]))
}

plot_mort <- function(country, ...){
  ggplot(mapping = aes(x = x, y = fit), data = fit_mort(country)) + 
    labs(x = "age", y = "f(age)", title = country) +
    geom_line() +
    geom_ribbon(aes(ymin = fit - se,
                    ymax = fit + se),
                alpha = 0.2) + 
    theme(plot.title = element_text(hjust = 0.5))
}

# ------------------------------------------------------------------------------

# Plot
plot_list <- lapply(use, plot_mort)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = here("Figures/mort.tif"), width = 800, height = 600)
grid.arrange(AT, BE, BG, CY, CZ, DE, DK, EE, EL, ES, FI, FR, HR, HU,
             IE, IS, IT, #LI,
             LT, LU, LV, MT, NL, NO, PL, PT, RO, SE, SI, SK, UK)
while(!is.null(dev.list())) dev.off()
