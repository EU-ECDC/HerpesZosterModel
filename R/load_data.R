# Clear workspace
rm(list = ls())

# Load required packages
library(tidyverse)
library(readxl)
library(magrittr)
library(data.table)
library(socialmixr)
library(eurostat)
library(akima)
library(gridExtra)
library(ggplot2)
theme_set(theme_classic() %+replace%
            theme(plot.title = element_text(hjust = 0.5))) # Ensure centred titles

# Load data
## ESEN2 Seroprevalence
esen <- read_csv("S:/HelenJohnson/Herpes Zoster/Force of infection/esen2vzv.csv",
                  col_names = TRUE)

# Create indicator variable for seroprevalence. Code 'EQI' as positive
esen <- esen %>%
  mutate(indic = ifelse(STDRES != "NEG", 1, 0))

## AGE
# Remove ages with + as we are unsure how to include them
esen <- esen[- which(esen$AGE == "60+"), ]

# In the case of age ranges, choose mid-point
esen$AGE <- sapply(strsplit(esen$AGE, split = "-"),
                    function(x) mean(as.numeric(x)))
# Assume mid point of year of age
esen <- esen %>% 
  mutate(AGE = AGE + 0.5)

# EU, EEA, and EU candidate countries
countries <- tibble(code = c("AT", "BE", "BG", "CY", "CZ", "DE", "DK",
                             "EE", "EL", "ES", "FI", "FR", "HR", "HU",
                             "IE", "IS", "IT", "LI", "LT", "LU", "LV",
                             "MT", "NL", "NO", "PL", "PT", "RO", "SE",
                             "SI","SK", "UK", "AL", "ME", "MK", "RS",
                             "TR"),
                    name = c("Austria", "Belgium", "Bulgaria", "Cyprus", 
                             "Czechia", "Germany", "Denmark", "Estonia",
                             "Greece", "Spain", "Finland", "France", "Croatia",
                             "Hungary", "Ireland", "Iceland", "Italy", 
                             "Liechtenstein", "Lithuania", "Luxembourg", 
                             "Latvia", "Malta", "Netherlands", "Norway", 
                             "Poland", "Portugal", "Romania", "Sweden",
                             "Slovenia", "Slovakia", "United Kingdom", 
                             "Albania", "Montenegro", "North Macedonia",
                             "Serbia", "Turkey"))

# Get and assign the data needed -----------------------------------------------
get_data <- function(code){
  if(!(code %in% countries$code)){stop("Invalid value")}
  # Clear previous values
  remove(sero, envir = .GlobalEnv)
  remove(contact_w, envir = .GlobalEnv)
  remove(demfit, envir = .GlobalEnv)
  remove(popSize, envir = .GlobalEnv)
  
  # Population data
  pop <- get_eurostat(id = "demo_pjan") %>% # Population size
    filter(geo == code) %>% # Country of interest
    filter(sex == "T") %>% # Total population (males & females)
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>% # Include only population values by age
    filter(time == "2003-01-01") # Determined from date ESEN2 gathered
  
  # Mortality data
  mort <- get_eurostat(id = "demo_magec") %>% # Number of deaths
    filter(geo == code) %>% # Country of interest
    filter(sex == "T") %>%   # Total number of deaths (males & females)
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>% # Include only deaths by age
    filter(time == "2003-01-01") # Determined from date ESEN2 gathered
  
  # Re-order population data by age
  popAge <- droplevels(pop$age) # Drop unused levels
  levels(popAge)[length(levels(popAge))] <- "0.5" # Relabel Y_LT1 as 0.5
  popAge <- as.numeric(gsub("[A-z]", "", popAge)) # Replace letters (here Y) with nothing
  # and convert to numeric
  
  # Population sizes by age
  popSize <- pop %>% 
    mutate(popAge = popAge) %>%
    arrange(popAge) %>% 
    mutate(age = factor(age, age)) %>%
    select(values)
  popSize <- unlist(popSize)
  
  # Re-order mortality data by age
  mortAge <- droplevels(mort$age) # Drop unused levels
  levels(mortAge)[length(levels(mortAge))] <- "0.5" # Relabel Y_LT1 as 0.5
  mortAge <- as.numeric(gsub("[A-z]", "", mortAge)) # Replace letters with nothing
  # and convert to numeric
  
  # Number of deaths by age
  nDeaths <- mort %>% 
    mutate(mortAge = mortAge) %>%
    arrange(mortAge) %>% 
    mutate(age = factor(age, age)) %>%
    select(values)
  nDeaths <- unlist(nDeaths)
  
  len <- min(length(popSize), length(nDeaths))
  nDeaths <- nDeaths[c(1 : len)] # Ensure same dimensions for later calculation
  popSize <- popSize[c(1 : len)] # Ensure same dimensions for later calculation
  
  popAge <- sort(popAge) # Sort population age such that it is ascending
  mortAge <- sort(mortAge) # Sort mortality age such that it is ascending
  
  # Fit mortality model
  ## Mortality rates (offset by population size) are estimated through a Poisson GAM
  ## with a log link
  demfit <- mgcv::gam(nDeaths ~ s(mortAge[c(1 : len)]), 
                      offset = log(popSize),
                      family = "poisson", link = "log")
  
  # Country name from code included in list of POLYMOD countries
  if(countries$name[which(countries$code %in% code)] %in%
     unique(polymod$participants$country)){
    # Obtain contact matrix from POLYMOD 
    cont <- contact_matrix(polymod,
                           countries = countries$name[which(countries$code %in% code)],
                           # Contacts lasting more than 15 minutes
                           filter = ("phys_contact" > 3),
                           quiet = TRUE)$matrix
    # quiet = TRUE supresses message about citing POLYMOD package
    
    # Select seroprevalence data
    if(code == "UK"){ # UK is coded differently in this data set
      sero <- esen %>% 
        filter(COUNTRY == code)
    } else {
      sero <- esen %>% 
        filter(COUNTRY == countries$name[which(countries$code %in% code)])
    }
    
    # Weight contact intensities by population sizes to obtain contact rates
    pop_weights <- data.frame(age = popAge, size = popSize)[1 : dim(cont)[1], ]
    weigh <- function(x){
      cont[x, ] / pop_weights[x, 2]
      # Person doing the contacting is weighted
    }
    tmp <- lapply(1 : dim(cont)[1], weigh)
    contact_w <- do.call(rbind, tmp)
    dimnames(contact_w) <- dimnames(cont)
  }
  
  if(code == "RS"){
    # Create Serbia seroprevalence data based on Table 1 from Medić et al 2018
    # Epidemiol Infect 146, 1593–1601
    # https://doi.org/10.1017/S0950268818001619
    zeros <- c(25, 235, 135, 64, 42, 12, 12, 7, 8, 6, 2)
    ones <- c(75, 165, 378, 448, 533, 188, 188, 193, 192, 194, 198)
    indic <- c(rep(0, sum(zeros)), rep(1, sum(ones)))
    age_vals <- c("0", "1–4", "5–9", "10–14", "15–19", "20–24", "25–29",
                  "30–34", "35–39", "40–49", "50–59")
    age_vals <- c(0, (1 + 4) / 2, (5 + 9) / 2, (10 + 14) / 2, (15 + 19) / 2, 
                  (20 + 24) / 2, (25 + 29) / 2, (30 + 34) / 2, (35 + 39) / 2,
                  (40 + 49) / 2, (50 + 59) / 2)
    serb <- data.frame(COUNTRY = rep("Serbia", length(indic)),
                       AGE = c(rep(age_vals, zeros),
                               rep(age_vals, ones)),
                       indic)
    ## Replace age values as age 0 causes error in the model
    age_rpl <- data.frame(AGE = c(0, seq(from = 1, to = 20, by = 1)),
                          to = c((1 + 0.5)/2, seq(from = 1.5, to = 20.5, by = 1)))
    setDT(serb)
    setDT(age_rpl)
    sero <- serb[age_rpl, on = c("AGE"), AGE := to]
    
    # Load contact matrix
    load("S:/HelenJohnson/Herpes Zoster/Data/prem.Rda")
    # Weight contact intensities by population sizes to obtain contact rates
    weigh <- function(x, dat, mat){
      AGE <- droplevels(dat$age)
      levels(AGE)[length(levels(AGE))] <- "0.5"
      AGE <- gsub("[A-z]","", AGE)
      AGE <- as.numeric(AGE)
      # Ensure both follow same age pattern
      AGE <- sort(AGE, index.return = TRUE)$`x`
      dat <- dat[sort(AGE, index.return = TRUE)$ix, ]$values
      return(mat[x, ] / dat[x])
    }
    # Use same weighting when combining as Fumanelli et al.
    mat <- 0.3 * prem$RS$home + 0.18 * prem$RS$school +
      0.19 * prem$RS$work + 0.33 * prem$RS$other
    tmp <- lapply(1 : dim(mat)[1],
                  function(x){weigh(x, 
                                    dat = pop,
                                    mat = mat)})
    RS <- as.matrix(do.call(rbind, tmp))
    dimnames(RS)[[1]] <- dimnames(RS)[[2]]
    #contact_w <- RS
    # Interpolate contact matrix
    tmp <- RS
    tmp <- melt(tmp)
    tmp %>% mutate_if(is.factor, as.character) -> tmp # Turn characters into factors
    # Remove age ranges with +
    tmp <- tmp[- which(tmp$Var1 == "75+"), ]
    tmp <- tmp[- which(tmp$Var2 == "75+"), ]
    # Take mid range value of age
    tmp$Var1 <- sapply(strsplit(tmp$Var1, split = "-"),
                       function(x) mean(as.numeric(x)))
    tmp$Var2 <- sapply(strsplit(tmp$Var2, split = "-"),
                       function(x) mean(as.numeric(x)))
    
    contact_w <- interp(tmp$Var1, tmp$Var2, tmp$value,
                        nx = 74, ny = 74)$z
    colnames(contact_w) <- rownames(contact_w) <- seq(0, 73, 1)
  }
  if(code == "SI"){
    # Create Slovenia data based on Figure 1 from Soca et al 2010
    # BMC Public Health 10:360
    # https://doi.org/10.1186/1471-2458-10-360
    tmp <- read_excel("S:/HelenJohnson/Herpes Zoster/Force of infection/old/slovenia_seroprev.xlsx")
    # Convert columns to numeric
    tmp[, c(2, 3, 4)] %<>% lapply(function(x) as.numeric(x))
    
    tmp <- na.omit(tmp) # Remove excess rows added due to additional calculations
    # further down in the sheet
    
    # Calculate how many zeros and ones based on seropositivy rate
    tmp <- tmp %>% 
      mutate(ones = round(n * mid / 100),
             zeros = n - ones)
    
    # Create sero data set
    age_vals <- c(rep(tmp$age, tmp$zeros), rep(tmp$age, tmp$ones))
    indic <- c(rep(rep(0, length(tmp$age)), tmp$zeros), 
               rep(rep(1, length(tmp$age)), tmp$ones))
    sero <- data.frame(COUNTRY = rep("Slovenia", length(indic)),
                       AGE = age_vals,
                       indic)
    # Contact matrix - same function used as that for Serbia
    weigh <- function(x, dat, mat){
      AGE <- droplevels(dat$age)
      levels(AGE)[length(levels(AGE))] <- "0.5"
      AGE <- gsub("[A-z]","", AGE)
      AGE <- as.numeric(AGE)
      AGE <- sort(AGE, index.return = TRUE)$`x`
      dat <- dat[sort(AGE, index.return = TRUE)$ix, ]$values
      return(mat[x, ] / dat[x])
    }
    load("S:/HelenJohnson/Herpes Zoster/Data/fumanelli.Rda")
    mat <- fumanelli$SI$total
    tmp <- lapply(1 : dim(mat)[1],
                  function(x){weigh(x, 
                                    dat = pop,
                                    mat = mat)})
    SI <- as.matrix(do.call(rbind, tmp))
    dimnames(SI)[[1]] <- dimnames(SI)[[2]]
    contact_w <- SI
  }
  
  contact_w[is.na(contact_w)] <- 0 # Replace missings with zeros
  sero <- sero[!is.na(sero$indic) & !is.na(sero$AGE), ] # remove serology results where no indication or age specified
  sero <- sero[order(sero$AGE), ] # order serology results by age
  
  if(dim(pop)[1] == 0)
    warning("Population data for year 2003 not available from Eurostat database demo_pjan")
  if(dim(mort)[1] == 0)
    warning("Mortality data for year 2003 not available from Eurostat database demo_magec")
  
  # Save objects
  assign("sero", sero, envir = .GlobalEnv)
  assign("contact_w", contact_w, envir = .GlobalEnv)
  assign("demfit", demfit, envir = .GlobalEnv)
  assign("popSize", popSize, envir = .GlobalEnv)
}

# Countries for which we have seroprevalence data
use <- c(countries$code[countries$name %in% 
                          unique(esen$COUNTRY)], "UK", "RS", "SI")

# Plots of data ----------------------------------------------------------------
## Plot mortality
plot_mort <- function(code, ...){
  get_data(code)
  vals <- plot(demfit) # Save output from demfit to plot in ggplot2 style
  ggplot(mapping = aes(x = x, y = fit), 
         data = as.data.frame(vals[[1]][c("x", "se", "fit")])
         # Requried parts from vals saved as data frame
         ) + 
    labs(x = "age", y = "f(age)", title = code) +
    geom_line() +
    geom_ribbon(aes(ymin = fit - se,
                    ymax = fit + se),
                alpha = 0.2)  + 
    theme(plot.title = element_text(hjust = 0.5))
}

# Save list of plots for the countries in use
plot_list <- lapply(use, plot_mort)
# Save the plots as their country codes
lapply(seq_along(1 : length(use)), function(x){assign(use[x], plot_list[[x]], 
                                                      envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/mortality.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
# TODO see if we can replace ^ with use somehow
while(!is.null(dev.list())) dev.off()

## Plot serological data
plot_sero <- function(code, ...){
  get_data(code)
  sero <- sero[(sero$AGE > 0.5) & (sero$AGE < 80) &
    (!is.na(sero$AGE)) & !is.na(sero$indic), ]
  y <- sero$indic[order(sero$AGE)]
  a <- sero$AGE[order(sero$AGE)]
  
  # Calculate seropositivity rates
  grid <- sort(unique(round(a)))
  neg <- table(y, round(a))[1, ]
  pos <- table(y, round(a))[2, ]
  tot <- neg + pos
  
  ggplot(data = as.data.frame(cbind(grid, neg, pos, tot)),
         mapping = aes(x = grid, y = pos / tot, size = 0.02 * tot)) +
    geom_point(pch = 1) + 
    labs(x = "age", y = "sero-prevalence", title = code) + 
    xlim(0, 72) + ylim(- 0.1, 1) + theme(legend.title = element_blank())
}

# Save list of plots and assign to the environment
plot_list <- lapply(use, plot_sero)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/serology.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
while(!is.null(dev.list())) dev.off()

## Plot contact rates
# Adapted from https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/old_scripts/plots/plot_polymod_matrix.R
plot_mat <- function(code, ...){
  get_data(code)
  tmp <- melt(contact_w)
  p <- ggplot(data = tmp, aes_string(x = names(tmp)[1], y = names(tmp)[2], 
                                fill = names(tmp)[3])) + 
    geom_tile() + 
    scale_fill_viridis_c(direction = - 1, option = "E", ...) +
    labs(x = "", y = "", title = code) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # Rotate x-axis labels
    if(is.factor(tmp[, 1])){ # Add spacing
      scale_x_discrete(breaks = levels(tmp[, 1])[c(TRUE, rep(FALSE, 9))])
    } else {
      NULL
    }
  p <- p +
    if(is.factor(tmp[, 2])){ # Add spacing
      scale_y_discrete(breaks = levels(tmp[, 2])[c(TRUE, rep(FALSE, 9))])
    } else {
      NULL
    }
  return(p)
}
plot_list <- lapply(use, plot_mat)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})
tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/contact_matrices.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
while(!is.null(dev.list())) dev.off()

## Plot population
plot_pop <- function(code, ...){
  get_data(code)
  ggplot(data = data.frame(age = 1 : length(popSize), popSize), 
         mapping = aes(x = age, 
                       y = popSize)) +
    geom_col(width = 1) +
    labs(y = "Population", x = "Age", title = code) +
    coord_flip() 
}
plot_list <- lapply(use, plot_pop)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/population.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
while(!is.null(dev.list())) dev.off()
