rm(list = ls()) # Clear workspace

# Calculating FoI and prevalence

# Source with specific line options, from https://stackoverflow.com/a/12215731
source2 <- function(file, start, end, ...){
  file.lines <- scan(file, what = character(), skip = start - 1,
                     nlines = end - start + 1, sep = "\n")
  file.lines.collapsed <- paste(file.lines, collapse = "\n")
  source(textConnection(file.lines.collapsed), ...)
}

# Load packages
# Internal FOI estimation package
#install.packages("S:\\HelenJohnson\\Herpes Zoster\\Force of infection\\FOIest_0.2.0.tar.gz", repos = NULL, type = "source")
library(FOIest)
library(magrittr)
library(tidylog)
library(dplyr)
library(ggplot2)
library(plotly)
library(socialmixr)
library(gridExtra)
library(eurostat)
source("https://raw.githubusercontent.com/EU-ECDC/EcdcSupportSnippets/master/R/eurostat_wrappers.R")

# Load serology data
data <- read.table("S:\\HelenJohnson\\Herpes Zoster\\Force of infection\\esen2vzv.csv",
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
data <- data %>% # Create indicator variable for having experienced infection
  mutate(indic = ifelse(STDRES != "NEG", 1, 0))
# Note - currently putting equivocal in negative set
data$AGE <- sapply(strsplit(data$AGE, split = "-"),
                   function(x) mean(as.numeric(x))) # Replace ranges with midpoint

# Countries we can do the analysis on
opts <- unique(data$COUNTRY)
opts[which(opts == "UK")] <- "United Kingdom"
opts <- intersect(opts, unique(polymod$participants$country))

# Country example --------------------------------------------------------------
EEA <- c("AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL",
         "HU", "IS", "IE", "IT", "LV", "LI", "LT", "LU", "MT", "NL", "NO", "PL", 
         "PT", "RO", "SK", "SI", "ES", "SC", "SE", "UK") 
countryList <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", 
                 "Czech Republic", "Denmark", "Estonia", "Finland", "France",
                 "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
                 "Latvia", "Liechtenstein", "Lithuania","Luxembourg", "Malta",
                 "Netherlands", "Norway", "Poland", "Portugal", "Romania",
                 "Slovakia", "Slovenia", "Spain", "Scotland", "Sweden", 
                 "United Kingdom" )
countries <- tibble(code = EEA, name = countryList)

## Generic example
## We consider Finland
thisCountry <- opts[2]
thisCode <- countries %>% filter(name == thisCountry) %>% select(code)
example_calc <- function(cont = contact_matrix(polymod,
                                               countries = thisCountry,
                                               filter = ("phys_contact" > 3))$matrix,
                         # Filtering for contacts > 15 mins
                         sero = data %>% 
                           filter(COUNTRY == thisCountry),
                         pop = get_eurostat(id = "demo_pjan") %>% # Population size
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2003-01-01"),
                         mort = get_eurostat(id = "demo_magec") %>% # Number of deaths
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2006-01-01"),
                         plot_inputs = TRUE, ...){
  if(isTRUE(plot_inputs)){
    source2("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/scripts/plot_polymod_matrix", 1, 15)
    source2("S:\\HelenJohnson\\Herpes Zoster\\Force of infection\\sero_europe.R", 13, 33)
    ## Plot
    theme_set(theme_classic())
    grid.arrange(plot_mat(na.omit(cont), option = "E", limits = c(0, 1)),
                 sero_plot(sero))
  }
  PS <- pop
  # Turn age into numerical variable
  AGE <- droplevels(PS$age)
  levels(AGE)[length(levels(AGE))] <- "0.5"
  # Remove Y-prefix
  AGE <- gsub("[A-z]", "", AGE)
  AGE <- as.numeric(AGE)
  AGE <- sort(AGE, index.return = TRUE)$`x`
  PS <- PS[sort(AGE, index.return = TRUE)$ix, ]$values
  ND <- mort
  ND <- ND[sort(AGE, index.return = TRUE)$ix, ]$values
  # Fit mortality data
  demfit <- mgcv::gam(ND ~ s(AGE), offset = log(PS),
                      family = "poisson", link = "log")
  ## Calculate contact rates from contact matrices using population weights
  pop_weights <- data.frame(age = AGE, size = PS)[1 : dim(cont)[1], ]
  
  weigh <- function(x){
    cont[x, ] / pop_weights[x, 2]
  }
  
  tmp <- lapply(1 : dim(cont)[1], weigh)
  contact_w <- do.call(rbind, tmp)
  dimnames(contact_w) <- dimnames(cont)
  contact_w[is.na(contact_w)] <- 0 # Replace missings with zeros
  
  # Remove NAs and order by age (pre-processing)
  sero <- sero[!is.na(sero$indic) & !is.na(sero$AGE), ]
  sero <- sero[order(sero$AGE), ]
  
  # Fit FOI with constant proportionality factor in social contact hypothesis
  res <- FOIest::contact(age = sero$AGE, y = sero$indic, rij = contact_w,
                         muy = predict(demfit, type = "response"),
                         N = sum(PS), ...)
  return(res)
}
vals <- invisible(example_calc(D = 6 / 365, A = 0.5,
                               Lmax = 85, plots = FALSE, startpar = 5e-2))
rbind(c("q", vals$qhat),
      c("R0", vals$R0), 
      c("R", vals$R))
# Plot results
graphics::layout(mat = matrix(cbind(c(1, 2)), ncol = 2, byrow = TRUE))
plot(vals$lambda, type = "l", ylim = c(0, 1), main = "FoI", xlab = "", ylab = "")
plot(vals$pi, type = "l", ylim = c(0, 1), main = "Prev", xlab = "", ylab = "")

# Wrap this in a function
run_calc <- function(country, ...){
  thisCountry <- opts[country]
  thisCode <- countries %>% filter(name == thisCountry) %>% select(code)
  if(!(thisCode == "UK")){
    vals <- example_calc(cont = contact_matrix(polymod,
                                               countries = thisCountry,
                                               filter = ("phys_contact" > 3))$matrix,
                         # Filtering for contacts > 15 mins
                         sero = data %>% 
                           filter(COUNTRY == thisCountry),
                         pop = get_eurostat(id = "demo_pjan") %>% # Population size
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2003-01-01"),
                         mort = get_eurostat(id = "demo_magec") %>% # Number of deaths
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2006-01-01"),
                         plot_inputs = FALSE,
                         Lmax = dim(contact_matrix(polymod,
                                               countries = thisCountry,
                                               filter = ("phys_contact" > 3))$matrix)[1], ...)}
  if(thisCode == "UK"){
    vals <- example_calc(cont = contact_matrix(polymod,
                                               countries = thisCountry,
                                               filter = ("phys_contact" > 3))$matrix,
                         # Filtering for contacts > 15 mins
                         sero = data %>% 
                           filter(COUNTRY == thisCode$code),
                         pop = get_eurostat(id = "demo_pjan") %>% # Population size
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2003-01-01"),
                         mort = get_eurostat(id = "demo_magec") %>% # Number of deaths
                           filter(geo == thisCode$code) %>%
                           filter(sex == "T") %>%
                           filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                           filter(time == "2006-01-01"),
                         plot_inputs = FALSE, 
                         Lmax = dim(contact_matrix(polymod,
                                               countries = thisCountry,
                                               filter = ("phys_contact" > 3))$matrix)[1], ...)}
  
  return(rbind(c("q", vals$qhat),
               c("R0", vals$R0), 
               c("R", vals$R)))
}

# Italy
run_calc(4, D = 6 / 365, A = 0.5, plots = FALSE, startpar = 5e-2)
# UK
run_calc(7, D = 6 / 365, A = 0.5, plots = FALSE, startpar = 5e-2)

# Range of starting values -----------------------------------------------------
runs <- setNames(lapply(1 : length(opts),
                        function(x){run_calc(x, D = 6 / 365, A = 0.5, 
                                             plots = FALSE, startpar = 5e-2)}), 
                 opts) 

# Which ones returned starting parameter value?
vals <- which(sapply(runs, function(x) x[1, 2]) == 5e-2)

# How different are the ones that didn't with other starting values?
suppressMessages(setNames(lapply(which(!(c(1:length(opts)) %in% vals)), 
                                 function(x){run_calc(x, D = 6 / 365, A = 0.5, 
                                                      plots = FALSE, startpar = 5e-2)}),
                          opts[which(!(c(1:length(opts)) %in% vals))]))
suppressMessages(setNames(lapply(which(!(c(1:length(opts)) %in% vals)), 
                                 function(x){run_calc(x, D = 6 / 365, A = 0.5, 
                                                      plots = FALSE, startpar = 1e-2)}),
                          opts[which(!(c(1:length(opts)) %in% vals))]))

# How different are the ones that did with other starting values?
suppressMessages(setNames(lapply(vals, 
                                 function(x){run_calc(x, D = 6 / 365, A = 0.5, 
                                                      plots = FALSE, startpar = 5e-2)}),
                          opts[vals]))
