# Required packages
library(dplyr)
library(readxl)
library(magrittr)
library(data.table)
library(socialmixr)
library(eurostat)
library(akima)
library(ggplot2)
library(gridExtra)

# Load data
## ESEN2 Seroprevalence
data <- read.table("S:/HelenJohnson/Herpes Zoster/Force of infection/esen2vzv.csv",
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)

data <- data %>%
  mutate(indic = ifelse(STDRES != "NEG", 1, 0))
# Make age numeric variable
data[which(data$AGE == "60+"), ]$AGE <- 60
data$AGE <- sapply(strsplit(data$AGE, split = "-"),
                   function(x) mean(as.numeric(x)))

# Replace age values
age_rpl <- data.frame(AGE = c(0, seq(from = 1, to = 20, by = 1)),
                      to = c((1 + 0.5)/2, seq(from = 1.5, to = 20.5, by = 1)))
setDT(data)
setDT(age_rpl)
data <- data[age_rpl, on = c("AGE"), AGE := to]

# Current options based on availability of data
opts <- cbind(c("Belgium", "Finland", "Germany", "Italy", "Luxembourg", 
                "Netherlands", "United Kingdom", "Serbia", "Slovenia"),
              c("Belgium", "Finland", "Germany", "Italy", "Luxembourg", 
                "Netherlands", "UK", "Serbia", "Slovenia"),
              c("BE", "FI", "DE", "IT", "LU", "NL", "UK", "RS", "SI"),
              c("BE", "FI", "DE", "IT", "LU", "NL", "UK", "RS", "SI"))

# Get and assign the data needed -----------------------------------------------
get_data <- function(i){
  if(i > dim(opts)[1]){stop("Invalid value")}
  # Remove objects
  remove(sero, envir = .GlobalEnv)
  remove(contact_w, envir = .GlobalEnv)
  remove(demfit, envir = .GlobalEnv)
  remove(PS, envir = .GlobalEnv)
  
  # Population data
  pop <- get_eurostat(id = "demo_pjan") %>% # Population size
    filter(geo == opts[i, 3]) %>%
    filter(sex == "T") %>%
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
    filter(time == "2003-01-01")
  # Mortality data
  mort <- get_eurostat(id = "demo_magec") %>% # Number of deaths
    filter(geo == opts[i, 4]) %>%
    filter(sex == "T") %>%
    filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
    filter(time == "2003-01-01")
  
  # Re-order population data by age
  popAge <- droplevels(pop$age)
  levels(popAge)[length(levels(popAge))] <- "0.5"
  popAge <- as.numeric(gsub("[A-z]", "", popAge))
  
  PS <- pop %>% mutate(popAge = popAge) %>%
    arrange(popAge) %>% 
    mutate(age = factor(age, age)) %>%
    select(values)
  PS <- unlist(PS)
  # Re-order mortality data by age
  mortAge <- droplevels(mort$age)
  levels(mortAge)[length(levels(mortAge))] <- "0.5"
  mortAge <- as.numeric(gsub("[A-z]", "", mortAge))
  ND <- mort %>% mutate(mortAge = mortAge) %>%
    arrange(mortAge) %>% 
    mutate(age = factor(age, age)) %>%
    select(values)
  ND <- unlist(ND)
  ND <- ND[c(1 : length(PS))]
  # Fit mortality model
  demfit <- mgcv::gam(ND ~ s(sort(mortAge)[c(1 : length(PS))]), offset = log(PS),
                      family = "poisson", link = "log")
  
  if(i < 8){
    cont <- contact_matrix(polymod,
                           countries = opts[i, 1],
                           filter = ("phys_contact" > 3))$matrix
    sero <- data %>% # Seroprevalence (ESEN2)
      filter(COUNTRY == opts[i, 2])
    
    # Weigh contact matrix by population sizes to obtain c(i, j) from m(i, j)
    pop_weights <- data.frame(age = popAge, size = PS)[1 : dim(cont)[1], ]
    weigh <- function(x){
      cont[x, ] / pop_weights[x, 2]
    }
    tmp <- lapply(1 : dim(cont)[1], weigh)
    contact_w <- do.call(rbind, tmp)
    dimnames(contact_w) <- dimnames(cont)
  }
  if(i == 8){
    # Create Serbia data based on Table 1 from Medić et al 2018
    # Epidemiol Infect 146, 1593–1601
    # https://doi.org/10.1017/S0950268818001619
    zeros <- c(25, 235, 135, 64, 42, 12, 12, 7, 8, 6, 2, 1)
    ones <- c(75, 165, 378, 448, 533, 188, 188, 193, 192, 194, 198, 269)
    indic <- c(rep(0, sum(zeros)), rep(1, sum(ones)))
    age_vals <- c("0", "1–4", "5–9", "10–14", "15–19", "20–24", "25–29",
                  "30–34", "35–39", "40–49", "50–59", "60+")
    age_vals <- c(0, (1 + 4) / 2, (5 + 9) / 2, (10 + 14) / 2, (15 + 19) / 2, 
                  (20 + 24) / 2, (25 + 29) / 2, (30 + 34) / 2, (35 + 39) / 2,
                  (40 + 49) / 2, (50 + 59) / 2, 60)
    serb <- data.frame(COUNTRY = rep("Serbia", length(indic)),
                       AGE = c(rep(age_vals, zeros),
                               rep(age_vals, ones)),
                       indic)
    ## Replace age values
    age_rpl <- data.frame(AGE = c(0, seq(from = 1, to = 20, by = 1)),
                          to = c((1 + 0.5)/2, seq(from = 1.5, to = 20.5, by = 1)))
    setDT(serb)
    setDT(age_rpl)
    serb <- serb[age_rpl, on = c("AGE"), AGE := to]
    sero <- serb
    
    load("S:/HelenJohnson/Herpes Zoster/Data/prem.Rda")
    weigh <- function(x, dat, mat){
      AGE <- droplevels(dat$age)
      levels(AGE)[length(levels(AGE))] <- "0.5"
      AGE <- gsub("[A-z]","", AGE)
      AGE <- as.numeric(AGE)
      AGE <- sort(AGE, index.return = TRUE)$`x`
      dat <- dat[sort(AGE, index.return = TRUE)$ix, ]$values
      return(mat[x, ] / dat[x])
    }
    # Use same weighting when combining as Fumanelli et al.
    mat <- 0.3 * prem$RS$home + 0.18 * prem$RS$school +
      0.19 * prem$RS$work + 0.33 * prem$RS$other
    tmp <- lapply(1 : dim(mat)[1],
                  function(x){weigh(x, 
                                    dat = get_eurostat(id = "demo_pjan") %>% # Population size
                                      filter(geo == opts[i, 3]) %>%
                                      filter(sex == "T") %>%
                                      filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                                      filter(time == "2003-01-01"),
                                    mat = mat)})
    RS <- as.matrix(do.call(rbind, tmp))
    dimnames(RS)[[1]] <- dimnames(RS)[[2]]
    #contact_w <- RS
    
    # Interpolate contact matrix
    tmp <- RS
    tmp <- melt(tmp)
    tmp %>% mutate_if(is.factor, as.character) -> tmp
    tmp[which(tmp$Var1 == "75+"), ]$Var1 <- rep(75, length(tmp[which(tmp$Var1 == "75+"), ]$Var1))
    tmp[which(tmp$Var2 == "75+"), ]$Var2 <- rep(75, length(tmp[which(tmp$Var2 == "75+"), ]$Var2))
    tmp[which(tmp$Var1 == "00-04"), ]$Var1 <- rep(0, length(tmp[which(tmp$Var1 == "00-04"), ]$Var1))
    tmp[which(tmp$Var2 == "00-04"), ]$Var2 <- rep(0, length(tmp[which(tmp$Var2 == "00-04"), ]$Var2))
    tmp$Var1 <- sapply(strsplit(tmp$Var1, split = "-"),
                       function(x) mean(as.numeric(x)))
    tmp$Var2 <- sapply(strsplit(tmp$Var2, split = "-"),
                       function(x) mean(as.numeric(x)))
    
    contact_w <- interp(tmp$Var1, tmp$Var2, tmp$value,
                        nx = 75, ny = 75)$z
    colnames(contact_w) <- rownames(contact_w) <- seq(0, 74, 1)
  }
  if(i == 9){
    # Create Slovenia data based on Figure 1 from Soca et al 2010
    # BMC Public Health 10:360
    # https://doi.org/10.1186/1471-2458-10-360
    tmp <- read_excel("S:/HelenJohnson/Herpes Zoster/Force of infection/old/slovenia_seroprev.xlsx")
    # Convert columns to numeric
    cols <- c(2, 3, 4)
    tmp[, cols] %<>% lapply(function(x) as.numeric(x))
    
    tmp <- na.omit(tmp) # Remove excess rows added due to additional calculations
    # further down in the sheet
    
    # Calculate how many zeros and ones based on seropositivy rate
    tmp <- tmp %>% 
      mutate(zeros = round(n * mid / 100),
             ones = n - zeros)
    
    # Create sero data set
    age_vals <- c(rep(tmp$age, tmp$zeros), rep(tmp$age, tmp$ones))
    indic <- c(rep(rep(0, length(tmp$age)), tmp$zeros), rep(rep(1, length(tmp$age)), tmp$ones))
    slov <- data.frame(COUNTRY = rep("Slovenia", length(indic)),
                       AGE = age_vals,
                       indic)
    sero <- slov
    
    # Contact matrix
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
                                    dat = get_eurostat(id = "demo_pjan") %>% # Population size
                                      filter(geo == opts[i, 3]) %>%
                                      filter(sex == "T") %>%
                                      filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
                                      filter(time == "2003-01-01"),
                                    mat = mat)})
    SI <- as.matrix(do.call(rbind, tmp))
    dimnames(SI)[[1]] <- dimnames(SI)[[2]]
    contact_w <- SI
  }
  
  contact_w[is.na(contact_w)] <- 0 # Replace missings with zeros
  sero <- sero[!is.na(sero$indic) & !is.na(sero$AGE), ]
  sero <- sero[order(sero$AGE), ]
  
  if(dim(pop)[1] == 0)
    warning("Population data for year 2003 not available from Eurostat database demo_pjan")
  if(dim(mort)[1] == 0)
    warning("Mortality data for year 2003 not available from Eurostat database demo_magec")
  
  # Save objects
  assign("sero", sero, envir = .GlobalEnv)
  assign("contact_w", contact_w, envir = .GlobalEnv)
  assign("demfit", demfit, envir = .GlobalEnv)
  assign("PS", PS, envir = .GlobalEnv)
}

# Plots of data ----------------------------------------------------------------

# Load data
countries_of_interest <- opts[, 3]
pop_data <- get_eurostat(id = "demo_pjan")
mort_data <- get_eurostat(id = "demo_magec")

use <- countries_of_interest[countries_of_interest %in% intersect(pop_data$geo, mort_data$geo)]

## Plot mortality
plot_mort <- function(i, ...){
  get_data(i)
  vals <- plot(demfit)
  ggplot(mapping = aes(x = x, y = fit), 
         data = as.data.frame(vals[[1]][c("x", "se", "fit")])) + 
    labs(x = "age", y = "f(age)", title = opts[i, 3]) +
    geom_line() +
    geom_ribbon(aes(ymin = fit - se,
                    ymax = fit + se),
                alpha = 0.2)  + 
    theme(plot.title = element_text(hjust = 0.5))
}

theme_set(theme_classic())
plot_list <- lapply(1 : length(use), plot_mort)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/mortality.tif",
     width = 800, height = 600)
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
while(!is.null(dev.list())) dev.off()

## Plot serological data
plot_sero <- function(i, ...){
  get_data(i)
  subset <- (sero$AGE > 0.5) & (sero$AGE < 80) &
    (!is.na(sero$AGE)) & !is.na(sero$indic)
  sero <- sero[subset, ]
  y <- sero$indic[order(sero$AGE)]
  a <- sero$AGE[order(sero$AGE)]

  grid <- sort(unique(round(a)))
  neg <- table(y, round(a))[1, ]
  pos <- table(y, round(a))[2, ]
  tot <- neg + pos
  
  data <- as.data.frame(cbind(grid, neg, pos, tot))
  
  ggplot(data = data, mapping = aes(x = grid,
                                    y = pos/tot, size = 0.02 * tot)) +
    geom_point(pch = 1) + 
    labs(x = "age", y = "sero-prevalence", title = opts[i, 3])) + 
    xlim(0, 72) + ylim(- 0.1, 1) + theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_list <- lapply(1 : length(use), plot_sero)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/serology.tif",
     width = 800, height = 600)
grid.arrange(BE, FI, DE, IT, LU, NL, UK, RS, SI)
while(!is.null(dev.list())) dev.off()