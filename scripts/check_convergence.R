# Data
data <- read.table("S:\\HelenJohnson\\Herpes Zoster\\Force of infection\\esen2vzv.csv",
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
library(dplyr)
data <- data %>%
  mutate(indic = ifelse(STDRES != "NEG", 1, 0))
data$AGE <- sapply(strsplit(data$AGE, split = "-"),
                   function(x) mean(as.numeric(x)))

library(socialmixr)
cont <- contact_matrix(polymod,
                          countries = "Finland",
                          filter = ("phys_contact" > 3))$matrix
sero <- data %>% 
  filter(COUNTRY == "Finland")
library(eurostat)
pop <- get_eurostat(id = "demo_pjan") %>% # Population size
  filter(geo == "FI") %>%
  filter(sex == "T") %>%
  filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
  filter(time == "2003-01-01")
mort <- get_eurostat(id = "demo_magec") %>% # Number of deaths
  filter(geo == "FI") %>%
  filter(sex == "T") %>%
  filter(!(age %in% c("TOTAL", "UNK", "Y_OPEN"))) %>%
  filter(time == "2006-01-01")

# Contact function with added print.level option
if(!packageVersion("FOIest") == "0.2.1"){
  install.packages("S:\\HelenJohnson\\Herpes Zoster\\Force of infection\\FOIest_0.2.1.tar.gz", repos = NULL, type = "source")
}
library(FOIest)

# Inner workings of example_calc
PS <- pop
AGE <- droplevels(PS$age)
levels(AGE)[length(levels(AGE))] <- "0.5"
AGE <- gsub("[A-z]", "", AGE)
AGE <- as.numeric(AGE)
AGE <- sort(AGE, index.return = TRUE)$`x`
PS <- PS[sort(AGE, index.return = TRUE)$ix, ]$values
ND <- mort
ND <- ND[sort(AGE, index.return = TRUE)$ix, ]$values
demfit <- mgcv::gam(ND ~ s(AGE), offset = log(PS),
                    family = "poisson", link = "log")
pop_weights <- data.frame(age = AGE, size = PS)[1 : dim(cont)[1], ]
weigh <- function(x){
  cont[x, ] / pop_weights[x, 2]
}
tmp <- lapply(1 : dim(cont)[1], weigh)
contact_w <- do.call(rbind, tmp)
dimnames(contact_w) <- dimnames(cont)
contact_w[is.na(contact_w)] <- 0 # Replace missings with zeros
sero <- sero[!is.na(sero$indic) & !is.na(sero$AGE), ]
sero <- sero[order(sero$AGE), ]
res <- FOIest::contact(a = sero$AGE, y = sero$indic, rij = contact_w,
               muy = predict(demfit, type = "response"),
               N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 85,
               startpar = 5e-2, print = 2)

res <- FOIest::contact(a = sero$AGE, y = sero$indic, rij = contact_w,
                       muy = predict(demfit, type = "response"),
                       N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 85,
                       startpar = c(5e-2, 0), print = 2, prop = "loglin")

# Check UK with two params
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/scripts/exmpl_calc.R",
                 ssl.verifypeer = FALSE)
invisible(suppressMessages(eval(parse(text = script))))

run_calc(7, D = 6 / 365, A = 0.5, plots = FALSE, startpar = c(5e-2, 0),
         print = 2, prop = "loglin")
