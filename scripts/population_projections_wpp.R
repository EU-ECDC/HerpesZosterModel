# Load required packages
library(dplyr)
library(ggplot2) # For plotting
library(tidyr)

# Load UN WPP data
file <- "https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2017_TotalPopulationBySex.csv"
data <- read.csv(file, header = TRUE)

# Select EU/EEA countries and create sex variable
data <- data %>%
  filter(LocID %in% c(40, 56, 100, 191, 196, 203, 208, 233, 246, 250, 276, 300,
                    348, 352, 372, 380, 428, 438, 440, 442, 470, 528, 578, 616,
                    620, 642, 703, 705, 724, 752, 826)) %>%
  gather(Sex, Population, PopMale:PopTotal, factor_key = TRUE) %>%
  filter(Sex != "PopTotal")
#data %>% glimpse()

# Plot smoothed projections
theme_set(theme_bw())
ggplot(data = data, 
       mapping = aes(x = Time, colour = Variant,
                     y = Population, linetype = Sex)) + 
geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  facet_wrap(. ~ Location, scales = "free_y")