# Load required packages
library(eurostat) # For downloading data
library(dplyr)
library(ggplot2) # For plotting
library(viridis)

# Load data
# Eurostat data "Children in formal childcare or education by age group and
# duration - % over the population of each age group - EU-SILC surveyp"

ilc_caindformal <- get_eurostat(id = "ilc_caindformal")
ilc_caindformal$age <- factor(ilc_caindformal$age,
                              levels = c("Y_LT3", "Y3-CSA", "CSA-Y12"),
                              labels = c("Less than 3 years",
                                         "From 3 yrs to MCSA",
                                         "From MCSA to 12 years"))
ilc_caindformal$duration <- factor(ilc_caindformal$duration,
                              levels = c("H0", "H1-29", "H_GE30"),
                              labels = c("Zero hours",
                                         "From 1 to 29 hours",
                                         "30 hours or over"))

theme_set(theme_bw())
# European level values
(eur <- ggplot(data = ilc_caindformal %>% filter(geo == "EU28"),
       mapping = aes(x = time, y = values,
                     group = interaction(time, age), colour = age)) +
  geom_boxplot() +
  facet_grid(. ~ duration) + 
  scale_colour_viridis(discrete = TRUE,  begin = 0, end = 0.7) +
  labs(title = "EU28", subtitle = "MCSA denotes minimum compulsory school age",
       y = "Percent") +
  theme(plot.title = element_text(hjust = 0.5)))

# Country level values
(cou <- ggplot(data = ilc_caindformal %>%
         filter(geo %in% c("AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "ES",
                           "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT",
                           "LI", "LT", "LU", "LV", "MT", "NL", "NO", "PL", "PT",
                           "RO", "SE", "SI", "SK")),
       mapping = aes(x = time, y = values,
                     group = geo, colour = geo)) +
  geom_line(stat = "summary", fun.y = "mean")  +
  facet_grid(age ~ duration) + 
  scale_colour_viridis(discrete = TRUE, option = "B") +
  labs(title = "EU28", subtitle = "MCSA denotes minimum compulsory school age",
       y = "Percent") +
  theme(plot.title = element_text(hjust = 0.5)))

# Save outputs
tiff(filename = "../figures/child_care_eu28_overall.tif",
     width = 1200, height = 800)
eur
dev.off()

tiff(filename = "../figures/child_care_eu28_per_country.tif",
     width = 1200, height = 800)
cou
dev.off()
