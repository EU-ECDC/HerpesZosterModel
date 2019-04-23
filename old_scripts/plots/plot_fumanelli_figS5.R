# Load required packages
library(eurostat) # For downloading data
library(dplyr)
library(ggplot2) # For plotting

# Load data
# Eurostat data "Population by sex, age group and country of citizenship"
cens_01nsctz <- get_eurostat(id = "cens_01nsctz")

# Recreate Figure S5 from Fumanelli et al 2012
theme_set(theme_bw())
# Calculate percentages
data <- cens_01nsctz %>%
  filter(!(age %in% c("TOTAL", "UNK")), sex == "T") %>%
  group_by(age, geo) %>%
  select(age, geo, values) %>% 
  summarise(n = sum(values)) %>%
  group_by(geo) %>%
  mutate(perc = (n / sum(n)) * 100) 
# Relevel age factor to be in correct order
data$age <- factor(data$age, levels = c("Y_LT5", "Y5-9", "Y10-14", "Y15-19",
                                        "Y20-24", "Y25-29", "Y30-34", "Y35-39",
                                        "Y40-44", "Y45-49", "Y50-54", "Y55-59",
                                        "Y60-64", "Y65-69", "Y70-74", "Y75-79",
                                        "Y80-84", "Y_GE85", "TOTAL", "UNK"))
# Plot
(p <- ggplot(data = data, mapping = aes(x = age, y = perc, group = 1)) +
  geom_line() +
  facet_wrap(. ~ geo) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Percent"))
# TODO: Simulate and add simulations to plots

# Save output
tiff(filename = "../figures/age_distr_percentage .tif",
     width = 1000, height = 600)
p
dev.off()
