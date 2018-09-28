# Load required packages
library(eurostat) # For downloading data
library(dplyr)
library(ggplot2) # For plotting

# Load data
# Eurostat data "Population by sex, age group and country of citizenship"
cens_01nsctz <- get_eurostat(id = "cens_01nsctz")
cens_01nsctz$age <- factor(cens_01nsctz$age,
                           levels = c("Y_LT5", "Y5-9", "Y10-14", "Y15-19",
                                      "Y20-24", "Y25-29", "Y30-34", "Y35-39",
                                      "Y40-44", "Y45-49", "Y50-54", "Y55-59",
                                      "Y60-64", "Y65-69", "Y70-74", "Y75-79",
                                      "Y80-84", "Y_GE85", "TOTAL", "UNK"))

#cens_01nsctz %>% glimpse()
cens_01nsctz <- cens_01nsctz %>%
  filter(sex != "T", age != "TOTAL", age != "UNK")

theme_set(theme_bw())
ggplot(data = cens_01nsctz, 
       mapping = aes(x = age, fill = sex, 
                     y = ifelse(test = sex == "M", 
                                yes = -values, no = values))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(cens_01nsctz$values) * c(-1, 1)) +
  labs(y = "population") +
  coord_flip() +
  facet_wrap(. ~ geo, scales = "free_x") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1))