# Load required packages
library(eurostat) # For downloading data
library(plyr)
library(dplyr)
library(ggplot2) # For plotting
library(gridExtra) 

country_codes <- c("AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES",
                   "FI", "FR", "HR", "HU", "IE", "IT", "IS", "LI", "LT", "LU",
                   "LV", "MT", "NL", "NO", "PL", "PT", "RO", "SE", "SI", "SK",
                   "UK")

# Load household data
cens_01nhsize <- get_eurostat(id = "cens_01nhsize")
cens_01nhsize <- as_tibble(transform(cens_01nhsize,
                                     n_person = plyr::revalue(n_person, c("GE7" = "7")),
                                     geo = plyr::revalue(geo, c("CY000" = "CY"))))
detach(package:plyr) # To avoid problems with summarise

hou <- cens_01nhsize %>%
  filter(!(n_person %in% c("TOTAL", "UNK")), !(age %in% c("TOTAL", "UNK")),
         sex == "T", citizen == "TOTAL", geo %in% country_codes) %>%
  group_by(geo) %>%
  select(n_person, age, geo, values)

# Remove values that aren't given in a range
hou_values <- hou %>%
  filter(n_person != "GE7", !(age %in% c("Y_LT5", "Y_GE100")))

## Reorder and convert to numeric
hou_values$n_person <- factor(hou_values$n_person, levels = c("1", "2", "3", "4", "5", "6", "7", "TOTAL", "UNK"))
hou_values$n_person <- as.numeric(as.character(hou_values$n_person))

## Reorder, then use medians, and convert to numeric
hou_values$age <- factor(hou_values$age, levels = c("Y_LT5", "Y5-9", "Y10-14", "Y15-19",
                                                    "Y20-24", "Y25-29", "Y30-34", "Y35-39",
                                                    "Y40-44", "Y45-49", "Y50-54", "Y55-59",
                                                    "Y60-64", "Y65-69", "Y70-74", "Y75-79",
                                                    "Y80-84", "Y85-89", "Y90-94", "Y95-99",
                                                    "Y_GE100", "TOTAL", "UNK"))
levels(hou_values$age) <- c("5", "7", "12", "17", "22", "27", "32", "37", "42", "47",
                            "52", "57", "62", "67", "72", "77", "82", "87", "92", "97",
                            "100", "TOTAL", "UNK")
hou_values$age <- as.numeric(as.character(hou_values$age))

# calculate average size and age, and combine
mean_size <- hou_values %>%
  group_by(geo, n_person) %>%
  summarise(pop = sum(values, na.rm = TRUE)) %>%
  mutate(val = pop * n_person) %>%
  group_by(geo) %>%
  summarise(mean = sum(val) / sum(pop))
  
mean_age <- hou_values %>%
  group_by(geo, age) %>%
  summarise(pop = sum(values, na.rm = TRUE)) %>%
  mutate(val = pop * age) %>%
  group_by(geo) %>%
  summarise(mean = sum(val) / sum(pop))

for_plot <- full_join(mean_size, mean_age, by = "geo")

# Plot points
theme_set(theme_bw())
ggplot(data = for_plot, mapping = aes(x = mean.x,
                                      y = mean.y, colour = geo)) + 
  geom_point() +
  ylim(34, 42) + 
  labs(x = "size", y = "age")
