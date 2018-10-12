# Load required packages
library(dplyr) # For data
library(tidyr)
library(ggplot2) # For plotting
library(gridExtra)
library(viridis) # For colours
library(scales)
#TODO: use library(gganimate)

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

# Plot smoothed projections
theme_set(theme_bw())
ggplot(data = data, 
       mapping = aes(x = Time, colour = Variant,
                     y = Population, linetype = Sex)) + 
geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  facet_wrap(. ~ Location, scales = "free_y") +
  scale_colour_viridis(option = "A", discrete = TRUE)

# Load UN WPP data by age
file <- "https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2017_PopulationByAgeSex_Medium.csv"
med <- read.csv(file, header = TRUE)
file <- "https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2017_PopulationByAgeSex_OtherVariants.csv"
others <- read.csv(file, header = TRUE)
names(others)[1] <- "LocID"
data <- full_join(med, others)

# Select EU/EEA countries and create sex variable
data <- data %>%
  filter(LocID %in% c(40, 56, 100, 191, 196, 203, 208, 233, 246, 250, 276, 300,
                      348, 352, 372, 380, 428, 438, 440, 442, 470, 528, 578, 616,
                      620, 642, 703, 705, 724, 752, 826)) %>%
  gather(Sex, Population, PopMale:PopTotal, factor_key = TRUE)
# Turn age group into factor
data$AgeGrp <- factor(data$AgeGrp,
                      levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                                 "30-34", "35-39", "40-44", "45-49", "50-54",
                                 "55-59", "60-64", "65-69", "70-74", "75-79",
                                 "80+", "80-84", "85-89", "90-94", "95-99", "100+"))

# Plot pyramids over time (three snapshots)
p1 <- ggplot(data = data %>% filter(Time == 2015) %>% filter(Sex == "PopTotal"), 
             mapping = aes(x = AgeGrp, fill = Sex, 
                           y = ifelse(test = Sex == "PopMale", 
                                      yes = -Population, no = Population),
                           colour = Variant)) +
  geom_bar(stat = "identity") +
  labs(y = "population") +
  coord_flip() +
  facet_wrap(. ~ Location, scales = "free_x") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(option = "A", discrete = TRUE) +
  scale_fill_grey(start = 0.6, end = 0.4) +
  guides(colour = guide_legend(override.aes = list(fill = NA)))
p2 <- ggplot(data = data %>% filter(Time == 2050) %>% filter(Sex == "PopTotal"), 
             mapping = aes(x = AgeGrp, fill = Sex, 
                           y = ifelse(test = Sex == "PopMale", 
                                      yes = -Population, no = Population),
                           colour = Variant)) +
  geom_bar(stat = "identity") +
  labs(y = "population") +
  coord_flip() +
  facet_wrap(. ~ Location, scales = "free_x") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(option = "A", discrete = TRUE) +
  scale_fill_grey(start = 0.6, end = 0.4) +
  guides(colour = guide_legend(override.aes = list(fill = NA)))
p3 <- ggplot(data = data %>% filter(Time == 2100) %>% filter(Sex == "PopTotal"), 
             mapping = aes(x = AgeGrp, fill = Sex, 
                           y = ifelse(test = Sex == "PopMale", 
                                      yes = -Population, no = Population),
                           colour = Variant)) +
  geom_bar(stat = "identity") +
  labs(y = "population") +
  coord_flip() +
  facet_wrap(. ~ Location, scales = "free_x") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(option = "A", discrete = TRUE) +
  scale_fill_grey(start = 0.6, end = 0.4) +
  guides(colour = guide_legend(override.aes = list(fill = NA)))
grid.arrange(p1, p2, p3, ncol = 3)

# Plot changes in age group 0-4 over time
ggplot(data = left_join(data %>%
                          filter(Sex == "PopTotal") %>%
                          filter(AgeGrp == "0-4") %>%
                          group_by(Location, Variant),
                        data %>%
                          filter(Sex == "PopTotal") %>%
                          filter(AgeGrp == "0-4", Time == 2015) %>%
                          group_by(Location, Variant) %>%
                          mutate(Subtr = Population)) %>%
         select(Time, Location, Variant, Population, Subtr) %>%
         mutate(Time =  factor(Time, levels = c(2015,
                                                seq(1950, 2015 - 1, 1),
                                                seq(2015 + 1, 2100, 1)))) %>%
         group_by(Location, Variant) %>%
         arrange(Time) %>%
         fill(Subtr, Subtr) %>%
         mutate(Time = factor(Time, levels = c(seq(1950, 2100, 1)))) %>%
         mutate (PopDiff = Population - Subtr),
       mapping = aes(x = Time, y = PopDiff,
                     colour = Location, group = Location)) +
  geom_line(aes(y = PopDiff)) +
  labs(title = "Changes from 2015 values", subtitle = "Population aged 0-4") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ Variant) +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  scale_x_discrete(breaks = c("1950", "2000", "2015", "2050", "2100"),
                   labels = c("1950", "2000", "2015", "2050", "2100"))

# Load indicators
file <- "https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2017_Period_Indicators_Medium.csv"
med <- read.csv(file, header = TRUE, na.strings = "NULL")
file <- "https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2017_Period_Indicators_OtherVariants.csv"
others <- read.csv(file, header = TRUE, na.strings = "NULL")
names(others)[1] <- "LocID"
others$TFR <- as.numeric(levels(others$TFR)[others$TFR])
others$NRR <- as.numeric(levels(others$NRR)[others$NRR])
others$CBR <- as.numeric(levels(others$CBR)[others$CBR])
others$Births <- as.numeric(levels(others$Births)[others$Births])
others$CDR <- as.numeric(levels(others$CDR)[others$CDR])
others$Deaths <- as.numeric(levels(others$Deaths)[others$Deaths])
others$DeathsMale <- as.numeric(levels(others$DeathsMale)[others$DeathsMale])
others$DeathsFemale <- as.numeric(levels(others$DeathsFemale)[others$DeathsFemale])
others$NatIncr <- as.numeric(levels(others$NatIncr)[others$NatIncr])
data <- full_join(med, others)

# Select EU/EEA countries
data <- data %>%
  filter(LocID %in% c(40, 56, 100, 191, 196, 203, 208, 233, 246, 250, 276, 300,
                      348, 352, 372, 380, 428, 438, 440, 442, 470, 528, 578, 616,
                      620, 642, 703, 705, 724, 752, 826))

# Plot indicators
ggplot(data = data,
       mapping = aes(x = Time,
                     y = NetMigrations,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Net migration") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE)

ggplot(data = data, 
       mapping = aes(x = Time,
                     y = IMR,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Infant mortality rate") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE)

ggplot(data = data, 
       mapping = aes(x = Time,
                     y = GrowthRate,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Population growth rate") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE)
