# Load required packages
library(dplyr) # For data
library(tidyr)
library(ggplot2) # For plotting
library(gridExtra)
library(viridis) # For colours
library(scales)
library(gganimate)

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
(p <- ggplot(data = data, 
            mapping = aes(x = Time, colour = Variant,
                          y = Population, linetype = Sex)) + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  facet_wrap(. ~ Location, scales = "free_y") +
  scale_colour_viridis(option = "A", discrete = TRUE) +
  geom_vline(mapping = aes(xintercept = 2015)) +
  labs(y = "Population in thousands"))

# Save output
tiff(filename = "../figures/pop_pred_curves.tif",
     width = 800, height = 500)
p
dev.off()

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

# Plot pyramids over time
# One country selected as example
use_data <- data %>%
  filter(LocID == 826) %>%
  filter(AgeGrp %in% c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                       "30-34", "35-39", "40-44", "45-49", "50-54",
                       "55-59", "60-64", "65-69", "70-74", "75-79")) %>%
  # droplevels didn't seem to work
  filter(Sex != "PopTotal") %>%
  filter(VarID == 2) # Medium variant

# Plot pyramids
p <- ggplot(data = use_data,
            mapping = aes(x = AgeGrp, fill = Sex)) +
  geom_bar(data = subset(use_data, Sex == "PopFemale"), 
           mapping = aes(y = Population * (- 1)),
           stat = "identity") + 
  geom_bar(data = subset(use_data, Sex == "PopMale"),
           mapping = aes(y = Population),
           stat = "identity") + 
  labs(y = "Population") +
  scale_y_continuous(breaks = seq(- 35e4, 35e4, by = 10e4),
                     labels = abs(seq(- 35e4, 35e4, by = 10e4))) +
  coord_flip() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_fill_viridis(option = "D", discrete = TRUE) +
  transition_time(Time) +
  labs(title = "{frame_time} (UN WPP)")

# It goes too fast by default slow it down
animate(plot = p, duration = 150)

# Save animation
anim_save(filename = "temp.gif",
          animation = last_animation(),
          path = NULL)

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
(p1 <- ggplot(data = data,
       mapping = aes(x = Time,
                     y = NetMigrations,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Net migration") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE))

(p2 <- ggplot(data = data, 
       mapping = aes(x = Time,
                     y = IMR,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Infant mortality rate") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE))

(p3 <- ggplot(data = data, 
       mapping = aes(x = Time,
                     y = GrowthRate,
                     group = Variant, colour = Variant)) +
  geom_line() +
  facet_wrap(. ~ Location, scales = "free_y") +
  labs(title = "Population growth rate") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 1)) +
  scale_colour_viridis(discrete = TRUE))

# Save outputs
tiff(filename = "../figures/changes_migration.tif",
     width = 1800, height = 1000)
p1
dev.off()

tiff(filename = "../figures/changes_inft_mortality.tif",
     width = 1800, height = 1000)
p2
dev.off()

tiff(filename = "../figures/changes_pop_growth.tif",
     width = 1800, height = 1000)
p3
dev.off()
