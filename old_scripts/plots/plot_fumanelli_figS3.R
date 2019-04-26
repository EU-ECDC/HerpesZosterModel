# Load required packages
library(eurostat) # For downloading data
library(plyr)
library(dplyr)
library(ggplot2) # For plotting
library(gridExtra)

# Load data
cens_01nhsize <- get_eurostat(id = "cens_01nhsize")
## Relabel GE7 as 7
cens_01nhsize <- transform(cens_01nhsize,
                           n_person = revalue(n_person, c("GE7" = "7")))

# Recreate Figure S3 from Fumanelli et al 2012
theme_set(theme_bw())
sizes <- cens_01nhsize %>%
  filter(!(n_person %in% c("TOTAL", "UNK")), !(age %in% c("TOTAL", "UNK")),
         sex == "T", citizen == "TOTAL",
         geo %in% c("DE", "UK", "IE", "IT"))

# Reorder age factor
sizes$age <- factor(sizes$age,
                    levels = c("Y_LT5", "Y5-9", "Y10-14", "Y15-19", "Y20-24",
                               "Y25-29", "Y30-34", "Y35-39", "Y40-44", "Y45-49",
                               "Y50-54", "Y55-59", "Y60-64", "Y65-69", "Y70-74",
                               "Y75-79", "Y80-84", "Y85-89", "Y90-94", "Y95-99",
                               "Y_GE100", "TOTAL", "UNK"))

p1 <- ggplot(data = sizes %>% filter(geo == "DE"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- ggplot(data = sizes %>% filter(geo == "IE"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
p3 <- ggplot(data = sizes %>% filter(geo == "UK"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
p4 <- ggplot(data = sizes %>% filter(geo == "IT"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
# TODO: Simulate and add simulations to plots

# Save output
tiff(filename = "../figures/hh_distribution.tif",
     width = 800, height = 600)
grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(1, 2), c(3, 4)))
dev.off()
