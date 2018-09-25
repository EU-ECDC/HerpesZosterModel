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

p1 <- ggplot(data = sizes %>% filter(geo == "DE"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p2 <- ggplot(data = sizes %>% filter(geo == "IE"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p3 <- ggplot(data = sizes %>% filter(geo == "UK"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p4 <- ggplot(data = sizes %>% filter(geo == "IT"),
             mapping = aes(x = age, y = values)) +
  geom_col() +
  facet_grid(n_person ~ geo) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
# TODO: Reorder factors
grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(1, 2), c(3, 4)))
# TODO: Simulate and add simulations to plots
