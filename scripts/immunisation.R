library(readxl)
library(purrr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(viridis)

# Load WHO/UNICEF estimates of national immunisation coverage
file <- "http://www.who.int/entity/immunization/monitoring_surveillance/data/coverage_estimates_series.xls"

sheets <- excel_sheets(file)
sheets <- sheets[-1]

data <- map_dfr(sheets, ~read_excel(path = file, sheet = .x)) %>%
  filter(ISO_code %in% c("AUT", "BEL", "BGR", "CYP", "CZE", "DEU", "DNK", "ESP",
                         "EST", "FIN", "FRA", "GBR", "GRC", "HRV", "HUN", "IRL",
                         "ISL", "ITA", "LIE", "LTU", "LUX", "LVA", "MLT", "NLD",
                         "NOR", "POL", "PRT", "ROU", "SVK", "SVN", "SWE"))

data$X__1 <- NA
data <- gather(data, year, value, -Region, -ISO_code, -Cname, -Vaccine)
data <- data %>% filter(year != "X__1")

# Plot mean immunisation coverage of 13 vaccines
## including mean of all countries
theme_set(theme_bw())
(p <- ggplot(data = data %>%
         group_by(Cname, year) %>%
         summarise(immunisation = mean(value, na.rm = TRUE)),
       mapping = aes(x = year, y = immunisation,
                     group = interaction(Cname, year), colour = Cname)) +
  stat_summary(aes(group = 1), fun.y = mean, geom = "line", lwd = 2) +
  geom_line(aes(group = Cname), alpha = 0.7) +
  labs(title = "Mean immunisation",
       subtitle = "Of BCG, DTP1, DTP3, HepB_BD, HepB3, Hib3, IPV1, MCV1, MCV2, PCV3, Pol3, RCV1, and RotaC, as well as average of these",
       y = "Percent") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_viridis(discrete = TRUE, option = "C"))

# Save output
tiff(filename = "../figures/mean_immunisation.tif",
     width = 1200, height = 800)
p
dev.off()
