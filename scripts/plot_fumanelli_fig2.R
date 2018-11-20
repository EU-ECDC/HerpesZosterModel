# Load contact matrices
load(file = "../data/fumanelli.Rda")

# Load required packages
library(reshape2) # For melting data for plotting
library(ggplot2) # For plotting
theme_set(theme_bw())
library(viridis) # Colour blind friendly and prints well in black and white
library(gridExtra) # For combining plots

# Recreate Figure 2 from Fumanelli et al 2012 (the UK)
uk_melt <- melt(as.matrix(fumanelli$GBR$household))
hm1 <- ggplot(data = uk_melt, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_viridis(option = "E", limits = c(-2, 3)) +
  labs(x = "Age", y = "Age of contact", title = "Household")
uk_melt <- melt(as.matrix(fumanelli$GBR$school))
hm2 <- ggplot(data = uk_melt, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_viridis(option = "E", limits = c(-2, 3)) +
  labs(x = "Age", y = "Age of contact", title = "School")
uk_melt <- melt(as.matrix(fumanelli$GBR$workplace))
hm3 <- ggplot(data = uk_melt, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_viridis(option = "E", limits = c(-2, 3)) +
  labs(x = "Age", y = "Age of contact", title = "Workplace")
uk_melt <- melt(as.matrix(fumanelli$GBR$general))
hm4 <- ggplot(data = uk_melt, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_viridis(option = "E", limits = c(-2, 3)) +
  labs(x = "Age", y = "Age of contact", title = "General community")
total <- fumanelli$GBR$household + fumanelli$GBR$school + 
  fumanelli$GBR$workplace + fumanelli$GBR$general
uk_melt <- melt(as.matrix(total))
hm5 <- ggplot(data = uk_melt, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_viridis(option = "E", limits = c(-2, 3)) +
  labs(x = "Age", y = "Age of contact", title = "Total")

# Save output
tiff(filename = "../figures/contact_uk.tif",
     width = 1000, height = 600)
grid.arrange(hm1, hm2, hm3, hm4, hm5,
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, NA)))
dev.off()

#TODO: Add plot f
#It would seem that the percentages are created with something like 
sum(total[1:4, 1:4]) / sum(total[, 1:4]) * 100
sum(total[5:9, 5:9]) / sum(total[, 5:9]) * 100
sum(total[10:14, 10:14]) / sum(total[, 10:14]) * 100
sum(total[15:19, 15:19]) / sum(total[, 15:19]) * 100
sum(total[20:24, 20:24]) / sum(total[, 20:24]) * 100
sum(total[25:29, 25:29]) / sum(total[, 25:29]) * 100
sum(total[30:34, 30:34]) / sum(total[, 30:34]) * 100
sum(total[35:39, 35:39]) / sum(total[, 35:39]) * 100
sum(total[40:44, 40:44]) / sum(total[, 40:44]) * 100
sum(total[45:49, 45:49]) / sum(total[, 45:49]) * 100
sum(total[50:54, 50:54]) / sum(total[, 50:54]) * 100
sum(total[55:59, 55:59]) / sum(total[, 55:59]) * 100
sum(total[60:64, 60:64]) / sum(total[, 60:64]) * 100
sum(total[65:69, 65:69]) / sum(total[, 65:69]) * 100
sum(total[70:dim(total)[1], 70:dim(total)[1]]) / 
  sum(total[, dim(total)[1]:70])  * 100