library(socialmixr)
library(reshape2)
library(ggplot2)
library(viridis)
library(gridExtra)

plot_mat <- function(data, ...){
  data <- melt(data)
  ggplot(data = data, aes_string(x = names(data)[1], y = names(data)[2], fill = names(data)[3])) + 
    geom_tile() + 
    scale_fill_viridis_c(...) +
    labs(x = "", y = "") + 
    theme(plot.title = element_text(hjust = 0.5))
}

polymodIT <- contact_matrix(survey = polymod, countries = "IT",
                            n = 2)

plot_mat(na.omit(polymodIT[[1]][1][[1]]$matrix), option = "E", limits = c(0, 2))
plot_mat(polymodIT[[1]][1][[1]]$matrix, option = "A", limits = c(0, 4))

grid.arrange(
  plot_mat(polymodIT[[1]][1][[1]]$matrix, option = "B", limits = c(0, 3.2)),
  plot_mat(polymodIT[[1]][2][[1]]$matrix, option = "B", limits = c(0, 3.2))
)

grid.arrange(
  plot_mat(polymodIT[[1]][1][[1]]$matrix, option = "B", limits = c(0, 2.3),
           begin = 0.3),
  plot_mat(polymodIT[[1]][2][[1]]$matrix, option = "B", limits = c(0, 2.3),
           begin = 0.3)
)
