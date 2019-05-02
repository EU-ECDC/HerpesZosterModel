source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(ggplot2)
theme_set(theme_classic() %+replace%
            theme(plot.title = element_text(hjust = 0.5))) # Ensure centred titles
library(gridExtra)

# Plot
plot_results <- function(code, ...){
  get_data(code)
  res1 <- FOI(age = seroData$AGE, y = seroData$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(popSize), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "constant", startpar = 0.5)
  res2 <- FOI(age = seroData$AGE, y = seroData$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(popSize), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "loglin", startpar = c(0.2, 0.3))
  
  seroData <- seroData[(seroData$AGE > 0.5) & (seroData$AGE < 80) &
                 (!is.na(seroData$AGE)) & !is.na(seroData$indic), ]
  htab <- table(floor(seroData$AGE[order(seroData$AGE)]), 
                seroData$indic[order(seroData$AGE)])
  
  ggplot() +
    geom_point(data = as.data.frame(cbind(age = as.numeric(row.names(htab)), 
                                          prop = htab[, 2] / rowSums(htab),
                                          tot = rowSums(htab))),
               mapping = aes(x = age, y = prop, size = tot),
               pch = 1) +
    geom_line(data = data.frame(age = 1 : length(res1$lambda), lambda = res1$lambda),
              mapping = aes(x = age, y = lambda), colour = 4) +
    geom_line(data = data.frame(unique(cbind(seroData$AGE, res1$pi))),
              mapping = aes(x = X1, y = X2), colour = 4) +
    geom_line(data = data.frame(age = 1 : length(res2$lambda), lambda = res1$lambda),
              mapping = aes(x = age, y = lambda), colour = 1) +
    geom_line(data = data.frame(unique(cbind(seroData$AGE, res2$pi))),
              mapping = aes(x = X1, y = X2), colour = 1) +
    annotate(geom = "text", x = res1$inputs$Lmax, y = 0.9, 
             label = "Constant", colour = 4) +
    annotate(geom = "text", x = res2$inputs$Lmax, y = 0.8, 
             label = "Log-linear") +
    annotate(geom = "text", x = 20, y = 0.7, 
             label = paste("Sum of abs diff:",
                           round(sum(abs(res1$lambda[as.numeric(names(htab[, 2] / 
                                                                  rowSums(htab)))] - 
                                     (htab[, 2] / rowSums(htab)))), 4)),
             colour = 4) +
    annotate(geom = "text", x = 20, y = 0.6, 
             label = paste("Sum of abs diff:",
                           round(sum(abs(res2$lambda[as.numeric(names(htab[, 2] / 
                                                                  rowSums(htab)))] - 
                                     (htab[, 2] / rowSums(htab)))), 4))) +
    labs(title = code, x = "age", y = "") + 
    theme(legend.position = "none") +
    if(code == "RS"){
      annotate(geom = "text", x = 20, y = 0.4, label = "NB Serbia's contact matrix is interpolated")
    } else {
      NULL
    }
}
plot_list <- lapply(use, plot_results)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

# Save plot
tiff("S:/HelenJohnson/Herpes Zoster/Figures/overview_all.tif",
     width = 1400, height = 800)
# Current options based on availability of data
grid.arrange(BE, FI, DE, IE, IT, LU, NL, SK, UK, RS, SI)
while(!is.null(dev.list())) dev.off()

