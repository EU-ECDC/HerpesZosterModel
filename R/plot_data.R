source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
# Plots of data ----------------------------------------------------------------

no_axes <- theme(#axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  #axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())

no_x_axis <- theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank())

no_y_axis <- theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank())

## Plot mortality
plot_mort <- function(code, ...){
  get_data(code)
  plotGAM(gamFit = demfit, smooth.cov = "mortAge") +
    labs(title = code, y = "Predicted deaths", x = "Age (smoothed)") +
    theme_classic() + theme(plot.title = element_text(hjust = 1,
                                                      margin = margin(t = 10, b = - 20)),
                            title = element_text(family = "serif")) +
    scale_y_continuous(breaks = seq(- 10, 0, 2.5),
                       labels = seq(- 10, 0, 2.5)) + 
    ylim(- 10, 0)
}

# Save list of plots for the countries in use
plot_list <- lapply(use, plot_mort)
# Save the plots as their country codes
lapply(seq_along(1 : length(use)), function(x){assign(use[x], plot_list[[x]], 
                                                      envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/mortality.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE + no_x_axis, FI + no_axes, DE + no_axes, 
             IE + no_x_axis, IT + no_axes, LU + no_axes,
             NL + no_x_axis, SK + no_axes, UK + no_axes, 
             RS, SI + no_y_axis,
             ncol = 3, nrow = 4)
# TODO see if we can replace ^ with use somehow
while(!is.null(dev.list())) dev.off()

## Plot deaths
plot_death <- function(code, ...){
  get_data(code)
  ggplot(mapping = aes(x = age, y = rate), 
         data = as.data.frame(cbind(age = 1 : length(popSize),
                                    rate = demfit$model$nDeaths / popSize * 1e+05))) + 
    labs(x = "Age", y = "Deaths per 100,000", title = code) +
    geom_line() +
    theme(plot.title = element_text(hjust = 1, # Ensure right-aligned titles
                                    margin = margin(t = 10, b = - 20)), # Title in panel
          title = element_text(family = "serif")) +
    scale_y_continuous(breaks = seq(5000, 5e04, 5000),
                       labels = seq(5000, 5e04, 5000)) + 
    ylim(0, 5e04)
}

# Save list of plots for the countries in use
plot_list <- lapply(use, plot_death)
# Save the plots as their country codes
lapply(seq_along(1 : length(use)), function(x){assign(use[x], plot_list[[x]], 
                                                      envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/deaths.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE + no_x_axis, FI + no_axes, DE + no_axes, 
             IE + no_x_axis, IT + no_axes, LU + no_axes,
             NL + no_x_axis, SK + no_axes, UK + no_axes, 
             RS, SI + no_y_axis,
             ncol = 3, nrow = 4)
while(!is.null(dev.list())) dev.off()

## Plot serological data
plot_sero <- function(code, ...){
  get_data(code)
  seroData <- seroData[(seroData$AGE > 0.5) & (seroData$AGE < 80) &
                         (!is.na(seroData$AGE)) & !is.na(seroData$indic), ]
  # Table with seroprevalence for each single year of age
  htab <- table(floor(seroData$AGE[order(seroData$AGE)]), 
                seroData$indic[order(seroData$AGE)])
  
  ggplot(data = as.data.frame(cbind(age = as.numeric(row.names(htab)), 
                                    prop = htab[, 2] / rowSums(htab),
                                    tot = rowSums(htab))),
         mapping = aes(x = age, y = prop, size = tot)) +
    geom_point(pch = 1) + 
    labs(x = "Age", y = "Sero-prevalence", title = code) + 
    xlim(0, 72) + ylim(- 0.1, 1) + 
    theme(legend.title = element_blank()) +
    scale_size_continuous(limits = c(1, 500),
                          breaks=c(50, 100, 150, 200,
                                   250, 300, 400, 500)) +
    theme(plot.title = element_text(hjust = 1,
                                    margin = margin(t = 10, b = - 20)),
          title = element_text(family = "serif"))
}

# Save list of plots and assign to the environment
plot_list <- lapply(use, plot_sero)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/serology.tif",
     width = 800, height = 800)
# Current options based on availability of data
grid_arrange_shared_legend(BE + no_x_axis, FI + no_axes, DE + no_axes, 
                           IE + no_x_axis, IT + no_axes, LU + no_axes,
                           NL + no_x_axis, SK + no_axes, UK + no_axes, 
                           RS, SI + no_y_axis,
                           ncol = 3, nrow = 4)
while(!is.null(dev.list())) dev.off()

## Plot contact rates
# Adapted from https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/old_scripts/plots/plot_polymod_matrix.R
plot_mat <- function(code, ...){
  get_data(code)
  tmp <- melt(contact_w)
  # Relabel socialmixr for panel plot
  if(is.factor(tmp$age.group)){
    levels(tmp$age.group) <- 1 : length(levels(tmp$age.group))
    tmp$age.group <- as.numeric(levels(tmp$age.group))[tmp$age.group]
  }
  if(is.factor(tmp$contact.age.group)){
    levels(tmp$contact.age.group) <- 1 : length(levels(tmp$contact.age.group))
    tmp$contact.age.group <- as.numeric(levels(tmp$contact.age.group))[tmp$contact.age.group]
  }
  ggplot(data = tmp, aes_string(x = names(tmp)[1], y = names(tmp)[2], 
                                fill = names(tmp)[3])) + 
    geom_tile() + 
    scale_fill_viridis_c(direction = - 1, option = "E", 
                         guide = guide_legend(label.theme = element_text(angle = 90)),
                         ...) +
    labs(x = "", y = "", title = code, fill = "Rate") +
    theme(plot.title = element_text(hjust = 1,
                                    margin = margin(t = 10, b = - 20)),
          title = element_text(family = "serif")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlim(0, 100) + ylim(0, 100)
}
plot_list <- lapply(use, plot_mat)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/contact_matrices.tif",
     width = 900, height = 900)
# Current options based on availability of data
grid_arrange_shared_legend(BE + no_x_axis, FI + no_axes, DE + no_axes, 
                           IE + no_x_axis, IT + no_axes, LU + no_axes,
                           NL + no_x_axis, SK + no_axes, UK + no_axes, 
                           RS, SI + no_y_axis,
                           ncol = 3, nrow = 4)
while(!is.null(dev.list())) dev.off()

## Plot population
plot_pop <- function(code, ...){
  get_data(code)
  ggplot(data = data.frame(age = 1 : length(popSize), popSize), 
         mapping = aes(x = age, 
                       y = popSize)) +
    geom_col(width = 1) +
    labs(y = "Population", x = "Age", title = code) +
    scale_y_continuous(breaks = seq(2500, 3e05, 5000),
                       labels = seq(2500, 3e05, 5000)) +
    ylim(0, 3e05) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 1,
                                    margin = margin(t = 10, b = - 20)),
          title = element_text(family = "serif"))
}
plot_list <- lapply(use, plot_pop)
lapply(seq_along(plot_list), function(x){assign(use[x], plot_list[[x]], 
                                                envir = .GlobalEnv)})

tiff(filename = "S:/HelenJohnson/Herpes Zoster/Figures/population.tif",
     width = 800, height = 600)
# Current options based on availability of data
grid.arrange(BE + no_x_axis, FI + no_axes, DE + no_axes, 
             IE + no_x_axis, IT + no_axes, LU + no_axes,
             NL + no_x_axis, SK + no_axes, UK + no_axes, 
             RS, SI + no_y_axis,
             ncol = 3, nrow = 4)
while(!is.null(dev.list())) dev.off()
