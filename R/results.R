source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

# Save plot
tiff("S:/HelenJohnson/Herpes Zoster/Figures/overview_all.tif",
     width = 2000, height = 1000)
layout(matrix(seq(1, dim(opts)[1]), nrow = 3, byrow = TRUE))
i <- 1
while(i <= dim(opts)[1]){
  get_data(i)
  res1 <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "constant", startpar = 0.5)
  res2 <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "loglin", startpar = c(0.5, 0.3))
  
  # Plot results
  plot(res1$lambda, ylim = c(0, 1), xlab = "age", ylab = "", type = "l", col = 4)
  points(unique(cbind(sero$AGE, res1$pi)), type = "l", col = 4)
  text(length(res1$lambda), tail(res1$lambda, 1) + 0.05, "FOI")#, col = 4)
  text(max(sero$AGE), tail(res1$pi, 1) - 0.05, "prev")#, col = 4)
  points(res2$lambda, ylim = c(0, 1), xlab = "age", ylab = "", type = "l")
  points(unique(cbind(sero$AGE, res2$pi)), type = "l")
  # Add seroprevalence points
  subset <- (sero$AGE > 0.5) & (sero$AGE < 80) &
    (!is.na(sero$AGE)) & !is.na(sero$indic)
  sero <- sero[subset, ]
  y <- sero$indic[order(sero$AGE)]
  a <- sero$AGE[order(sero$AGE)]
  
  grid <- sort(unique(round(a)))
  neg <- table(y, round(a))[1, ]
  pos <- table(y, round(a))[2, ]
  tot <- neg + pos
  points(grid, pos / tot, cex = 0.02 * tot)
  # Add right axis, legend, and title
  axis(4, at = pretty(range(c(res1$lambda, res2$lambda))))
  legend("topright", col = c(1, 4), c("Log-linear", "Constant"), lty = 1)
  title(main = opts[i, 1])
  
  text(20, 0.7, paste("Sum of abs diff:",
                      sum(abs(res1$lambda[as.numeric(names(pos / tot))] - (pos / tot)))),
       col = 4)
  text(20, 0.6, paste("Sum of abs diff:",
                      sum(abs(res2$lambda[as.numeric(names(pos / tot))] - (pos / tot)))))
  
  if(i == 8){
    text(20, 0.4, "NB Serbia's contact matrix is interpolated")
  }
  i <- i + 1
}
dev.off()

if(FALSE){ #TODO ggplot2 version
  ggplot() +
    geom_point(data = as.data.frame(cbind(grid, neg, pos, tot)), 
               mapping = aes(x = grid, y = pos / tot, size = 0.02 * tot),
               pch = 1) +
    geom_line(data = data.frame(age = 1 : length(res1$lambda), lambda = res1$lambda),
              mapping = aes(x = age, y = lambda), colour = 4) +
    geom_line(data = data.frame(unique(cbind(sero$AGE, res1$pi))),
              mapping = aes(x = X1, y = X2), colour = 4) +
    geom_line(data = data.frame(age = 1 : length(res2$lambda), lambda = res1$lambda),
              mapping = aes(x = age, y = lambda), colour = 1) +
    geom_line(data = data.frame(unique(cbind(sero$AGE, res2$pi))),
              mapping = aes(x = X1, y = X2), colour = 1) +
    labs(title = opts[i, 1])
}
