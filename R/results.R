source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

# Save plot
tiff("S:/HelenJohnson/Herpes Zoster/Force of infection/Figures/overview_all.tif",
     width = 2000, height = 1000)
layout(matrix(seq(1, dim(opts)[1]), nrow = 3, byrow = TRUE))
i <- 1
while(i < dim(opts)[1]){
  get_data(i)
  res1 <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "constant", startpar = 0.5)
  res2 <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
              muy = predict(demfit, type = "response"),
              N = sum(PS), D = 6 / 365, A = 0.5, Lmax = 70, 
              prop = "loglin", startpar = c(0.5, 0.3))
  
  plot(res1$lambda, ylim = c(0, 1), xlab = "age", ylab = "", type = "l", col = 4)
  points(unique(cbind(sero$AGE, res1$pi)), type = "l", col = 4)
  text(length(res1$lambda), tail(res1$lambda, 1) + 0.05, "FOI")#, col = 4)
  text(max(sero$AGE), tail(res1$pi, 1) - 0.05, "prev")#, col = 4)
  if(!(opts[i, 3] %in% c("FI", "UK", "RS"))){
    # Currently experiencing FOI issues for Finland, Serbia, UK
    points(res2$lambda, ylim = c(0, 1), xlab = "age", ylab = "", type = "l")
    points(unique(cbind(sero$AGE, res2$pi)), type = "l")
    #text(length(res2$lambda), tail(res2$lambda, 1) + 0.05, "FOI")
    #text(max(sero$AGE), tail(res2$pi, 1) - 0.05, "prev")
  }
  axis(4, at = pretty(range(c(res1$lambda, res2$lambda))))
  legend("topright", col = c(1, 4), c("Log-linear", "Constant"), lty = 1)
  title(main = opts[i, 1])
  if(i == 8){
    text(20, 0.4, "NB Serbia's contact matrix is interpolated")
  }
  i <- i + 1
}
dev.off()