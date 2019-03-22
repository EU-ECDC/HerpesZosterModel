#source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/scripts/exmpl_calc.R")
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/scripts/exmpl_calc.R",
                 ssl.verifypeer = FALSE)
eval(parse(text = script))

################################################################################
starts <- seq(0.01, 0.1, length.out = 100)
ests <- suppressMessages(sapply(starts,
                                function(x){run_calc(7, D = 6 / 365, A = 0.5, 
                                                     plots = FALSE, startpar = x)[1, 2]}))
ests <- as.numeric(ests)
plot(starts, ests, main = "q", ylab = "", type = "l", col = 1)
abline(h = mean(ests), col = 2)

# Attempts at parallelising - doesn't currently work :(
library(parallel)
par(mfrow = c(2, 4))
for(i in 1:length(opts)){
  num_cores <- detectCores()
  #num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, list(library(dplyr), library(socialmixr)))
  clusterExport(cl, list("get_eurostat", "example_calc", 
                         "run_calc", "opts", "countries",
                         "starts", "data", "i"))
  ests <- suppressMessages(parSapply(cl, starts,
                                     function(x){run_calc(i, D = 6 / 365, A = 0.5, 
                                                          plots = FALSE, startpar = x)[1, 2]}))
  plot(starts, ests, main = opts[i], ylab = "", type = "l", col = 1)
  stopCluster(cl)
}

#library(doParallel)
##library(doMC)
#library(snow)
#cl <- makeCluster(num_cores, type = "SOCK")
##registerDoMC(num_cores)
#registerDoParallel(cl)
#clusterExport(cl, list("%>%", "get_eurostat", "example_calc", 
#                       "run_calc", "countries"))
#
#ests <- foreach(i = starts) %dopar% 
#  run_calc(7, D = 6 / 365, A = 0.5,
#                       plots = FALSE, startpar = i)[1, 2]
#stopCluster(cl)