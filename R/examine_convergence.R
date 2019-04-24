source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(stringr)
library(ggplot2)
library(reshape2)

# Checking convergence ---------------------------------------------------------

check_conv <- function(i, ...){
  get_data(i)
  # Capture printed output from nlm and put into table ---------------------------
  tmp <- capture.output(res <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
                                   muy = predict(demfit, type = "response"),
                                   N = sum(PS), ..., print = 2))
  tmp <- tolower(tmp)
  # Determine params based on length of starting values for parameters
  no_par <- length(unlist(lapply(str_split(tmp[which(str_detect(tmp, "parameter")) + 1][1], " "),
                          function(x) x[!(x %in% c("[1]", ""))])[[1]]))
  no_par <- as.numeric(no_par)
  
  iter_table <- function(tmp, no_par){
    its <- which(str_detect(tmp, "iteration"))
    stp <- which(str_detect(tmp, "step"))
    par <- which(str_detect(tmp, "parameter"))
    fnc <- which(str_detect(tmp, "function"))
    grad <- which(str_detect(tmp, "gradient"))
    iteration <- seq(from = 1, to = length(its)) - 1
    
    get_vals <- function(input){
      input <- input[1 : stop]
      values <- as.numeric(sapply(str_split(string = tmp[input + 1], pattern = " "), 
                                  function(x) x[[2]]))
      return(values)
    }
    get_vals2 <- function(input){
      input <- input[1 : stop]
      str <- str_split(string = tmp[input + 1], pattern = " ")
      str <- lapply(str, function(x) x[!(x %in% c("[1]", ""))])
      values <- cbind(as.numeric(sapply(str, function(x) x[[1]])),
                      as.numeric(sapply(str, function(x) x[[2]])))
      return(values)
    }
    
    if(no_par == 1){
      stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
      
      step <- get_vals(stp)
      params <- get_vals(par)
      func <- get_vals(fnc)
      gradient <- get_vals(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
    if(no_par == 2){
      stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
      
      step <- get_vals2(stp)
      params <- get_vals2(par)
      func <- get_vals(fnc)
      gradient <- get_vals2(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
  }
  if(no_par == 1)
  {return(as_tibble(cbind(iteration = iter_table(tmp, no_par)$iteration,
                          step = iter_table(tmp, no_par)$step,
                          params = iter_table(tmp, no_par)$params,
                          func = iter_table(tmp, no_par)$func,
                          gradient = iter_table(tmp, no_par)$gradient)
                    [1 : length(which(!is.na(iter_table(tmp, no_par)$gradient))), ]))
  }
  
  if(no_par == 2)
  {return(as_tibble(cbind(iteration = iter_table(tmp, no_par)$iteration,
                          step1 = iter_table(tmp, no_par)$step[, 1],
                          step2 = iter_table(tmp, no_par)$step[, 2],
                          params1 = iter_table(tmp, no_par)$params[, 1],
                          params2 = iter_table(tmp, no_par)$params[, 2],
                          func = iter_table(tmp, no_par)$func,
                          gradient1 = iter_table(tmp, no_par)$gradient[, 1],
                          gradient2 = iter_table(tmp, no_par)$gradient[, 2])
                    [1 : length(which(!is.na(iter_table(tmp, no_par)$gradient[, 1]))), ]))
  }
}

# Example
check_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "constant", startpar = 0.5)
check_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "loglin", startpar = c(0.5, 0.3))

# Plots of convergence for each country
plot_nlm_print <- function(i, ...){
  tmp <- check_conv(i, ...)
  if("params" %in% names(tmp)){
    grid.arrange(ggplot(mapping = aes(x = iteration, y = params, group = 1), 
         data = tmp) + 
    geom_point() + 
    geom_line() +
    labs(title = opts[i, 1]))
  }
  if("params1" %in% names(tmp)){
    grid.arrange(
      ggplot(mapping = aes(x = iteration, y = value, colour = variable), 
             data = melt(tmp, id = "iteration", c("params1", "params2"))) + 
        geom_point() + 
        geom_line() +
        labs(title = opts[i, 1]),
      ggplot(mapping = aes(x = iteration, y = params1, group = 1), 
             data = tmp[, c("iteration", "params1")]) + 
        geom_point(colour = hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[1]) + 
        geom_line(colour = hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[1]),
      ggplot(mapping = aes(x = iteration, y = params2, group = 1), 
             data = tmp[, c("iteration", "params2")]) + 
        geom_point(colour = hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[2]) + 
        geom_line(colour = hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[2]),
      layout_matrix = matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
  }
}

# Example
plot_nlm_print(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "constant", startpar = 0.5)
plot_nlm_print(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "loglin", startpar = c(0.5, 0.3))

rates_conv <- function(i, ...){
  tmp <- check_conv(i, ...)
  # Needed for calculations
  age <- sero$AGE
  y <- sero$indic
  rij <- contact_w
  muy <- predict(demfit, type = "response")
  N <- sum(PS)
  muy <- muy[1 : Lmax]
  L <- Lmax * mean(exp(- cumsum(muy)))
  My <- exp(- cumsum(muy))
  My <- My[1 : Lmax]
  htab <- table(floor(age), y)
  if("params" %in% names(tmp)){
    dat <- sapply(tmp$params, function(qpar){
      bij <- 365 * qpar * (rij)[1 : Lmax, 1 : Lmax]
      R0ij <- (N / L) * D * bij[1 : Lmax, 1 : Lmax]
      Mij <- diag(c(My[1 : Lmax]))
      NGM <- Mij %*% R0ij
      NGM[is.na(NGM)] <- 0 # Replace missings with zeros
      R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
      suscp <- htab[, 1] / rowSums(htab)
      len <- min(dim(NGM)[1], length(suscp))
      Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
      return(list(R0 = max(as.double(R0vec)), R = max(as.double(Rvec))))
    })
  }
  if("params1" %in% names(tmp)){
    dat <- mapply(function(q1, q2){
      qpar <- c(q1, q2)
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
      R0ij <- (N / L) * D * bij[1 : Lmax, 1 : Lmax]
      Mij <- diag(c(My[1 : Lmax]))
      NGM <- Mij %*% R0ij
      NGM[is.na(NGM)] <- 0 # Replace missings with zeros
      R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
      suscp <- htab[, 1] / rowSums(htab)
      len <- min(dim(NGM)[1], length(suscp))
      Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
      return(list(R0 = max(as.double(R0vec)), R = max(as.double(Rvec))))
    }, tmp$params1, tmp$params2)
  }
  return(data.frame(iteration = tmp$iteration,
                    R0 = unlist(t(dat)[, 1]), R = unlist(t(dat)[, 2])))
}

# Example
rates_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "constant", startpar = 0.5)
ggplot(data = rates_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, 
                         prop = "constant", startpar = 0.5),
       mapping = aes(x = iteration, y = R0)) + 
  geom_point() + geom_line()
rates_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "loglin", startpar = c(0.5, 0.3))