source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(stringr)
library(ggplot2)

# Checking convergence ---------------------------------------------------------

check_conv <- function(i, ...){
  get_data(i)
  # Capture printed output from nlm and put into table ---------------------------
  tmp <- capture.output(res <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
                                   muy = predict(demfit, type = "response"),
                                   N = sum(PS), ..., print = 2))
  params <- length(startpar)
  iter_table <- function(tmp, params){
    tmp <- tolower(tmp)
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
    
    if(params == 1){
      stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
      
      step <- get_vals(stp)
      params <- get_vals(par)
      func <- get_vals(fnc)
      gradient <- get_vals(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
    if(params == 2){
      stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
      
      step <- get_vals2(stp)
      params <- get_vals2(par)
      func <- get_vals(fnc)
      gradient <- get_vals2(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
  }
  if(params == 1)
  {return(as_tibble(cbind(iteration = iter_table(tmp, params)$iteration,
                          step = iter_table(tmp, params)$step,
                          params = iter_table(tmp, params)$params,
                          func = iter_table(tmp, params)$func,
                          gradient = iter_table(tmp, params)$gradient)
                    [1 : length(which(!is.na(iter_table(tmp, params)$gradient))), ]))
  }
  
  if(params == 2)
  {return(as_tibble(cbind(iteration = iter_table(tmp, params)$iteration,
                          step1 = iter_table(tmp, params)$step[, 1],
                          step2 = iter_table(tmp, params)$step[, 2],
                          params1 = iter_table(tmp, params)$params[, 1],
                          params2 = iter_table(tmp, params)$params[, 2],
                          func = iter_table(tmp, params)$func,
                          gradient1 = iter_table(tmp, params)$gradient[, 1],
                          gradient2 = iter_table(tmp, params)$gradient[, 2])
                    [1 : length(which(!is.na(iter_table(tmp, params)$gradient[, 1]))), ]))
  }
}

# Example
check_conv(i, D = 6 / 365, A = 0.5, Lmax = 80, prop = "constant", startpar = 0.5)

# Plots of convergence for each country
plot_nlm_print <- function(i, ...){
  ggplot(mapping = aes(x = iteration, y = params, group = 1), data = check_conv(i, ...)) + 
    geom_point() + 
    geom_line() +
    labs(title = opts[i, 1])
}

