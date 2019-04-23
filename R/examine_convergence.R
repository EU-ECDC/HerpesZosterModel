source("https://raw.githubusercontent.com/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/HerpesZosterModel/master/R/model.R")

library(stringr)
library(ggplot2)

# Checking convergence ---------------------------------------------------------

check_conv <- function(i, ...){
  pre_proc(i)
  # Capture printed output from nlm and put into table ---------------------------
  tmp <- capture.output(res <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
                                   muy = predict(demfit, type = "response"),
                                   N = sum(PS), ..., print = 2))
  iter_table <- function(tmp, params = 1){
    tmp <- tolower(tmp)
    its <- which(str_detect(tmp, "iteration"))
    stp <- which(str_detect(tmp, "step"))
    par <- which(str_detect(tmp, "parameter"))
    fnc <- which(str_detect(tmp, "function"))
    grad <- which(str_detect(tmp, "gradient"))
    iteration <- seq(from = 1, to = length(its)) - 1
    if(params == 1){
      stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
      get_vals <- function(input){
        input <- input[1 : stop]
        values <- as.numeric(sapply(str_split(string = tmp[input + 1], pattern = " "), 
                                    function(x) x[[2]]))
        return(values)
      }
      step <- get_vals(stp)
      params <- get_vals(par)
      func <- get_vals(fnc)
      gradient <- get_vals(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
    # Currently only able to use this with constant proportionality factor
    #  if(params > 1){
    #    step <- tmp[stp + 1]
    #    params <- tmp[par + 1]
    #    func <- tmp[fnc + 1]
    #    gradient <- tmp[grad + 1]
    #     TODO For more than one parameter consider something like
    #    get_val <- function(input){
    #      str_tmp <- str_split(string = step, pattern = "\\[")[2]
    #      out <- str_split(string = sapply(str_tmp, function(x) x[2]), pattern = "\\] ")
    #     index <- sapply(out, function(x) x[[1]])
    #      value <- sapply(out, function(x) x[[2]])
    #    }
  }
  return(as_tibble(cbind(iteration = iter_table(tmp)$iteration,
                         step = iter_table(tmp)$step,
                         params = iter_table(tmp)$params,
                         func = iter_table(tmp)$func,
                         gradient = iter_table(tmp)$gradient)[1 : length(which(!is.na(iter_table(tmp)$gradient))), ]))
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

