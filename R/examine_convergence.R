source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/load_data.R")
source("https://raw.githubusercontent.com/EU-ECDC/HerpesZosterModel/master/R/model.R")

library(stringr)
library(ggplot2)
theme_set(theme_classic() %+replace%
            theme(plot.title = element_text(hjust = 0.5))) # Ensure centred titles
library(reshape2)
library(pbapply)
library(beepr)
library(gridExtra)

# Checking convergence ---------------------------------------------------------

check_conv <- function(code, ...){
  get_data(code)
  # Capture printed output from nlm and put into table ---------------------------
  # Capture the output and make it lowercase for regex
  tmp <- capture.output(res <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
                                   muy = predict(demfit, type = "response"),
                                   N = sum(popSize), ..., print = 2))
  tmp <- tolower(tmp)
  # Determine params based on length of starting values for parameters
  no_par <- length(unlist(lapply(str_split(tmp[which(str_detect(tmp, "parameter")) + 1][1], " "),
                                 function(x) x[!(x %in% c("[1]", ""))])[[1]]))
  no_par <- as.numeric(no_par)
  
  iter_table <- function(tmp, no_par){
    # Locate the following terms in the captured text
    its <- which(str_detect(tmp, "iteration"))
    stp <- which(str_detect(tmp, "step"))
    par <- which(str_detect(tmp, "parameter"))
    fnc <- which(str_detect(tmp, "function"))
    grad <- which(str_detect(tmp, "gradient"))
    iteration <- seq(from = 1, to = length(its)) - 1 # it starts at 0
    
    # Stop value to help us obtain matching length of inputs in the tables returned
    stop <- min(c(length(stp), length(its), length(par), length(fnc), length(grad)))
    get_vals <- function(input){
      input <- input[1 : stop]
      # Variables with one value
      values <- as.numeric(sapply(str_split(string = tmp[input + 1], pattern = " "), 
                                  function(x) x[[2]])) 
      return(values)
    }
    get_vals2 <- function(input){
      input <- input[1 : stop]
      # Variables with two values
      str <- str_split(string = tmp[input + 1], pattern = " ")
      str <- lapply(str, function(x) x[!(x %in% c("[1]", ""))])
      values <- cbind(as.numeric(sapply(str, function(x) x[[1]])),
                      as.numeric(sapply(str, function(x) x[[2]])))
      return(values)
    }
    # Get the values for the proportionality factor with one parameter
    if(no_par == 1){
      step <- get_vals(stp)
      params <- get_vals(par)
      func <- get_vals(fnc)
      gradient <- get_vals(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
    # Get the values for the proportionality factor with two parameters
    if(no_par == 2){
      step <- get_vals2(stp)
      params <- get_vals2(par)
      func <- get_vals(fnc) # Function does not return two values
      gradient <- get_vals2(grad)
      return(list(iteration = iteration, step = step, params = params,
                  func = func, gradient = gradient))
    }
  }
  # Save the values as tables
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
check_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
check_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "loglin", startpar = c(0.5, 0.3))

# Plots of convergence for each country
plot_conv <- function(code, ...){
  # Get the table of convergence
  tmp <- check_conv(code, ...)
  if("params" %in% names(tmp)){ # One parameter
    grid.arrange(ggplot(mapping = aes(x = iteration, y = params, group = 1), 
                        data = tmp) + 
                   geom_point() + 
                   geom_line() +
                   labs(title = code))
  }
  if("params1" %in% names(tmp)){ # Two parameters
    grid.arrange(
      ggplot(mapping = aes(x = iteration, y = value, colour = variable), 
             data = melt(tmp, id = "iteration", c("params1", "params2"))) + 
        geom_point() + 
        geom_line() +
        labs(title = code),
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
plot_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
plot_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "loglin", startpar = c(0.5, 0.3))

rates_conv <- function(code, ...){
  # Get the table of convergence
  tmp <- check_conv(code, ...) # This runs get_data so objects below exist
  # Calculations of rates - inner workings from the model
  ## Data
  age <- sero$AGE
  y <- sero$indic
  rij <- contact_w
  muy <- predict(demfit, type = "response")
  N <- sum(popSize)
  ## Values as provided in model
  inputs <- FOI(age, y, rij, muy, N, ...)$inputs
  Dur <- inputs$Dur
  A <- inputs$A
  Lmax <- inputs$Lmax
  ## Inner workings from model
  muy <- muy[1 : Lmax]
  L <- Lmax * mean(exp(- cumsum(muy)))
  My <- exp(- cumsum(muy))
  My <- My[1 : Lmax]
  htab <- table(floor(age), y)
  # For each value of the parameter(s) estimated during the minimisation procedure
  # calculate the corresponding R0 and R values
  if("params" %in% names(tmp)){ # One parameter
    dat <- sapply(tmp$params, function(qpar){
      bij <- 365 * qpar * (rij)[1 : Lmax, 1 : Lmax]
      R0ij <- (N / L) * Dur * bij[1 : Lmax, 1 : Lmax]
      Mij <- diag(c(My[1 : Lmax]))
      NGM <- Mij %*% R0ij
      NGM[is.na(NGM)] <- 0 # Replace missings with zeros
      R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
      suscp <- htab[, 1] / rowSums(htab)
      len <- min(dim(NGM)[1], length(suscp))
      Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
      return(list(R0 = max(as.double(R0vec)), R = max(as.double(Rvec))))
    })
    return(data.frame(iteration = tmp$iteration,
                      R0 = unlist(t(dat)[, 1]), R = unlist(t(dat)[, 2]),
                      params = tmp$params))
  }
  if("params1" %in% names(tmp)){ # Two parameters
    dat <- mapply(function(q1, q2){
      qpar <- c(q1, q2)
      q.f <- function(x, y){exp(qpar[1] + qpar[2] * x)}
      qij <- outer(c(1 : Lmax), c(1 : Lmax), q.f)
      bij <- 365 * qij * (rij)[1 : Lmax, 1 : Lmax]
      R0ij <- (N / L) * Dur * bij[1 : Lmax, 1 : Lmax]
      Mij <- diag(c(My[1 : Lmax]))
      NGM <- Mij %*% R0ij
      NGM[is.na(NGM)] <- 0 # Replace missings with zeros
      R0vec <- eigen(NGM, symmetric = FALSE, only.values = TRUE)$values
      suscp <- htab[, 1] / rowSums(htab)
      len <- min(dim(NGM)[1], length(suscp))
      Rvec <- NGM[1 : len, 1 : len] %*% diag(suscp[1 : len])
      return(list(R0 = max(as.double(R0vec)), R = max(as.double(Rvec))))
    }, tmp$params1, tmp$params2)
    return(data.frame(iteration = tmp$iteration,
                      R0 = unlist(t(dat)[, 1]), R = unlist(t(dat)[, 2]),
                      params1 = tmp$params1, params2 = tmp$params2))
  }
}
# TODO Find out why this only works when Lmax already defined in the global environment

# Example
rates_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
ggplot(data = rates_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, 
                         prop = "constant", startpar = 0.5),
       mapping = aes(x = iteration, y = R0)) + 
  geom_point() + geom_line()
rates_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "loglin", startpar = c(0.5, 0.3))

# Look at all the R and R0 values through loops which run rates_conv
# for all the countries we have available data for and save this as capt

## One parameter
capt <- lapply(1 : length(use), function(x){
  tryCatch(rates_conv(use[x], Dur = 6 / 365, A = 0.5, Lmax = 70, 
                      prop = "constant", startpar = 0.5), 
           error = function(e) NULL)}) # Ignore errors for now
names(capt) <- use
capt

## Two parameters
capt <- lapply(1 : length(use), function(x){
  tryCatch(rates_conv(use[x], Dur = 6 / 365, A = 0.5, Lmax = 70, 
                      prop = "loglin", startpar = c(0.5, 0.3)), 
           error = function(e) NULL)}) # Ignore errors for now
names(capt) <- use
capt

# Examining various starting parameters ----------------------------------------
# Generate random starting parameters for age-dependent proportionality factor
set.seed(13)
# All combinations of two values from n generated uniformly distributed random
# variables between 0 and 1
vals <- t(combn(runif(n = 4, min = 0, max = 1), 2))
# Adding previously used value as a sanity check
vals <- rbind(vals, c(0.2, 0.1))
vals <- rbind(vals, c(0.5, 0.3))

# Save the outputs 
capt <- sapply(1 : length(use), function(y){
  sapply(1 : dim(vals)[1], function(x){
    tryCatch(cbind(rates_conv(code = use[y], Dur = 6 / 365, A = 0.5, Lmax = 70, 
                              prop = "loglin", startpar = c(vals[x, 1], vals[x, 2])),
                   id = rep(use[y], 
                            dim(rates_conv(code = use[y], Dur = 6 / 365, A = 0.5, Lmax = 70,
                                           prop = "loglin", startpar = c(vals[x, 1], vals[x, 2])))[1])),
             error = function(e) NULL)})}) # Ignore errors for now
names(capt) <- rep(1 : (length(capt) / length(use)),
                   times = length(use))
colnames(capt) <- use
capt

# Rework data into a format that is useful
dat <- bind_rows(capt, .id = "column_label")

ggplot(data = dat,
       mapping = aes(x = iteration, y = R0, 
                     group = column_label)) +
  geom_line() + facet_wrap(. ~ id, scales = "free")

# Compare inputs and outputs (starting values and estimates)
run_model <- function(code, ...){
  res <- FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
             muy = predict(demfit, type = "response"),
             N = sum(popSize), ...)
  if(length(res$inputs$start) == 1){ # One parameter
    out <- list(start = res$inputs$start,
                est = res$qhat,
                est_R0 = res$R0,
                est_R = res$R,
                id = code)
  }
  if(length(res$inputs$start) == 2){ # Two parameters
    out <- list(start_gamma1 = res$inputs$start[1],
                start_gamma2 = res$inputs$start[2], 
                est_gamma1 = res$qhat[1],
                est_gamma2 = res$qhat[2],
                est_R0 = res$R0,
                est_R = res$R,
                id = code)
  }
  return(out)
}
# Example
capt <- lapply(1 : dim(vals)[1], function(x){
  tryCatch(t(run_model(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, 
                       startpar = c(vals[x, 1], vals[x, 2]), prop = "loglin")),
           error = function(e) NULL)}) # Ignore errors
## Turn into table
nam <- colnames(capt[[1]]) # Save the names
capt <- matrix(unlist(capt), ncol = 7, byrow = TRUE) # Turn the list into a matrix
# Note that ncol is 7 because we are considering prop = "loglin"
colnames(capt) <- nam # Add the names back
capt <- as.data.frame(capt)
# Ensure all columns except ID are numeric
capt[, - 7] %<>% lapply(function(x) as.numeric(levels(x))[x])

# Reformat for use with ggplot2
tmp <- data.frame(start = c(capt[, 1], capt[, 2]),
                  est = c(capt[, 3], capt[, 4]),
                  par = rep(c("gamma1", "gamma2"), c(dim(capt)[1], dim(capt)[1])),
                  R0 = c(capt[, 5], capt[, 5]),
                  R = c(capt[, 6], capt[, 6]),
                  id = c(capt[, 7], capt[, 7]))
# Plot
ggplot() +
  geom_count(data = tmp,
             mapping = aes(x = start, y = est, colour = par), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Estimated value", x = "Starting value")

ggplot() +
  geom_histogram(data = tmp,
                 mapping = aes(x = R0))
ggplot() +
  geom_histogram(data = tmp,
                 mapping = aes(x = R))

# Do the same as above for all countries
## Remove "pb" if progress bar not desired
{capt <- pbsapply(1 : length(use), function(y){lapply(1 : dim(vals)[1], function(x){
  tryCatch(t(run_model(code = use[y], Dur = 6 / 365, A = 0.5, Lmax = 70, 
                       startpar = c(vals[x, 1], vals[x, 2]), prop = "loglin")),
           error = function(e) NULL)})}) # Ignore errors
# Uncomment below if you want the script to announce when done
#beep(sound = 4, expr = NULL)
}
nam <- colnames(capt[[1]])
capt <- matrix(unlist(capt), ncol = 7, byrow = TRUE)
colnames(capt) <- nam
capt <- as.data.frame(capt)
capt[, - 7] %<>% lapply(function(x) as.numeric(levels(x))[x])

dat <- data.frame(start = c(capt[, 1], capt[, 2]),
                  est = c(capt[, 3], capt[, 4]),
                  par = rep(c("gamma1", "gamma2"), c(dim(capt)[1], dim(capt)[1])),
                  R0 = c(capt[, 5], capt[, 5]),
                  R = c(capt[, 6], capt[, 6]),
                  id = levels(capt[, 7])[c(capt[, 7], capt[, 7])])

(fig <- ggplot() +
  geom_count(data = dat,
             mapping = aes(x = start, y = est, colour = par), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Estimated value", x = "Starting value") +
  facet_grid(id ~ par, scales = "free_y"))

tiff("P:/temp.tif",
     width = 800, height = 1000)
grid.arrange(fig,
             ggplot() +
               geom_count(data = dat,
                          mapping = aes(x = start, y = est, colour = par), alpha = 0.3) +
               geom_abline(intercept = 0, slope = 1) +
               labs(y = "Estimated value", x = "Starting value"), ncol = 2)
dev.off()

# Plot logs of R and R0
grid.arrange(ggplot(data = dat,
                    mapping = aes(x = log(R0), fill = id)) +
               geom_histogram() + scale_fill_viridis_d(direction = - 1) +
               theme_linedraw(),
             ggplot(data = dat,
                    mapping = aes(x = log(R), fill = id)) +
               geom_histogram() + scale_fill_viridis_d(direction = - 1) +
               theme_linedraw())
