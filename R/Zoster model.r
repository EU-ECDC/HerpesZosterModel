# Load required packaged
library(deSolve)
library(ggplot2)

# Function to solve model given initial conditions
model <- function(initial_states, time, parameters, func = func,
                  event = FALSE, time2 = NULL, func2 = lambda){
  if(!is.logical(event))
    stop("event must be TRUE or FALSE")
  if(event == FALSE){
    out <- lsoda(y = initial_states, times = time, func = func, parms = parameters)
    out <- as.data.frame(out)}
  if(event == TRUE){
    if(is.null(time2)){
      stop("event times must be given")
    }
    out <- lsoda(y = initial_states, times = time, func = func, parms = parameters,
                 events = list(func = func2, time = time2))
    out <- as.data.frame(out)}
  return(out)}

# Garnett and Grenfell
gg <- function(time, state, params){
  with(as.list(c(state, params)), {
    # Maternal protection
    dM <- - (delta + mu)
    # Susceptible
    dS <- delta * M - (lambda + mu) * S
    # Exposed (latent)
    dE <- lambda * S - (sigma + mu) * E
    # Infectious
    dI <- sigma * E - (nu + mu) * I
    # Virus carriers
    dV <- - mu * rho * V
    # Zoster
    #dZ <- integrate(V * rho - (mu + gamma) * Z, lower = 0, upper = Inf)
    
    # Force of infection
    lambda <- beta * I
    
    out <- list(c(dM, dS, dE, dI, dV), c(lambda = lambda, mu = mu))
    return(out)})}

t <- 50
time <- seq(0, t, by = t / (2 * length(1 : t)))

initial_states <- c(M = integrate(function(a){1 / ((a + 1) * sqrt(a))},
                                  #function(a){N(a) * mu(a)}, 
                                  # Placeholder from help page until we have data
                                  lower = 0, upper = Inf)$value,
                    S = 0, E = 0, I = 0, V = 0)
parameters <- c(lambda = 0.1,
                beta = 0.5,
                delta = 0.74, 
                mu = 0.14,
                sigma = 0.01,
                nu = 0.61,
                rho = 1.18,
                gamma = 0.45,
                rho = 0.345)

ggplot(data = model(initial_states, time, parameters, func = gg),
       mapping = aes(x = time, y = S)) + geom_line()

# Horn et al
horn <- function(time, state, params){
  with(as.list(c(state, params)), {
    # Maternal protection
    dM <- - b * M + v * M
    # Susceptible
    dS <- b * M - (lambda + v) * S
    # Exposed (latent)
    dE <- lambda * S - c * E
    # Infectious
    dI <- c * E - d * I
    # Resistant regarding natural varicella
    dR <- d * I + e * lambda * SZN - f * R
    # Susceptible vaccinated regarding herpes zoster (natural varicella)
    dSZN <- - e * lambda * SZN - g * SZN + h * RZN + f * R
    # Infectious vaccinated regarding herpes zoster (natural varicella)
    dIZN <- g * SZN - i * IZN
    # Resistant vaccinated regarding herpes zoster (natural varicella)
    dRZN <- i * IZN - h * RZN
    # Susceptible regarding breakthrough varicella
    dSB <- v1e * v * (M + S) - (f + lambda) * SB + vw1 * V1 + vw2 * V2
    # Susceptible to herpes zoster due to vaccine virus 
    dSBV <- - lambda * SBV + f * SB - k * g * SBV + vw1 * V1V + vw2 * V2V
    # Vaccinated against varicella with 1 dose
    dV1 <- (1 - v1e) * v * (M + S) - (f + vw1 + lambda) * V1
    # Vaccinated against varicella with 2 doses
    dV2 <- (1 - v2e) * v2 * (SB + V1) - (f + vw1 + lambda) * V2
    # Varicella vaccine-related protection boosted to lifelong duration
    dVL <- lambda * (V1 + V2 + V1V + V2V + SZV) - f * VL
    # Exposed regarding breakthrough varicella;
    dEB <- lambda * (SB + SBV) - c * EB
    # Infectious regarding breakthrough varicella;
    dIB <- c * EB - d2 * IB
    # Resistant regarding breakthrough varicella;
    dRB <- d2 * IB + e * lambda * SZB - f * RB
    # Susceptible regarding herpes zoster (breakthrough varicella)
    dSZB <- - e * lambda * SZB - j * g * SZB + h * RZB + f * RB
    # Infectious regarding herpes zoster (breakthrough varicella)
    dIZB <- j * g * SZB - i * IZB
    # Resistant regarding herpes zoster (breakthrough varicella)
    dRZB <- i * IZB - h * RZB
    # Vaccinated against varicella (with one dose) and susceptible to HZ due to vaccine virus
    dV1V <- f * V1 - (k * g + lambda + vw1) * V1V
    # Vaccinated against varicella (with two doses) and susceptible to HZ due to vaccine virus
    dV2V <- f * V2 - (k * g + lambda + vw2) * V2V
    # Susceptible regarding herpes zoster (vaccine virus)
    dSZV <- e * lambda * SZV - k * g * SZV + h * RZV + f * VL
    # Infectious regarding herpes zoster (vaccine virus)
    dIZV <- k * g * (SBV + V1V + V2V + SZV) - i * IZV
    # Resistant regarding herpes zoster (vaccine virus)
    dRZV <- i * IZV - h * RZV
    # Combined in a list
    out <- list(c(dM, dS, dE, dI, dR, dSZN, dIZN, dRZN, dSB, dSBV, dV1, dV2, dVL,
                 dEB, dIB, dRB, dSZB, dIZB, dRZB, dV1V, dV2V, dSZV, dIZV, dRZV))
    return(out)})}

# Example without force of infection function (values chosen arbitrarily)
t <- 50
time <- seq(0, t, by = t / (2 * length(1 : t)))
initial_states <- c(M = 164, S = 38, E = 115, I = 206, R = 42, SZN = 65,
                    IZN = 3, RZN = 29, SB = 36, SBV = 44, V1 = 170,
                    V2 = 109, VL = 129, EB = 163, IB = 29, RB = 108,
                    SZB = 129, IZB = 65, RZB = 63, V1V = 31, V2V = 11, 
                    SZV = 156, IZV = 8, RZV = 5)
parameters <- c(b = 0.74, v = 0.14, lambda = 0.01, c = 0.61, d = 1.18, e = 1.09,
                f = 1.06, g = 0.02, h = 1.42, i = 1.26, v1e = 0.69, vw1 = 0.47,
                vw2 = 0.45, k = 3.32, v2e = 1.13, v2 = 0.11, d2 = 2.59, j = 0.43)

ggplot(data = model(initial_states, time, parameters, func = horn),
       mapping = aes(x = time, y = S)) + geom_line()
