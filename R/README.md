# Baseline varicella force of infection

The suggested order of running these scripts is
1. [load_data.R](load_data.R)
2. [model.R](model.R)
  1. [examine_convergence.R](examine_convergence.R)
3. [results.R](results.R)

The script for loading data [1](load_data.R) creates the inputs required for use in the force of infection calculation. The calculation itself can be found in the script called model [2](model.R). Convergence can be examined [2a](examine_convergence.R). Results from running the model after loading the data are also available [3](results.R).

## load_data.R
The function `get_data` retrieves and attaches to the global envioronment 
- The serological data
- The contact matrix
- Fitted mortality data
- Population size
for a given country. To retrieve data for Italy we call 
```R
get_data(code = "IT")
```
The script further contains plots of the serological data, the mortality data, the population data, and the contact matrices for the countries considered.

## model.R
The function `FOI` calculates the force of infection. It takes the following arguments

| Variable | Description | Default value |
|:--------:|-------------|---------------|
| `age` | Age of serological information | |
| `y` | Serological information given as indicator variable | |
| `rij` | Contact matrix | |
| `muy` | Mortality | |
| `N` | Population size at demographic equilibrium | |
| `Dur` | Duration | |
| `Lmax` | Life expectancy | |
| `A` | Duration of maternal immunity | |
| `startpar` | Starting parameters | |
| `prop` | Desired proportionality function. Must be one of `"constant"`, `"loglin"`, `"ext_loglin"`, or `"logpoly"` | Defaults to `"constant"` |
| print | Information on non-linear minimisation interations | Defaults to `0` (not printing) |

It returns

| Variable | Description |
|:--------:|-------------|
| `qhat` | The estimate of parameters making up $q(a, a')$ |
| `deviance` | The deviance of the fitting of `qhat` |
| `aic` | The AIC of the fitting of `qhat` |
| `bic` | The BIC of the fitting of `qhat` |
| `bij` | The transmission matrix |
| `R0` |  The basic reproduction number |
| `lambda` | The force of infection |
| `R` | The reproduction number |
| `pi` | The seroprevalence |
| `inputs` | A list of the inputs used |
| `iterations` | The number of iterations used to fit the estimate(s) of the proportionality factor |

To fit the two models we are currently considering for Italy, run
```R
get_data(code = "IT")
FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
    muy = predict(demfit, type = "response"),
    N = sum(PS), Dur = 6 / 365, A = 0.5, Lmax = 70, 
    prop = "constant", startpar = 0.5)
FOI(age = sero$AGE, y = sero$indic, rij = contact_w,
    muy = predict(demfit, type = "response"),
    N = sum(PS), Dur = 6 / 365, A = 0.5, Lmax = 70, 
    prop = "loglin", startpar = c(0.5, 0.3))
```

## examine_convergence.R

This script is used to examine the performance of the model fitting. The functions are created such that they wrap around the functions of the previous two scripts. 

### Minimisation 
The function `check_conv` returns a table containing information regarding the minimisation performed in `FOI`. To obtain this information for the two models above, we use
```R
check_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
check_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "loglin", startpar = c(0.5, 0.3))
```
A function providing a plot of the output of `check_conv` is also included and is called `plot_conv` and takes the same arguments as above, i.e.
```R
plot_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
```
We have further included the function `rates_conv` which calculates the reproduction numbers for each of the values of the proportionality factor estimate obtained during the optimisation.
```R
rates_conv(code = "IT", Dur = 6 / 365, A = 0.5, Lmax = 70, prop = "constant", startpar = 0.5)
```
### Minimisation 
We generate starting parameters to examine the difference between the starting value and the estimated value. For this we have created a function `run_model` which wraps around `FOI` and returns starting values, estimated values, reproduction numbers, and an indicator for the country.

## results.R
This script contains a loop plotting the force of infection for each country. It currently takes the same time inputs for each country, i.e.:
```R
Dur <- 6 / 365
A <- 0.5
Lmax <- 70
```
as well as the same starting values, though this could be updated to consider different inputs for each country. The plot includes the sum of the absolute difference between the seroprevalence point and the values of the curve for the same age as a rough measure of performance.
