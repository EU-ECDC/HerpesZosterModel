# Baseline varicella force of infection

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

## plot_data.R
This script contains plots of data provided by load_data.R

## MCMC.R
This contains the function `FoI` which calculates the force of infection. It takes the following arguments

| Variable | Description | Default value |
|:--------:|-------------|---------------|
| `age` | Age of serological information | |
| `y` | Serological information given as indicator variable | |
| `rij` | Contact matrix | |
| `muy` | Mortality | |
| `N` | Population size at demographic equilibrium | |
| `D` | Duration | |
| `Lmax` | Life expectancy | |
| `A` | Duration of maternal immunity | |
| `propFac` | Desired proportionality function. Must be one of `"constant"`, `"loglin"`, `"ext_loglin"`, or `"logpoly"` | Defaults to `"constant"` |