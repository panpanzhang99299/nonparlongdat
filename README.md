Nonparametric Estimation Methods for Longitudinal Data Analysis
================

The $\mathtt{nonparlongdat}$ package intends to be a user-friendly
package implementing a nonparametric estimation method for missing
covariates in longitudinal data analysis.

## Features

The package provides functions to deal with

1.  time-invariant missing covariates with time-invariant complete
    auxiliary variables
2.  time-varying missing covariates with time-invariant complete
    auxiliary variables
3.  time-varying missing covariates with time-varying complete auxiliary
    variables

## Installation

You can install the $\mathtt{nonparlongdat}$ package from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("panpanzhang99299/nonparlongdat")
```

## Import Dataset

``` r
mock.df <- readr::read_csv(system.file("extdata", "mock_data.csv", package = "nonparlongdat"))
head(mock.df)
```

    ## # A tibble: 6 × 5
    ##      id  resp   covX  covZ   aux
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     1  14.0  0.847     0  2.95
    ## 2     2  16.5 NA         0  2.95
    ## 3     3  16.4  0.644     0  2.88
    ## 4     4  12.2 NA         0  3.01
    ## 5     5  15.5  0.744     0  3.05
    ## 6     6  14.3 NA         0  2.98

## Initial Values

There can be different ways of deriving the initial values for the
target parameters. We present an example based on an **available case
analysis** (ACA).

``` r
mock.aca.df <-
  mock.df %>% dplyr::filter(complete.cases(.)) %>% as.data.frame()

lme.aca <-
  nlme::lme(resp ~ covX + covZ, data = mock.aca.df, random = ~ 1 |
              id)
est.aca <-
  c(as.numeric(nlme::fixed.effects(lme.aca)),
    as.numeric(nlme::VarCorr(lme.aca)[, 2]))
```

## Preparation

``` r
n.sub <- length(unique(mock.df$id))
m.visit <- sum(mock.df$id == 1)
data.y <- mock.df$resp
data.x <- mock.df$covX
data.z <- mock.df$covZ
data.aux <- mock.df$aux
para.ini <- est.aca
```

## Implementation

``` r
est.prop <- nonparlongdat::nonpar_X_timeinv_S_timeinv(
  n.sub,
  m.visit,
  data.y,
  data.x,
  data.z,
  data.aux,
  para.ini,
  cov.cont = TRUE,
  aux.cont = TRUE
)

print(est.prop)
```

    ## $est
    ## [1] 11.8395443  3.1906061 -0.2922046  0.8285213  1.1653613
    ## 
    ## $se
    ## [1] 0.52932028 0.72602025 0.03968970 0.05151625 0.02836309
    ## 
    ## $converge
    ## [1] 0
