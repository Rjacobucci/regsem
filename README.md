
<!-- README.md is generated from README.Rmd. Please edit that file -->
regsem: Regularized Structural Equation Modeling
================================================

This package is based on the paper by Jacobucci, Grimm, & McArdle (2016) in the Structural Equation Modeling journal.

Installing the package
----------------------

``` r
devtools::install_github("RJacobucci/regsem")
library(regsem)
```

It is best to install directly from Github as submissions to CRAN take place only once every couple of months.

### Things that regsem does not currently support

-   missing data -- don't use missing="fiml" in lavaan code
-   outputting standard errors or p-values
-   multiple group models
-   categorical indicators -- will be treated as continuous
