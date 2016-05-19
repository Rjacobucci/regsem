
<!-- README.md is generated from README.Rmd. Please edit that file -->
regsem: Regularized Structural Equation Modeling
================================================

This package is based on the paper by Jacobucci, Grimm, & McArdle (2016), available at: <http://www.tandfonline.com/doi/full/10.1080/10705511.2016.1154793>. The package allows for using ridge or lasso penalties on specific parameters in general structural equation models. It is currently an extremely developmental version in both its implementation and theory behind the method.

First: Installing the package
-----------------------------

``` r
devtools::install_github("RJacobucci/regsem")
library(regsem)
```

It is best to install directly from Github as submissions to CRAN take place only once every couple of months.

As a simple example, variable selection can be done in a simple one factor model.

### Example 1: Variable selection in a One Factor CFA

``` r
library(regsem,quietly=T) # Depends on lavaan
#> This is lavaan 0.5-20
#> lavaan is BETA software! Please report any bugs.
HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
mod <- '
f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
'
# Recommended to specify meanstructure in lavaan
outt = cfa(mod,HS,meanstructure=TRUE)

fit1 <- regsem(outt,lambda=0.1,type="lasso")
summary(fit1)
#> $call
#> regsem(model = outt, lambda = 0.1, type = "lasso")
#> 
#> $estimates
#>   f -> x2 f -> x3 f -> x4 f -> x5 f -> x6 f -> x7 f -> x8 f -> x9 1 -> x1
#> 1   0.113   0.118   0.969   0.951   0.947   0.042   0.068    0.21       0
#>   1 -> x2 1 -> x3 1 -> x4 1 -> x5 1 -> x6 1 -> x7 1 -> x8 1 -> x9 x1 ~~ x1
#> 1       0       0       0       0       0       0       0       0    0.928
#>   x2 ~~ x2 x3 ~~ x3 x4 ~~ x4 x5 ~~ x5 x6 ~~ x6 x7 ~~ x7 x8 ~~ x8 x9 ~~ x9
#> 1    0.966    0.964    0.282    0.304    0.308    0.987     0.98    0.927
#>   f ~~ f
#> 1  0.703
#> 
#> $returnVals
#>       convergence df       fit   rmsea      BIC
#> rmsea           0 27 0.6538811 0.21275 7307.671
#> 
#> attr(,"class")
#> [1] "summary.regsem"
```

regsem uses the lavaan package for model specification and grabbing specific objects. Actual model computation is all done in regsem, with most of the heavy computation is done in C++. Different arguments used in any of lavaan's main functions (cfa,sem,lavaan,growth) will influence the model run in regsem. If you are familiar with RAM notation, you can check the matrix specification using extractMatrices(lavaan model name). This will give you detail to what the A and S matrices look like when used by regsem.

By default, regsem penalizes only the factor loadings (all) unless otherwise specified. If specific factor loadings, or other parameters are desired for regularization, use pars\_pen().

### Example 2: Penalizing specific parameters

``` r
# find parameter numbers
parTable(outt)
# or
extractMatrices(outt) # from regsem, can look at par #'s directly in A or S matrix

# only penalize first 3 loadings
fit2 <- regsem(outt,lambda=0.1,type="lasso",pars_pen=c(1:3))
summary(fit2)
```

Multiple Starting Values
------------------------

As currently implemented, regsem() encounters difficulty in convergence, particularly as lambda increases. Therefore, it is almost always recommended to use multi\_optim() with at least 10 different starting values.

### Example 3: Multiple starting values

``` r
?multi_optim
fit.mult <- multi_optim(outt,lambda=0.15,max.try=100,type="lasso",
                        max.iter=100000,tol=1e-4)
summary(fit.mult) 

# get fit of model -- using bootstrapping (naive)
fit_indices(fit.mult,CV="boot")
# to use CV=T, need to pre=partition the sample and provide a test covariance matrix
```

This isn't a perfect solution to the convergence difficulties and will be further addressed with different types of optimization in the future. With ridge penalties, it currently works best to not use a gradient for optimization (gradFun="none"). However, for lasso penalties, it is currently recommended to use the gradient specification from von Oertzen & Brick (2014) sparse RAM matrix derivatives paper (gradFun="ram").

### Example 4: Testing multiple lambda values

``` r
cv.out = cv_regsem(outt,type="ridge",gradFun="none",fit.ret2="boot",
                   n.lambda=40,mult.start=TRUE)

# see the parameter estimates for each value of lambda
cv.out[[1]]
# see rmsea and BIC for each value of lambda
cv.out[[2]]

# lowest

min.bic <- min(cv.out[[2]][,"BIC"])
loc <- which(cv.out[[2]][,"BIC"] ==  min.bic)
cv.out[[2]][1,]
```

A plotting function to visualize the changes in fit and parameter estimates across values of lambda will be added in the future.

### Things that regsem does not currently support

-   missing data -- don't use missing="fiml" in lavaan code
-   outputting standard errors or p-values
-   multiple group models
-   categorical indicators -- will be treated as continuous

### Important things to note

-   generally, ridge penalties converge faster than lasso
-   fit indices should be interpreted with a grain of salt
-   using lasso, the degrees of freedom reflect the number of parameters not equal to zero
-   using ridge, the degrees of freedom don't change (will be modified in future)

### Questions, Suggestions, Errors?

send an email to: <rcjacobuc@gmail.com>
