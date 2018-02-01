## ----results='hide',warning=FALSE,message=FALSE--------------------------
# install.packages("regsem") # lavaan is a dependency
# install.packages("semPlot")
library(semPlot) # for plotting the model
library(lavaan)
library(regsem)

## ------------------------------------------------------------------------
sim.mod <- "
f1 =~ 1*y1 + 1*y2 + 1*y3+ 1*y4 + 1*y5
f1 ~ 0*x1 + 0*x2 + 0*x3 + 0*x4 + 0*x5 + 0.2*x6 + 0.5*x7 + 0.8*x8
f1~~1*f1"
dat.sim = simulateData(sim.mod,sample.nobs=60,seed=12)

## ------------------------------------------------------------------------
run.mod <- "
f1 =~ NA*y1 + y2 + y3+ y4 + y5
f1 ~ c1*x1 + c2*x2 + c3*x3 + c4*x4 + c5*x5 + c6*x6 + c7*x7 + c8*x8
f1~~1*f1
"
lav.out <- sem(run.mod,dat.sim,fixed.x=FALSE)
#summary(lav.out)
parameterestimates(lav.out)[6:13,] # just look at regressions

## ----message=FALSE,warning=FALSE,fig.width=5,fig.height=5----------------
semPaths(lav.out)

## ----results='hide'------------------------------------------------------
reg.out <- cv_regsem(lav.out,n.lambda=50,type="lasso",jump=0.03,
                     pars_pen=c("c1","c2","c3","c4","c5","c6","c7","c8"))

## ----eval=FALSE----------------------------------------------------------
#  # not run
#  extractMatrices(lav.out)["A"]

## ------------------------------------------------------------------------
summary(reg.out)

## ----fig.width=5,fig.height=5--------------------------------------------

plot(reg.out)

## ------------------------------------------------------------------------
head(reg.out$fits,10)

## ------------------------------------------------------------------------
reg.out$final_pars[1:13] # don't show variances/covariances

