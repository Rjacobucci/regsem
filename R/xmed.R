#'
#'
#' Function to performed exploratory mediation with continuous and categorical variables
#'
#' @param data Name of the dataset
#' @param iv Name (or vector of names) of independent variable(s)
#' @param mediators Name of mediators
#' @param dv Name of dependent variable
#' @param covariates Name of covariates to be included in model.
#' @param type What type of penalty. Options include lasso, ridge, and enet.
#' @param nfolds Number of cross-validation folds.
#' @param show.lambda Displays lambda values in output
#' @param epsilon Threshold for determining whether effect is 0 or not.
#' @param seed Set seed to control CV results
#' @return Coefficients from best fitting model
#' @export
#' @examples
#' \donttest{
#'# example
#'library(ISLR)
#'College1 = College[which(College$Private=="Yes"),]
#'Data = data.frame(scale(College1[c("Grad.Rate","Accept","Outstate","Room.Board","Books","Expend")]))
#'Data$Grad.Rate <- ifelse(Data$Grad.Rate > 0,1,0)
#'Data$Grad.Rate <- as.factor(Data$Grad.Rate)
#'#lavaan model with all mediators
#'model1 <-
#'  ' # direct effect (c_prime)
#'Grad.Rate ~ c_prime*Accept
#'# mediators
#'Outstate ~ a1*Accept
#'Room.Board ~ a2*Accept
#'Books ~ a3*Accept
#'Expend ~ a6*Accept
#'Grad.Rate ~ b1*Outstate + b2*Room.Board + b3*Books + b6*Expend
#'# indirect effects (a*b)
#'a1b1 := a1*b1
#'a2b2 := a2*b2
#'a3b3 := a3*b3
#'a6b6 := a6*b6
#'# total effect (c)
#'c := c_prime + (a1*b1) + (a2*b2) + (a3*b3) + (a6*b6)
#''
#'#p-value approach using delta method standard errors
#'fit.delta = sem(model1,data=Data,fixed.x=TRUE,ordered="Grad.Rate")
#'summary(fit.delta)
#'
#'#xmed()
#'
#'iv <- "Accept"
#'dv <- "Grad.Rate"
#'mediators <- c("Outstate","Room.Board","Books","Expend")
#'
#'out <- xmed(Data,iv,mediators,dv)
#'out
#'}


xmed = function (data, iv, mediators, dv, covariates = NULL, type = "lasso",
                 nfolds = 10, show.lambda = F, epsilon = 0.001, seed = NULL)
{
  Data <- data
  if (type == "lasso") {
    alpha = 1
  }
  else if (type == "ridge") {
    alpha = 0
  }
  else if (type == "enet") {
    alpha = 0.5
  }
  var.check = function(data) {
    data = as.data.frame(data)
    num.response.options = flag = integer(ncol(data))
    for (i in 1:ncol(data)) {
      num.response.options[i] = flag[i] = NA
      num.response.options[i] = length(unique(data[, i]))
      if (is.factor(data[, i]) & num.response.options[i] >
          2) {
        flag[i] = 2
      }
      else if (num.response.options[i] == 2) {
        flag[i] = 1
      }
      else if (num.response.options[i] != 2) {
        flag[i] = 0
      }
    }
    return(flag)
  }
  check.out <- var.check(Data[, c(iv, mediators, dv, covariates)])
  if (any(check.out == 2)) {
    stop("Factor variables with > 2 response options need to be recoded as integer or numeric variables")
  }
  data.proc <- caret::preProcess(Data[, c(iv, mediators, dv, covariates)])
  data2 <- predict(data.proc, Data[, c(iv, mediators, dv, covariates)])
  iv.mat <- as.matrix(data2[, iv])
  mediators.mat <- as.matrix(data2[, mediators])
  mediv.mat <- as.matrix(data2[, c(mediators,iv)])
  dv.mat <- as.matrix(data2[, dv])
  if (sum(is.na(iv.mat)) > 0 | sum(is.na(mediators.mat)) >
      0 | sum(is.na(dv.mat)) > 0) {
    stop("Missing values are not allowed")
  }
  if (var.check(dv.mat) == 0) {
    dv.class = "gaussian"
  }
  else if (var.check(dv.mat) == 1) {
    dv.class = "binomial"
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #for 1 IV
  if((length(iv)==1) & (is.null(covariates))){
    b.cv.lasso = glmnet::cv.glmnet(mediv.mat,
                                   dv.mat, alpha = alpha,
                                   family = dv.class, standardize = FALSE,
                                   nfolds = nfolds,
                                   penalty.factor = c(rep(1,ncol(mediv.mat) - 1),0))
    b.coefs = coef(b.cv.lasso,
                   s = b.cv.lasso$lambda.min)[-c(1,ncol(mediv.mat)+1),1]
    direct = coef(b.cv.lasso,
                  s = b.cv.lasso$lambda.min)[(ncol(mediv.mat)+1),1]
    a.cv.lasso = a.fit.lasso = vector("list", ncol(mediators.mat))
    a.lambda = numeric(ncol(mediators.mat))
    for (i in 1:ncol(mediators.mat)) {
      if (var.check(mediators.mat[, i]) == 0) {
        med.class = "gaussian"
      }
      else if (var.check(mediators.mat[, i]) == 1) {
        med.class = "binomial"
      }
      if (!is.null(seed)) {
        set.seed(seed)
      }
      a.cv.lasso[[i]] = glmnet::cv.glmnet(
        as.matrix(cbind(rnorm(nrow(data),1, 1e-04), iv.mat)),
        mediators.mat[, i], alpha = alpha, family = med.class,
        standardize = FALSE, nfolds = nfolds,
        intercept = F, penalty.factor = c(0, 1))
      a.lambda[i] = a.cv.lasso[[i]]$lambda.min
    }
    a.coefs = numeric(length(b.coefs))
    for (i in 1:length(a.coefs)) {
      if (!is.null(a.cv.lasso[[i]])) {
        a.coefs[i] = coef(a.cv.lasso[[i]],
                          s = a.cv.lasso[[i]]$lambda.min)[-1,1][2]
      }
    }
    names(a.coefs) = mediators
    res <- list()
    res$a.coefs <- a.coefs
    res$b.coefs <- b.coefs
    indirect = a.coefs * b.coefs
    selected <- names(indirect[abs(indirect) > epsilon])
    indirect <- as.data.frame(indirect)
    indirect = round(indirect, 4)
    indirect[abs(indirect) < epsilon] = 0
    indirect <- t(indirect)
    indirect[abs(indirect) >= epsilon] =
      as.numeric(indirect[abs(indirect) >= epsilon])
    if(show.lambda==T){
      res$a.lambda = a.lambda
      res$b.lambda = b.cv.lasso$lambda.min
    }
    res$selected = selected
    res$indirect <- indirect
  }
  #for multiple IVs
  else if((length(iv)>1) | !is.null(covariates)){
    medivcov.mat <- as.matrix(data2[, c(mediators,iv,covariates)])
    ivcov.mat <- as.matrix(data2[, c(iv,covariates)])
    b.cv.lasso = glmnet::cv.glmnet(medivcov.mat, dv.mat, alpha = alpha,
                                   family = dv.class, standardize = F, nfolds = nfolds, penalty.factor =
                                     c(rep(1,length(mediators)), rep(0,(length(iv)+length(covariates)))))
    b.coefs = coef(b.cv.lasso,
                   s = b.cv.lasso$lambda.min)[(2:(length(mediators)+1)),1]
    direct = coef(b.cv.lasso, s = b.cv.lasso$lambda.min)[
      (length(mediators)+2):(ncol(mediv.mat)+1),1]

    a.cv.lasso = a.fit.lasso = vector("list", ncol(mediators.mat))
    a.lambda = numeric(ncol(mediators.mat))
    for (i in 1:ncol(mediators.mat)) {
      if (var.check(mediators.mat[, i]) == 0) {
        med.class = "gaussian"
      }
      else if (var.check(mediators.mat[, i]) == 1) {
        med.class = "binomial"
      }
      if (!is.null(seed)) {
        set.seed(seed)
      }
      a.cv.lasso[[i]] = glmnet::cv.glmnet(ivcov.mat,
                                          mediators.mat[, i], alpha = alpha, family = med.class,
                                          standardize = FALSE, nfolds = nfolds)
      a.lambda[i] = a.cv.lasso[[i]]$lambda.min
    }
    a.coefs = matrix(nrow=length(iv),ncol=length(b.coefs))
    colnames(a.coefs) = mediators
    rownames(a.coefs) = iv
    for (i in 1:ncol(a.coefs)) {
      if (!is.null(a.cv.lasso[[i]])) {
        a.coefs[,i] = coef(a.cv.lasso[[i]],s=a.cv.lasso[[i]]$lambda.min)[
          (2:(length(iv)+1)),1]
      }
    }
    if(!is.null(covariates)){
      dvcov.coefs = coef(b.cv.lasso, s = b.cv.lasso$lambda.min)[
        ((ncol(mediv.mat)+2):(ncol(medivcov.mat)+1)),1]
      medcov.coefs = matrix(nrow=length(covariates),ncol=length(mediators))
      colnames(medcov.coefs) = mediators
      rownames(medcov.coefs) = names(dvcov.coefs) = covariates
      for(i in 1:ncol(medcov.coefs)){
        medcov.coefs[,i] = coef(a.cv.lasso[[i]],s=a.cv.lasso[[i]]$lambda.min)[
          ((ncol(iv.mat)+2):(ncol(ivcov.mat)+1)),1]
      }
    }
    res <- vector("list", length(iv))
    names(res) = iv
    for(i in 1:length(iv)){
      res[[i]]$a.coefs <- a.coefs[i,]
      res[[i]]$b.coefs <- b.coefs
      indirect = a.coefs[i,] * b.coefs
      selected <- names(indirect[abs(indirect) > epsilon])
      indirect <- as.data.frame(indirect)
      indirect = round(indirect, 4)
      indirect[abs(indirect) < epsilon] = 0
      indirect <- t(indirect)
      indirect[abs(indirect) >= epsilon] =
        as.numeric(indirect[abs(indirect) >= epsilon])
      if(show.lambda==T){
        res[[i]]$a.lambda = a.lambda
        res[[i]]$b.lambda = b.cv.lasso$lambda.min
      }
      res[[i]]$selected = selected
      res[[i]]$indirect <- indirect
    }
  }
  res$direct = direct
  if(!is.null(covariates)){
    res$medcov.coefs = medcov.coefs
    res$dvcov.coefs = dvcov.coefs
  }
  res$call <- match.call()
  class(res) <- "xmed"
  return(res)
}
