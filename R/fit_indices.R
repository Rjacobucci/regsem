#'
#'
#' Calculates the fit indices
#' @param model regsem model object.
#' @param CV cross-validation. Note that this requires splitting the dataset
#'           into a training and test set prior to running the model. The
#'           model should be run on the training set, with the test set
#'           held out and then passed to CovMat=.
#' @param CovMat If CV=T then test covariance matrix must be supplied. Note
#'        That this should be done before running the lavaan model and should
#'        not overlap with the data or covariance matrix used to run the model.
#' @param data supply the dataset?
#' @param n.obs Number of observations in the test set for CV.
#' @return fits Full set of fit indices
#' @keywords fit chisq ncp rmsea
#' @export
#' @examples
#' \dontrun{
#' fit_indices()
#' }


fit_indices =  function(model,CV=FALSE,CovMat=NULL,data=NULL,n.obs=NULL){

  res <- list()
  ret <- as.vector(rep(0,24))

  names(ret) <- c("Fmin","p","chisq","p.chisq","nfac","df","npar","N","baseline.chisq","baseline.df",
                  "logl","ncp","rmsea","rmsea.lower","rmsea.upper","rmsea.pval","srmr",
                  "CFI","TLI","BIC","AIC",
                  "CAIC","EBIC.5","EBIC.25")





  p=model$nvar
  ret["p"] = p
  nfac = model$nfac
  ret["nfac"] = nfac
  df = model$df
  ret["df"] = df
  npar = model$npar
  ret["npar"] = npar

  if(CV==FALSE){
    N = model$N
    ret["N"] = N
  }else{
    if(is.null(n.obs)){
      stop("For CV, need to provide sample size of test data")
    }
    N = n.obs
    ret["N"] = N
  }


  if(CV==FALSE){
    fit = model$fit
    res$Data_Type = "Train"
    SampCov = model$SampCov
    if(model$mean==TRUE){
      ImpCov = model$Imp_Cov1
      SampCov2 = model$SampCov2
    }else{
      ImpCov = model$Imp_Cov
      SampCov2 = SampCov
    }
  }else if(CV==TRUE){
    if(is.null(CovMat) ==TRUE){
      stop("Need to Provide Test CovMat")
    }
    res$Data_Type = "Test"

    #dat <- model$data
  # ids <-  sample(nrow(dat),nrow(dat)/2)
    SampCov <- SampCov2 <- CovMat


    if(model$mean==TRUE){
      ImpCov = model$Imp_Cov1
    }else{
      ImpCov = model$Imp_Cov
    }

   # SampCov=CovMat
    fit = 0.5*(log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov))  - p)
  }#else if(CV=="boot"){
   # fit.rep <- rep(NA,n.boot)
   # ImpCov = model$Imp_Cov
   # SampCov = model$SampCov
   # data <- model$data
    #data = model@Data@X[[1]]
   #   for(i in 1:n.boot){

   #     ids <- sample(nrow(data),nrow(data),replace=TRUE)
   #     new.dat <- data[ids,]
   #     SampCov.boot <- cov(new.dat)
   #     fit.rep[i] = 0.5*(log(det(ImpCov)) + trace(SampCov.boot %*% solve(ImpCov)) - log(det(SampCov.boot))  - p)

   #   }
   # fit = mean(fit.rep)
   # varFit = var(fit.rep)
 # }



  ret["Fmin"] = fit

 # if(CV=="boot"){
 #   ret["varFit"] = varFit
 # }else{
 #   ret["varFit"] = 0
 # }

  chisq = fit*N*2
  ret["chisq"] = chisq




  if(model$lav.model@Fit@converged == FALSE){
    baseline.chisq = NA
    baseline.df = NA
  }else{
    lav.fits <- fitmeasures(model$lav.model)
    baseline.chisq = lav.fits["baseline.chisq"]
    baseline.df = lav.fits["baseline.df"]
  }


  #baseline.chisq = model$baseline.chisq
  ret["baseline.chisq"] = baseline.chisq
  #baseline.df = model$baseline.df
  ret["baseline.df"] = baseline.df


#  c <- N*p/2 * log(2 * pi)



  if(model$lav.model@Fit@converged == FALSE){
    logl_sat = NA
  }else{
    logl_sat <- as.numeric(lav.fits["unrestricted.logl"])
  }
  #logl_sat = model$logl_sat# -c -(N/2) * log(det(SampCov)) - (N/2)*p
  logl = -N * (fit- logl_sat/N)
  ret["logl"] = logl




    d = function(chisq,df,N) max(0,(chisq -df)/(N-1))
    ncp = d(chisq,df,N)
    ret["ncp"] = ncp


    pp <- try(pchisq(chisq,df))
    if(inherits(pp, "try-error")) {
      ret["p.chisq"] <- NA
    }else{
      ret["p.chisq"] <- 1-pp
      }




    rmsea = function(ncp,df) sqrt(ncp/df)
    ret["rmsea"] = max(rmsea(ncp,df),0)


    lower.lambda <- function(lambda) {
      (pchisq(chisq, df=df, ncp=lambda) - 0.95)
    }

    lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=chisq)$root,
                    silent=TRUE)
    if(inherits(lambda.l, "try-error")) {
      ret["rmsea.lower"] = NA
    }else{
      ret["rmsea.lower"] = sqrt( lambda.l/(N*df))
      }


    upper.lambda <- function(lambda) {

      ppp <- try(pchisq(chisq, df=df, ncp=lambda))
      if(inherits(pp, "try-error")) {
        NA
      }else{
        ppp-0.05
      }
    }

    lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=100)$root,
                    silent=TRUE)
    if(inherits(lambda.u, "try-error")){
      ret["rmsea.upper"] = NA
    }else{
      ret["rmsea.upper"] = sqrt( lambda.u/(N*df))
    }




    pppp <- try(pchisq(chisq, df=df, ncp=(N*df*0.05^2)))
    if(inherits(pp, "try-error")) {
      ret["rmsea.pval"] <- NA
    }else{
      ret["rmsea.pval"] <- 1-pppp
    }


    #srmr
    imp = cov2cor(ImpCov);obs = cov2cor(SampCov2)
    lobs <-  obs[!lower.tri(obs)]
    limp <-  imp[!lower.tri(imp)]
    ret["srmr"] <- sqrt(mean((limp - lobs)^2))

    BIC<- function(logl,N,df,p){
      -2*(logl) + log(N)*npar
    }
    ret["BIC"]= BIC(logl,N,df,p)


    AIC<- function(logl,df,p){
      -2*(logl) + 2*npar
    }
    ret["AIC"] = AIC(logl,df,p)


    CAIC<- function(logl,df,N){
      -2*(logl) + log(N+1) * df
    }
    ret["CAIC"] = CAIC(logl,df,N)


    EBIC<- function(logl,df,N,p,nfac,delta){
      -2*(logl) + log(N) * df + 2*df * delta * log(p * nfac)
    }
    ret["EBIC.5"] = EBIC(logl,df,N,p,nfac,0.5)
    ret["EBIC.25"] = EBIC(logl,df,N,p,nfac,0.25)

    ncp.null = d(baseline.chisq,baseline.df,N)
    CFI <- function(ncp.null,ncp){
      (ncp.null - ncp) / ncp.null
    }
    ret["CFI"] = CFI(ncp.null,ncp)

    TLI <- function(baseline.chisq,baseline.df,chisq,df){
      (baseline.chisq/baseline.df - chisq/df) / (baseline.chisq/baseline.df - 1)
    }
    ret["TLI"] = TLI(baseline.chisq,baseline.df,chisq,df)

  res$fits <- round(ret,5)
    #ret
  res
}



