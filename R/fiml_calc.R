
fiml_calc = function(ImpCov,mu.hat,h1,Areg,lambda,alpha,type,pen_vec,nvar,lav.miss){

# look at lav_objective.R from lavaan
  m = dim(ImpCov)[1]
  IntCol = which(colnames(Areg) == "1")
  IntCol2 = colnames(Areg[,1:nvar])
  fit = 0



  inv.chol <- function(S, logdet=FALSE) {
    cS <- chol(S)
    #if( inherits(cS, "try-error") ) {
    #    print(S)
    #    warning("lavaan WARNING: symmetric matrix is not positive symmetric!")
    #}
    S.inv <- chol2inv( cS )
    if(logdet) {
      diag.cS <- diag(cS)
      attr(S.inv, "logdet") <- sum(log(diag.cS*diag.cS))
    }
    S.inv
  }






  if(type=="none"){

    estimator.FIML <- function(Sigma.hat=NULL, Mu.hat=NULL, M=NULL, h1=NULL) {

      npatterns <- length(M)

      fx.p <- numeric(npatterns)
      w.p <- numeric(npatterns)

      # for each missing pattern, combine cases and compute raw loglikelihood
      for(p in 1:npatterns) {

        SX <- M[[p]][["SX"]]
        MX <- M[[p]][["MX"]]
        w.p[p] <- nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]

        # note: if a decent 'sweep operator' was available (in fortran)
        # we might win some time by 'updating' the inverse by sweeping
        # out the changed patterns... (but to get the logdet, we need
        # to do it one column at a time?)

        #cat("FIML: pattern ", p, "\n")
        #print(Sigma.hat[var.idx, var.idx])
        #print(Sigma.hat[var.idx,var.idx])

        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx], logdet=TRUE)
        #Sigma.inv <- chol2inv(chol(Sigma.hat[var.idx, var.idx]))
        Sigma.log.det <- attr(Sigma.inv, "logdet")
        Mu <- Mu.hat[var.idx]

        TT <- SX + tcrossprod(MX - Mu)
        trace <- sum(Sigma.inv * TT)

        fx.p[p] <- Sigma.log.det + trace
      }

      fx <- weighted.mean(fx.p, w=w.p)

      # ajust for h1
      if(!is.null(h1)) {
        fx <- fx - h1

        # no negative values
        if(fx < 0.0) fx <- 0.0
      }

      fx
    }

    fit = estimator.FIML(Sigma.hat=ImpCov,Mu.hat=mu.hat,M=lav.miss,h1=h1)*0.5


  }else{
    stop("Only type==none is currently supported")
  }


  fit

}
