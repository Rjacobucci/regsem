fiml_calc4 <- function(ImpCov,F,mats2,type,lambda,model,pen_vec){
  A <- mats2$A
  S <- mats2$S
  A_est <- mats2$A_est
  S_est <- mats2$S_est
  F <- mats2$F
  I <- diag(nrow(A_est))

  nncol = which(colnames(A_est) == "1")
  m.pars = A[-nncol,"1"]
  A.pars = A[-nncol,-nncol]
  S.pars = S[-nncol,-nncol]

  F0<-F[rownames(F)!="1",colnames(F)!="1"]
  A0<-A_est[rownames(A_est)!="1",colnames(A_est)!="1"]
  I0<-diag(nrow(A0))
  S0<-S_est[rownames(S_est)!="1",colnames(S_est)!="1"]
  Imp.mean<-A_est[-which(colnames(A_est) == "1"),colnames(A_est)=="1"]
  Imp.cov<-F0 %*% solve(I0-A0) %*% S0 %*% t(solve(I0-A0)) %*% t(F0)
  mu = F0 %*% solve(I0 - A0) %*%Imp.mean


  lav.miss <-  model@SampleStats@missing[[1]]

  nobs <- model@Data@nobs[[1]]
  npatterns <- length(lav.miss)

  fit <- rep(NA,npatterns)

  samps <- rep(NA,npatterns)

  for(i in 1:npatterns){

    var.miss <- lav.miss[[i]]$var.idx  #missing pattern, T=exist, F=missing?

    X <-model@Data@X[[1]][model@Data@Mp[[1]]$case.idx[[i]],var.miss]

    nvar <- sum(var.miss)

    ImpCov<-Imp.cov[var.miss,var.miss]

    samps[i]<- lav.miss[[i]]$freq
    ###########################

    ll_part2 = matrix(NA,samps[i],1)

    X_minus_mu<-X-matrix(rep(mu[var.miss],each=samps[i]),ncol=nvar)

    for(j in 1:samps[i]){
      ll_part2[j]<-1/2*t(X_minus_mu[j,])%*%solve(ImpCov)%*%X_minus_mu[j,]
    }
    ##########################################

    ll_part2a = colSums(ll_part2)

    ll = - (samps[i] * nvar)/2 * log(2 * pi) - samps[i]/2 * log(det(ImpCov)) - ll_part2a

    fit[i] <- ll
  }

  #loglikelihood to fiml discrepancy  function (fit.fiml)
  #ref: SAS User Guide - estimation criteria - fiml
  fit.fiml = -2*sum(fit)/nobs

  #add penalty (ref:regsem::rcpp_fit_fun)
  add<-0

  if(type %in% c("lasso","alasso")){
    #print(pen_vec)
    add <- lambda * sum(abs(pen_vec))
  }else if(type=="ridge"){
    add <- lambda * sum(pen_vec)
  }#else if(type=="enet"){
    #elastic net
   # add <- lambda * sum(alpha*(Areg*Areg)  + (1- alpha)*abs(Areg))
  #}

  fit.sum <- fit.fiml + add

  ret<-list()
  ret$lik <- fit.sum
  ret$ImpCov <- Imp.cov
  ret$S2 <- S0
  ret$A2 <- A0
  ret$m <- Imp.mean
  ret$m.pars <- m.pars
  ret$A.pars <- A.pars
  ret$S.pars <- S.pars
  ret$F <- F0
  ret$mu <- mu
  ret
}
