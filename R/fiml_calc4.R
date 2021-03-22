fiml_calc4 <- function(ImpCov,F,mats2,type,lambda,model,sat.lik,pen_vec){
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
  
  
  lav.miss <-  fits@SampleStats@missing[[1]]
  
  nobs <- model@Data@nobs[[1]]
  npatterns <- length(lav.miss)
  
  fit <- rep(NA,npatterns)
  # fit2 <- rep(NA,npatterns)
  # fit3 <- rep(NA,npatterns)
  samps <- rep(NA,npatterns)
  
  for(i in 1:npatterns){
    
    var.miss <- lav.miss[[i]]$var.idx  #missing pattern, T=exist, F=missing?
    
    #X <- lav.miss[[i]]$X
    #################################################
    #  does not exist
    # is this the raw data corresponding to the missing pattern (with missing var removed)?
    # if yes, may be changed to:
    X <-fits@Data@X[[1]][fits@Data@Mp[[1]]$case.idx[[i]],var.miss]
    #################################################
    
    #Scov <- lav.miss[[i]]$SX
    #mean <- lav.miss[[i]]$MX
    
    #F2 <- F[var.miss,c(var.miss,rep(T,ncol(F)-nrow(F)))]
    #A2 <- A_est[c(var.miss,rep(T,ncol(F)-nrow(F))),c(var.miss,rep(T,ncol(F)-nrow(F)))]
    #S2 <- S_est[c(var.miss,rep(T,ncol(F)-nrow(F))),c(var.miss,rep(T,ncol(F)-nrow(F)))]
    #############################################
    #rep(T,ncol(F)-nrow(F)) to include latent vars
    #intercepts not considered; var.miss only consider observed vars, however in F there may be a "1" for the intercept (regression for example), and since the n.col of F does not match the length of var.miss, 1 will be overwrite by the missing pattern of the 1st var
    #maybe change to:
    #if (sum(!var.miss)==0){ # no missing
    #  F2<-F
    #  A2<-A_est
    #  S2<-S_est
    #}else{
    #  which.miss<-which(!var.miss)
    #  F2<-F[-which.miss,-which.miss]
    #  A2<-A_est[-which.miss,-which.miss]
    #  S2<-S_est[-which.miss,-which.miss]
    #}
    
    #############################################
    
    #I = diag(nrow(A2))
    nvar <- sum(var.miss)
    
    #ImpCov = F2 %*% solve(I-A2) %*% S2 %*% t(solve(I-A2)) %*% t(F2)
    #ImpCov=lav.miss[[i]]$SY
    #ImpCov <- fits@implied$cov[[1]][var.miss,var.miss] #Correct! but need to use ram matrixinstead
    ImpCov<-Imp.cov[var.miss,var.miss]
    
    
    #samps[i] <- lav.miss[[i]]$nobs
    #####################no such info found:
    # in lav_SampleStats.R the returns are:
    #Yp[[p]] <- list(SY = SY, MY = MY, var.idx = Mp$pat[p,],freq = FREQ)
    
    #change to :
    samps[i]<- lav.miss[[i]]$freq
    ###########################
    
    
    
    # fit[i] <- (log(det(ImpCov)) + trace(Scov %*% solve(ImpCov)) - log(det(Scov)) - nvar)
    # fit2[i] <- (-nvar * log(2*pi) + log(det(ImpCov)) + t(mean -
    #                 rep(0,nvar)) %*% solve(ImpCov) %*% (mean-rep(0,nvar)))
    
    #   fit2[i] <- log(det(ImpCov)) + t(mean -rep(0,nvar)) %*% solve(ImpCov) %*% (mean-rep(0,nvar))
    
    
    
    
    ll_part2 = matrix(NA,samps[i],1)
    
    
    
    #for (j in 1:samps[i]){
    #ll_part2[j] = ((t(as.matrix(X[j,])) - rep(0,nvar)) %*% solve(ImpCov)  %*% (as.matrix(X[j,]) - rep(0,nvar)))*.5
    #}
    
    
    ##########################################
    #should be mean instead of 0, and dimension won't match if intercept is included
    #the above part from "for (j in 1:samps[i])" change to: (using matrix algebra)
    #X_mi_mu<-X-matrix(rep(colMeans(X),each=dim(X)[1]),ncol=dim(X)[2])
    #X_minus_mu<-cbind(X_mi_mu,rep(0,dim(ImpCov)[1]-dim(X)[2]))
    #X_minus_mu<-X-matrix(rep(fits@implied$mean[[1]][var.miss],each=samps[i]),ncol=nvar)#correct, but need to use RAM matrices
    X_minus_mu<-X-matrix(rep(mu[var.miss],each=samps[i]),ncol=nvar)
    #fits@implied$mean[[1]][var.miss]
    #X_minus_mu<-X-matrix(rep(lav.miss[[i]]$MY,each=dim(X)[1]),ncol=dim(X)[2])
    
    for(j in 1:samps[i]){
      ll_part2[j]<-1/2*t(X_minus_mu[j,])%*%solve(ImpCov)%*%X_minus_mu[j,]
    }
    ##########################################
    
    ll_part2a = colSums(ll_part2)
    
    ll = - (samps[i] * nvar)/2 * log(2 * pi) - samps[i]/2 * log(det(ImpCov)) - ll_part2a
    
    add <- 0
    
    
    if(type=="lasso"){
      print(pen_vec)
      add <- 2*lambda * sum(abs(pen_vec))
    }
    #print(add * samps[i]/nobs)
    #print(samps[i]/nobs)
    #print(nobs)
    #print(samps[i])
    #print(add * samps[i]/nobs)
    fit[i] <- ll - (add * samps[i]/nobs)
    
    #  TT <- Scov + tcrossprod(mean - rep(0,nvar))
    #  Sigma.inv <- log(det(ImpCov))
    #  trace2 <- sum(Sigma.inv * TT)
    
    #  fit3[i] <- (-nvar * log(2*pi)) + Sigma.inv + trace2
    
  }
  fit.sum = sum(fit)
  
  fit2 = -2 * fit.sum
  
  ret<-list()
  ret$lik <- fit2
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
