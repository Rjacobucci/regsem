
fiml_calc2 = function(ImpCov,F,mats2,type,lambda,model,sat.lik,pen_vec){

# look at lav_objective.R from lavaan
  m = dim(ImpCov)[1]
 # IntCol = which(colnames(Areg) == "1")
 # IntCol2 = colnames(Areg[,1:nvar])
 # fit = 0
  lav.miss <-  model@SampleStats@missing[[1]]
  A_est <- mats2$A_est22
  S_est <- mats2$S_est22


#  if(type=="none"){
    nobs <- model@Data@nobs[[1]]
    npatterns <- length(lav.miss)
    lav.miss <- model@SampleStats@missing[[1]]

    fit <- rep(NA,npatterns)
   # fit2 <- rep(NA,npatterns)
   # fit3 <- rep(NA,npatterns)
    samps <- rep(NA,npatterns)
    for(i in 1:npatterns){

      var.miss <- lav.miss[[i]]$var.idx

      X <- lav.miss[[i]]$X

      Scov <- lav.miss[[i]]$SX
      #mean <- lav.miss[[i]]$MX

      F2 <- F[var.miss,c(var.miss,rep(T,ncol(F)-nrow(F)))]
      A2 <- A_est[c(var.miss,rep(T,ncol(F)-nrow(F))),c(var.miss,rep(T,ncol(F)-nrow(F)))]
      S2 <- S_est[c(var.miss,rep(T,ncol(F)-nrow(F))),c(var.miss,rep(T,ncol(F)-nrow(F)))]

      I = diag(nrow(A2))
      nvar <- sum(var.miss)

      ImpCov = F2 %*% solve(I-A2) %*% S2 %*% t(solve(I-A2)) %*% t(F2)
      #ImpCov <- outt@implied$cov[[1]][var.miss,var.miss]

      samps[i] <- lav.miss[[i]]$nobs
     # fit[i] <- (log(det(ImpCov)) + trace(Scov %*% solve(ImpCov)) - log(det(Scov)) - nvar)
      # fit2[i] <- (-nvar * log(2*pi) + log(det(ImpCov)) + t(mean -
      #                 rep(0,nvar)) %*% solve(ImpCov) %*% (mean-rep(0,nvar)))

   #   fit2[i] <- log(det(ImpCov)) + t(mean -rep(0,nvar)) %*% solve(ImpCov) %*% (mean-rep(0,nvar))




      ll_part2 = matrix(NA,samps[i],1)

      for (j in 1:samps[i]){
        ll_part2[j] = ((t(as.matrix(X[j,])) - rep(0,nvar)) %*% solve(ImpCov)  %*% (as.matrix(X[j,]) - rep(0,nvar)))*.5


      }

      ll_part2a = colSums(ll_part2)

      ll = - (samps[i] * nvar)/2 * log(2 * pi) - samps[i]/2 * log(det(ImpCov)) - ll_part2a

      add <- 0


      if(type=="lasso"){
        #print(pen_vec)
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


 # }else{
 #   stop("Only type==none is currently supported")
 # }


  fit = sum(fit)
 # print(fit)
  fit2 = -2 * fit
  diff <- fit2 - (-2*sat.lik)
 # print(diff)
 # fit.ret <- diff/301
 # print(fit.ret)
 # fit.ret
 # print(diff)
  diff
}
