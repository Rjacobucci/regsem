
coord_desc <- function(start,func,grad,hess,pars_pen,model,lambda,mats){
  count = 0
  ret <- list()
  max.iter = 200
  tol=1e-6
  type2=TRUE # experimental
  dif <- 0.5

  # mats
 # mats <- extractMatrices(model)


  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 100
  new.pars <- matrix(NA,max.iter+1,length(start))
  new.pars[1,] <- start


  while(count < max.iter){
    count=count+1

    if(count < 5){
      alpha <- .5
    }else if(count < 50 & count >= 5 ){
      alpha <- .5
    }else{
      alpha <- .1
    }
    alpha <- .5


    if(is.null(hess)==TRUE & type2==FALSE){
      if(count < 150){
        new.pars[count+1,] = new.pars[count,] - 0.5 * grad(new.pars[count,])
      }else if(count >= 150){
        new.pars[count+1,] = new.pars[count,] - 0.1 * grad(new.pars[count,])
      }else if(count >= 500){
        new.pars[count+1,] = new.pars[count,] - 0.000001 * grad(new.pars[count,])
      }
    }else if(type2 == TRUE){

      update.pars <- new.pars[count,]
      #a.pars <- update.pars[1:max(mats$A)]
      #s.pars <- update.pars[min(mats$S != 0):max(mats$S)]
    # gg <- grad(new.pars[count,])
      for(j in 1:length(update.pars)){ # update A
#print(lambda)
        gg <- grad(update.pars)
        nn.par <- update.pars[j] - alpha*gg[j]#*update.pars[1:max(mats$A)]
        if(any(j == pars_pen) & lambda > 0){
          update.pars[j] <- sign(nn.par)*max(abs(nn.par)-lambda,0)
        }else{
          update.pars[j] <- nn.par
        }
       # update.pars[j] <- sign(nn.par)*max(nn.par-lambda,0)
        #if(abs(nn.par) < lambda) update.pars[j] = 0

      }
    # print(update.pars)

   # for(j in 1:length(update.pars)){
    #  gg <- grad(update.pars)

     # if(abs(grad(new.pars[count,]))[j] < 0.0001){
      #  update.pars[j] <- new.pars[count,j]
     # }else{
    #    update.pars[j] <- new.pars[count,j] - alpha*gg[j]
     # }
   # }

    #  ind <- which(abs(gg) == max(abs(gg)))
    #  print(gg)
     # print(ind)
    #  update.pars[ind] <- update.pars[ind] - .5*gg[ind]

    # for(j in min(mats$S):max(mats$S)){
   #     update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <- update.pars[min(mats$S[mats$S !=0]):max(mats$S)] - .5*grad(update.pars)[min(mats$S[mats$S !=0]):max(mats$S)]#*update.pars[min(mats$S):max(mats$S)]
   #   }

      #  new.pars[count+1,] <- update.pars

      #  for(j in 1:length(update.pars)){
      #    if(abs(grad(new.pars[count,j])) < 0.00001){
      #      new.pars[count+1,j] <- new.pars[count,j]
      #    }
       # }

    new.pars[count+1,] <- update.pars
   # print(update.pars)
    #print(update.pars)

   #   new.pars[count+1,pars_pen] <- new.pars[count,pars_pen]

    #  for(j in 1:length(pars_pen)){
    #    S <- sign(new.pars[count+1,pars_pen[j]]) * max(new.pars[count+1,pars_pen[j]]-lambda,0)
    #    new.pars[count+1,pars_pen[j]] <- S
   #   }




      }else{
      hh = try(solve(hess(new.pars[count,])))

     # if(inherits(hh, "try-error")){
        delta = - hh %*% (grad(new.pars[count,]))
        new.pars[count+1,] = new.pars[count,] + delta
     # }else{
     #   if(count < 100){
     #     new.pars[count+1,] = new.pars[count,] - delta
     #   }else if(count >= 200){
     #     new.pars[count+1,] = new.pars[count,] - delta
     #   }else if(count >= 500){
     #     new.pars[count+1,] = new.pars[count,] - delta
     #   }
     # }
    }


    vals[count+1] = func(new.pars[count+1,])


    st.crit = try(abs(vals[count+1] - vals[count])<tol)

    st.crit2 <- all(abs(gg) < .01)
    dif <- abs(vals[count+1] - vals[count])

    #print(round(dif,5))
   # print(as.vector(round(gg,3)))
    #print(dif)
  #  print(round(gg,3))
  #  print(convergence)

    if(inherits(st.crit, "try-error")){
      convergence=99
    }else if(is.na(st.crit)==TRUE){
      convergence=99
    }else{
      if(st.crit==TRUE){
        convergence = 0
        print(convergence)
       # print(round(gg,3))
        break

      }
    }
  }
  ret$iterations <- count
  ret$value <- vals[count+1]
  ret$pars <- new.pars[count+1,] #+ rnorm(length(new.pars[count+1,]),0,0.00001)
  ret$convergence <- convergence
  ret

}
