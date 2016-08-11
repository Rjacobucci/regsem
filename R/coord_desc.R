
coord_desc <- function(start,func,grad,hess,pars_pen,model,lambda,mats,block,max.iter,tol){
  count = 0
  ret <- list()
  max.iter = max.iter
  tol=tol

  # mats
 # mats <- extractMatrices(model)


  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 0
  new.pars <- matrix(NA,max.iter+1,length(start))
  new.pars[1,] <- start


  while(count < max.iter){
    count=count+1

    if(count < 50){
      alpha <- .5
    }else{
      alpha=.5
    }





      update.pars <- new.pars[count,]
      #a.pars <- update.pars[1:max(mats$A)]
      #s.pars <- update.pars[min(mats$S != 0):max(mats$S)]
    # gg <- grad(new.pars[count,])

      if(block == FALSE){
        for(j in 1:length(update.pars)){ # update A

          gg <- grad(update.pars)
          nn.par <- update.pars[j] - alpha*gg[j]

          if(any(j == pars_pen) & lambda > 0){
              update.pars[j] <- sign(nn.par)*max(abs(nn.par)-lambda,0)
          }else{
              update.pars[j] <- nn.par
          }
        }
      }else if(block==TRUE){
        # A
        gg <- grad(new.pars[count,])
        update.pars[1:max(mats$A)] <- update.pars[1:max(mats$A)] - alpha*gg[1:max(mats$A)]

        if(lambda > 0){
          for(j in pars_pen){
            update.pars[j] <- sign(update.pars[j])*max(abs(update.pars[j])-lambda,0)
          }
        }


        # S
        gg2 <- grad(update.pars)

        update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
                 update.pars[min(mats$S[mats$S !=0]):max(mats$S)] - alpha*gg2[min(mats$S[mats$S !=0]):max(mats$S)]

      }


    new.pars[count+1,] <- update.pars



    vals[count+1] = func(new.pars[count+1,])


    st.crit = try(abs(vals[count+1] - vals[count])<tol)

    st.crit2 <- all(abs(gg) < .01)
    dif <- abs(vals[count+1] - vals[count])

    #print(round(dif,5))
   # print(as.vector(round(gg,3)))
   # print(dif)
  #  print(round(gg,3))
  #  print(convergence)

    if(inherits(st.crit, "try-error")){
      convergence=99
    }else if(is.na(st.crit)==TRUE){
      convergence=99
    }else{
      if(st.crit==TRUE){
        convergence = 0
        #print(convergence)
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
