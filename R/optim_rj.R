
optim_rj <- function(start,func,grad,hess,pars_pen,model,lambda){
  count = 0
  ret <- list()
  max.iter = 200
  tol=1e-4
  type2=FALSE # experimental
  dif <- 0.5
  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 100
  new.pars <- matrix(NA,max.iter+1,length(start))
  new.pars[1,] <- start


  while(count < max.iter){
    count=count+1


    if(is.null(hess)==TRUE & type2==FALSE){
      if(count < 150){
        new.pars[count+1,] = new.pars[count,] - 0.5 * grad(new.pars[count,])
      }else if(count >= 150){
        new.pars[count+1,] = new.pars[count,] - 0.1 * grad(new.pars[count,])
      }else if(count >= 500){
        new.pars[count+1,] = new.pars[count,] - 0.000001 * grad(new.pars[count,])
      }
    }else if(type2 == TRUE){

      if(count < 150){
        new.pars[count+1,] = new.pars[count,] - 1 * grad(new.pars[count,])
      }else if(count >= 150){
        new.pars[count+1,] = new.pars[count,] - 1 * grad(new.pars[count,])
      }else if(count >= 500){
        new.pars[count+1,] = new.pars[count,] - 1 * grad(new.pars[count,])
      }

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
    dif <- abs(vals[count+1] - vals[count])
    print(dif)

    if(inherits(st.crit, "try-error")){
      convergence=99
    }else if(is.na(st.crit)==TRUE){
      convergence=99
    }else{
      if(st.crit==TRUE){
        convergence = 0
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
