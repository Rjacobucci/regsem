
coord_desc <- function(start,func,type,grad,hess,hessFun,pars_pen,model,lambda,mats,
                       block,max.iter,tol,full,solver,solver.maxit,alpha.inc,momentum,step,
                       step.ratio,diff_par,pen_vec,e_alpha){
  count = 0
  ret <- list()
  max.iter = max.iter
  tol=tol
  #solver=TRUE
  line.search=FALSE


  if(step.ratio == TRUE){
    alpha1 <- .01*step
    alpha2 <- step
  }else if(step.ratio == FALSE){
    alpha <- alpha1 <- alpha2 <- step
  }

  # mats
 # mats <- extractMatrices(model)


  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 0
  new.pars <- matrix(NA,max.iter+1,length(start))
  new.pars[1,] <- start
 # print(new.pars[1,])

  while(count < max.iter){
    count=count+1

    if(alpha.inc==FALSE){
      alpha <- step
    }else if(alpha.inc==TRUE){
      alpha <- 0.01 + 0.01*count
    }



    if(type=="diff_lasso"){
      pen_vec = new.pars[count,pars_pen]
      pen_diff = pen_vec - diff_par
    }else{
      pen_diff=0
    }

      update.pars <- new.pars[count,]
      #print(round(update.pars,4))
      #a.pars <- update.pars[1:max(mats$A)]
      #s.pars <- update.pars[min(mats$S != 0):max(mats$S)]
    # gg <- grad(new.pars[count,])

    if(hessFun=="none" & solver==FALSE){
      if(block == FALSE){
        for(j in 1:length(update.pars)){ # update A

          gg <- grad(update.pars)
          nn.par <- update.pars[j] - alpha*gg[j]

          if(any(j == pars_pen) & type=="lasso" & lambda > 0){
              update.pars[j] <- sign(nn.par)*max(abs(nn.par)-lambda,0)
          }else{
              update.pars[j] <- nn.par
          }
        }
      }else if(block==TRUE){

        if(full==TRUE & solver == FALSE){

          gg <- grad(new.pars[count,])


         # print(func(new.pars[count,]))
          #update.pars2 <- new.pars[count,]

          update.pars <- new.pars[count,] - alpha*gg



          if(type!="none" | type!="ridge" | type!="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha)
            }
          }else if(type=="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- update.pars[j] + sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,0)
            }
          }else if(type=="alasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha)
            }
          }



        }else if(full==TRUE & solver == TRUE){
          out <- nlminb(new.pars[count,],func,grad,control=list(iter.max=solver.maxit))
         # print(out$objective)
          update.pars <- out$par

          if(type!="none" | type!="ridge" | type!="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha)
            }
          }else if(type=="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- update.pars[j] + sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,0)
            }
          }else if(type=="alasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha)
            }
          }

        }else if(full==FALSE & line.search==FALSE){

          gg <- grad(new.pars[count,])

          if(type=="diff_lasso"){
            update.temp <- update.pars[pars_pen]
          }
          update.pars[1:max(mats$A)] <- update.pars[1:max(mats$A)] - alpha1*t(gg[1:max(mats$A)])


          if(type=="lasso" | type=="alasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha1,e_alpha)
            }
          }else if(type=="diff_lasso" & lambda > 0){
            for(j in 1:length(pen_diff)){
              ind <- pars_pen[j]
              print(ind)
            #  print(round(sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,0),3))
             # print(ind); print(sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,0))
              #update.pars[j] <- sign(update.pars[j])*max(abs(update.pars[j])-alpha*lambda,0)
             # print(sign(pen_diff[j]));print(max(abs(pen_diff[j])-alpha*lambda,0))
           #   print(sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,abs(pen_diff[j])))

             # print(ind)
             # update.pars[ind] <- update.temp[j] - sign(pen_diff[j])*max(abs(pen_diff[j]) - alpha1*lambda,0)
              pen  <- pen_diff[j] - sign(pen_diff[j])*max(abs(pen_diff[j]) - alpha*lambda,0)
              update.pars[ind] <- update.temp[j] - pen
             # print(round(update.pars[ind],3));print(round(pen_diff[j],3))
             # update.pars[ind] <- update.temp[j] - pen_diff[j]

             # print(round(update.pars[ind] - pen_diff[j],3));print(round(update.pars[ind],3))
             # print(round(update.pars[ind],3))

            }
          }#else if(type=="alasso" & lambda > 0){
           # for(j in pars_pen){
           #   update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha1,e_alpha)
          #  }
          #}

          ####################################################
          ####################################################

          # one problem with convergence seems to be going to far with A
          # parameters. This is why it seems to work well to set alpha=1
          # for S, but use small step sizes for A
         # print(round(update.pars,3))

          # S
          gg2 <- grad(update.pars)
         # print(round(rbind(t(gg),t(gg2))),4)

          update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
            update.pars[min(mats$S[mats$S !=0]):max(mats$S)] - alpha2*gg2[min(mats$S[mats$S !=0]):max(mats$S)]




        }else if(full == FALSE & line.search==TRUE){
          gg <- grad(new.pars[count,])


         # print(new.pars[count,])
         # print(round(gg,3))
        #  print(func(new.pars[count,]-.001*gg))
          delta1 <- function(step){
            func(new.pars[count,]-step*gg)
          }

        #  s <- try(uniroot(f=delta1, c(0,1),f.lower=0),silent=TRUE)

        #  if(inherits(s, "try-error")) {
        #    s <- 0.01
        #  }else{
        #    s <- s$root
        #  }

          s <- 0.1


          update.pars[1:max(mats$A)] <- update.pars[1:max(mats$A)] - s*t(gg[1:max(mats$A)])


          if(type=="lasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- sign(update.pars[j])*max(abs(update.pars[j])-s*lambda,0)
            }
          }





          # S
          gg2 <- grad(update.pars)
          #print(rbind(t(gg),t(gg2)))

          delta2 <- function(step){
            func(new.pars[count,]-step*gg2)
          }


          s <- try(uniroot(f=delta2, c(0,1),f.lower=0),silent=TRUE)

        if(inherits(s, "try-error")) {
            s <- 0.01
        }else{
            s <- s$root
        }


          update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
            update.pars[min(mats$S[mats$S !=0]):max(mats$S)] - s*gg2[min(mats$S[mats$S !=0]):max(mats$S)]
        }
      }
    }else if(hessFun!="none" & solver==FALSE){
      #alpha <- .1 + .01*count
     # alpha <- 1
      #print(new.pars[count,])

      # A
      gg <- grad(new.pars[count,])
      hh <- hess(new.pars[count,])
      #print(round(solve(hh)%*%gg,3))
     update.pars[1:max(mats$A)] <- update.pars[1:max(mats$A)] - alpha1*solve(hh[1:max(mats$A),1:max(mats$A)]) %*% gg[1:max(mats$A)]
      #update.pars <- update.pars - alpha*(solve(hh) %*% gg)


      if(lambda > 0 & type=="lasso"){
        for(j in pars_pen){
          update.pars[j] <- sign(update.pars[j])*max(abs(update.pars[j])-alpha1*lambda,0)
        }
      }

      # S
      gg2 <- grad(update.pars)
      hh2 <- hess(update.pars)
      minS <- min(mats$S[mats$S !=0])
      maxS <- max(mats$S)
      print(eigen(hh2)$values)

      update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
        update.pars[minS:maxS] - alpha2*(solve(hh2[minS:maxS,minS:maxS])%*% gg2[minS:maxS])


     # print(min(mats$S[mats$S !=0]))
     # print(max(mats$S))


    }else if(solver==TRUE){


      out <- nlminb(new.pars[count,],func,grad,control=list(eval.max=1))
      #print(out$objective)
      update.pars <- out$par

      if(type!="none" | type!="ridge" | type!="diff_lasso" & lambda > 0){
        for(j in pars_pen){
          update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha)
        }
      }else if(type=="diff_lasso" & lambda > 0){
        for(j in pars_pen){
          #print(update.pars[j])
          update.pars[j] <- update.pars[j] + sign(pen_diff[j])*max(abs(pen_diff[j])-alpha*lambda,0)
        }
      }else if(type=="alasso" & lambda > 0){
        for(j in pars_pen){
          #print(update.pars[j])
          update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha)
        }
      }

      # S
      out <- nlminb(update.pars,func,control=list(eval.max=solver.maxit))
      pp.pars <- out$par

      update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <- pp.pars[min(mats$S[mats$S !=0]):max(mats$S)]


    }

    if(momentum==FALSE){
      new.pars[count+1,] <- update.pars
    }else if(momentum==TRUE){
      new.pars[count+1,] <- update.pars + (count/(count+3))*(update.pars-new.pars[count,])
    }



    if(type != "diff_lasso") pen_diff <- 0
    vals[count+1] = func(new.pars[count+1,])
  #  print(round(vals[count+1],3))

    st.crit = try(abs(vals[count+1] - vals[count])<tol)
  #  st.crit2 <- all(abs(gg) < .01)
    dif <- abs(vals[count+1] - vals[count])

   # print(rbind(update.pars,t(gg2)))
    #print(round(dif,5))
   # print(as.vector(round(gg,3)))
   # print(dif)
  #  print(round(gg,3))
  #  print(convergence)

    if(inherits(st.crit, "try-error")){
      convergence=99
    }else if(is.na(st.crit)==TRUE){
      convergence=99
    }else if(any(new.pars[count+1,] > 500)){
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
