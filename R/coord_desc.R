
coord_desc <- function(start,func,type,grad,hess,hessFun,pars_pen,model,lambda,mats,
                       block,max.iter,tol,full,solver,solver.maxit,alpha.inc,step,
                       step.ratio,diff_par,pen_vec,e_alpha,gamma,momentum,par.lim,quasi){
  count = 0
  ret <- list()
  max.iter = max.iter
  tol=tol
  #solver=TRUE
  phi_func <- rep(1,max.iter)
  phi_grad <- rep(1,max.iter)


  #if(type=="enet"){
 #   step=step*2
 # }
  line.search=FALSE

  if(step.ratio == TRUE){
    alpha1 <- .01*step
    alpha2 <- step
  }else if(step.ratio == FALSE){
    alpha <- alpha1 <- alpha2 <- step
  }

  # mats
 # mats <- extractMatrices(model)

  alpha.vec <- rep(NA,max.iter)
  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 0
  new.pars <- matrix(NA,max.iter+1,length(start))
  grad.vec <- matrix(NA,max.iter+1,length(start))
 # B <- rep(NA,max.iter+1)
 # B <- matrix(NA,2,2)
  new.pars[1,] <- start
 # print(new.pars[1,])

  while(count < max.iter){
    count=count+1
    if(alpha.inc==FALSE){
      alpha <- step
    }else if(alpha.inc==TRUE){
      alpha <- 0.01 + 0.01*count
    }else if(alpha.inc=="dec"){
      alpha <- step - 0.01*count
    }



    if(type=="diff_lasso"){
      pen_vec = new.pars[count,pars_pen]
      pen_diff = pen_vec - diff_par
    }else{
      pen_diff=0
    }

   # print(pen_diff)

      update.pars <- new.pars[count,]
      #print(round(update.pars,4))
      #a.pars <- update.pars[1:max(mats$A)]
      #s.pars <- update.pars[min(mats$S != 0):max(mats$S)]
    # gg <- grad(new.pars[count,])

    if(solver==FALSE){
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

        if(full==TRUE & quasi == FALSE){

          #print(round(new.pars[count,],3))
          gg <- try(grad(new.pars[count,]),silent=TRUE)
          #print(round(gg,2))
          if(inherits(gg, "try-error")) {
            gg <- rnorm(length(new.pars[count,]),0,.01)
          }else{
            gg <- gg
          }


          delta <- function(alpha){

              update.pars <- new.pars[count,] - alpha*gg

            # print(round(t(alpha*gg),3))
             if(type == "ridge" | type=="none"){
              update.pars <- update.pars
             }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
               for(j in pars_pen){
                  update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
                }
             }else if(type=="diff_lasso" & lambda > 0){
                cc=0
                for(j in pars_pen){
                  cc <- cc + 1
                  update.pars[j] <- update.pars[j] -
                   pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
               }
             }else if(type=="alasso" & lambda > 0){
                for(j in pars_pen){
                 #print(update.pars[j])
                 update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
               }
             }

              func(update.pars)
          } #end delta


         c=.5
          p=update.pars-new.pars[count,]#cbind(rep(0.5,length(new.pars[count,])))

          alpha=1
          if(count==1){
            alpha =step =1
            #
          }else{
            while(delta(alpha) > func(update.pars)+c*alpha*(t(gg)%*%p)){
              alpha = 0.5*alpha
            }
          }


          # works well with momentum -- delta{only contains func(update.pars)}
         #alpha <- optimize(delta,interval=c(0,step),maximum=FALSE)$minimum





           update.pars <- new.pars[count,] - alpha*gg

           # print(round(t(alpha*gg),3))
           if(type == "ridge" | type=="none"){
             update.pars <- update.pars
           }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
             for(j in pars_pen){
               update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
             }
           }else if(type=="diff_lasso" & lambda > 0){
             cc=0
             for(j in pars_pen){
               cc <- cc + 1
               update.pars[j] <- update.pars[j] -
                 pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
             }
           }else if(type=="alasso" & lambda > 0){
             for(j in pars_pen){
               #print(update.pars[j])
               update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
             }
           }




        }else if(full==TRUE & quasi == TRUE){

          #out <- nlminb(new.pars[count,],func,control=list(iter.max=1,eval.max=1,step.min=alpha,step.max=alpha)) # iter.max=1,eval.max=1
         # out <- nlminb(new.pars[count,],func,control=list(iter.max=1,eval.max=1,step.min=alpha,step.max=alpha))
          #out <- optim(new.pars[count,],func,method="BFGS")#,control=list(maxit=1,ndeps=rep(alpha,length(new.pars[count,]))))
          grad.vec[count,] <- grad(new.pars[count,])



          # if not using hessian for first iteration, best to take small step sizes (0.1)
          # step should be less than 1


          if(count == 1){
            s = cbind(new.pars[count,])
            y = cbind(grad.vec[count,])
            #alpha = 1 # always use as first step length
            alpha.vec[count] <- s1 <- alpha<- 1

            if(hessFun != "none"){
              H <- solve(hess(new.pars[count,]))
            }else{
              H <-  diag(length(new.pars[count,]))#as.numeric((t(y)%*%s)/t(y)%*%y)
            }
            dir <- -H %*% grad.vec[count,]
          }else{
            s = cbind(new.pars[count,] - new.pars[count-1,])
            y =  cbind(grad.vec[count,] - grad.vec[count-1,])

            #http://terminus.sdsu.edu/SDSU/Math693a/Lectures/18/lecture.pdf
            p <- as.numeric(1/(t(y)%*%s))
            H <- (diag(length(new.pars[count,])) - p*s%*%t(y))%*%H%*%(diag(length(new.pars[count,]))-p*y%*%t(s)) + p*s%*%t(s)
            dir <- -H %*% grad.vec[count,]
          }






          delta <- function(alpha){

            update.pars <- new.pars[count,] + alpha*dir

            # print(round(t(alpha*gg),3))
            if(type == "ridge" | type=="none"){
              update.pars <- update.pars
            }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
              for(j in pars_pen){
                update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
              }
            }else if(type=="diff_lasso" & lambda > 0){
              cc=0
              for(j in pars_pen){
                cc <- cc + 1
                update.pars[j] <- update.pars[j] -
                  pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
              }
            }else if(type=="alasso" & lambda > 0){
              for(j in pars_pen){
                #print(update.pars[j])
                update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
              }
            }

            func(update.pars)
          } #end delta


          # s1 <- try(uniroot(f=delta,interval=c(-1,1),f.lower=0,f.upper=100),silent=TRUE)
          #  print(s1)
          # if(inherits(s1, "try-error")) {
          #  alpha <- 0.01
          # }else{
          #   alpha <- s1$root
          # }

          if(count==1){
            alpha <- 1
          }else{
            alpha <- optimize(delta,interval=c(0.01,step),maximum=FALSE)$minimum
          }

          update.pars <- new.pars[count,] + alpha*dir

          # print(round(t(alpha*gg),3))
          if(type == "ridge" | type=="none"){
            update.pars <- update.pars
          }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
            }
          }else if(type=="diff_lasso" & lambda > 0){
            cc=0
            for(j in pars_pen){
              cc <- cc + 1
              update.pars[j] <- update.pars[j] -
                pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
            }
          }else if(type=="alasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
            }
          }



        }else if(full==FALSE & line.search==FALSE){

          gg <- grad(new.pars[count,])
          #print(round(gg,2))

          update.pars[1:max(mats$A)] <- update.pars[1:max(mats$A)] - alpha1*t(gg[1:max(mats$A)])


          if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
            for(j in pars_pen){
              update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha1,e_alpha,gamma)
            }
          }else if(type=="diff_lasso" & lambda > 0){
            cc=0
            for(j in pars_pen){
              cc <- cc + 1
              pen_vec = update.pars[pars_pen]
              pen_diff = pen_vec - diff_par

              print(pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma))
              update.pars[j] <- update.pars[j] -
                  pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
            }
          }else if(type=="alasso" & lambda > 0){
            for(j in pars_pen){
              #print(update.pars[j])
              update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha1,e_alpha,gamma)
            }
          }

          # S
          gg2 <- grad(update.pars)
          #print(round(rbind(t(gg),t(gg2))),4)

          update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
            update.pars[min(mats$S[mats$S !=0]):max(mats$S)] - alpha2*gg2[min(mats$S[mats$S !=0]):max(mats$S)]

          #print(round(update.pars,3))


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
     update.pars <- update.pars - alpha1*MASS::ginv(hh) %*% gg
      #update.pars <- update.pars - alpha*(solve(hh) %*% gg)


     if(type == "ridge" | type=="none"){
       update.pars <- update.pars
     }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
       for(j in pars_pen){
         update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
       }
     }else if(type=="diff_lasso" & lambda > 0){
       cc=0
       for(j in pars_pen){
         cc <- cc + 1
         #print(update.pars[j])
         #print(pen_diff[j])
         #print(soft(pen_diff[j],lambda,type="lasso",step=alpha1,e_alpha,gamma))
         update.pars[j] <- update.pars[j] -
           pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
       }
     }else if(type=="alasso" & lambda > 0){
       for(j in pars_pen){
         #print(update.pars[j])
         update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
       }
     }

      # S
     # gg2 <- grad(update.pars)
    #  hh2 <- hess(update.pars)
    #  minS <- min(mats$S[mats$S !=0])
    #  maxS <- max(mats$S)
   #   print(eigen(hh2)$values)

     # update.pars[min(mats$S[mats$S !=0]):max(mats$S)] <-
     #   update.pars[minS:maxS] - alpha2*(solve(hh2[minS:maxS,minS:maxS])%*% gg2[minS:maxS])


     # print(min(mats$S[mats$S !=0]))
     # print(max(mats$S))


    }else if(solver==TRUE){

      out <- nlminb(new.pars[count,],func,grad,control=list(iter.max=1,step.min=alpha,step.max=alpha))
      #out <- Rsolnp::solnp(new.pars[count,],func,control=list(trace=0))
      #print(out$objective)
      update.pars <- out$par #new.pars[count,] - alpha*(new.pars[count,] - out$par)


      if(type == "ridge" | type=="none"){
        update.pars <- update.pars
      }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
        for(j in pars_pen){
          update.pars[j] <- soft(update.pars[j],lambda,type,step=alpha,e_alpha,gamma)
        }
      }else if(type=="diff_lasso" & lambda > 0){
        cc=0
        for(j in pars_pen){
          cc <- cc + 1
          #print(update.pars[j])
          #print(pen_diff[j])
          #print(soft(pen_diff[j],lambda,type="lasso",step=alpha1,e_alpha,gamma))
          update.pars[j] <- update.pars[j] -
            pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
        }
      }else if(type=="alasso" & lambda > 0){
        for(j in pars_pen){
          #print(update.pars[j])
          update.pars[j] <- soft(pen_vec[j],lambda,type,step=alpha,e_alpha,gamma)
        }
      }

      # S



    }

    if(momentum==FALSE){
      new.pars[count+1,] <- update.pars
    }else if(momentum==TRUE){
      new.pars[count+1,] <- update.pars + (count/(count+3))*(update.pars-new.pars[count,])
    }



    if(type != "diff_lasso") pen_diff <- 0
    vals[count+1] = func(new.pars[count+1,])
    #print(round(vals[count+1],3))

    st.crit = try(abs(vals[count+1] - vals[count])<tol)
  #  st.crit2 <- all(abs(gg) < .01)
    dif <- abs(vals[count+1] - vals[count])
    #print(dif)
  #print(count)
   # print(rbind(update.pars,t(gg2)))
    #print(vals[count+1])
    #print(round(dif,5))
   # print(as.vector(round(gg,3)))
   # print(dif)
   # print(count)
  #  print(round(gg,3))
  #  print(convergence)

    if(inherits(st.crit, "try-error")){
      convergence=99
      print(9999)
    }else if(is.na(st.crit)==TRUE){
      convergence=99
      break
      print(8888)

    }else if(any(new.pars[count+1,] > par.lim[2]) | any(new.pars[count+1,] < par.lim[1])){
      break
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
