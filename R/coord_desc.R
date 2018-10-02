
coord_desc <- function(start,func,type,grad,hess,hessFun,pars_pen,model,lambda,mats,
                       block,max.iter,tol,full,solver,solver.maxit,alpha.inc,step,
                       step.ratio,diff_par,pen_vec,e_alpha,gamma,momentum,par.lim,quasi,
                       line.search,pen_vec_ml,prerun){
  count = 0
  ret <- list()
  max.iter = max.iter
  tol=tol
  #solver=TRUE
  phi_func <- rep(1,max.iter)
  phi_grad <- rep(1,max.iter)





  ## remove !!!!!!!!!!!!!!

  # !!!!!!!!!!!!!!!!!!


 # if(type=="enet"){
 #   step=step*2
 # }

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




  if(prerun==TRUE){

     out <- try(Rsolnp::solnp(start,func,control=list(trace=0)),silent=TRUE)
     if(inherits(out, "try-error")){
       out.solution <- start
     }else{
       out.solution <- out$pars
      # print(func(out.solution))
     }



   # out <- hydroPSO::hydroPSO(start,fn=func,lower=start-5,upper=start+5,
   #                           control=list(verbose=FALSE,write2disk=FALSE))
   # out.solution <- out$par
  }






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

    if(solver==FALSE & prerun==FALSE){

        if(full==TRUE & quasi == FALSE & hessFun=="none"){

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
                 update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha,e_alpha,gamma)
               }
             }

              func(update.pars)
          } #end delta

          if(line.search==TRUE){
            c=.5
            p=cbind(rep(0.5,length(new.pars[count,])))#update.pars-new.pars[count,]

            alpha=1
            if(count==1){
              alpha =step =1
              #
            }else{
              while(delta(alpha) > func(update.pars)+c*alpha*(t(gg)%*%p)){
                alpha = 0.8*alpha
              }
            }
          }else{
            alpha=alpha
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
               update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha,e_alpha,gamma)
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
            #alpha.vec[count] <- s1 <- alpha<- step

            if(hessFun != "none"){
              H <- solve(hess(new.pars[count,]))
            }else{
              H <-  diag(length(new.pars[count,])) #as.numeric((t(y)%*%s)/t(y)%*%y)
            }
            dir <- -H %*% grad.vec[count,]
          }else{
            s = cbind(new.pars[count,] - new.pars[count-1,])
            y =  cbind(grad.vec[count,] - grad.vec[count-1,])

            #http://terminus.sdsu.edu/SDSU/Math693a/Lectures/18/lecture.pdf


            #p <- as.numeric(1/(t(y)%*%s)) # add to rcpp
            Imat = diag(length(new.pars[count,]))

            #H <- (Imat - p*s%*%t(y))%*%H%*%(Imat-p*y%*%t(s)) + p*s%*%t(s)

            H <- rcpp_quasi_calc(Imat,s,y,H)$H2



            dir <- -H %*% grad.vec[count,]

            #alpha=step

          }


          alpha = 1
          update.pars <- new.pars[count,] + alpha*dir

          # print(out$objective)
          #update.pars <- out$par# - new.pars[count,]) + new.pars[count,]



            if(type == "ridge" | type=="none"){
              update.pars <- update.pars
            }else if(type!="none" & type!="ridge" & type!="diff_lasso" & lambda > 0){
              for(j in pars_pen){
                update.pars[j] <- soft(update.pars[j],lambda,type,step=1,e_alpha,gamma)
              }
            }else if(type=="diff_lasso" & lambda > 0){
              cc=0
              for(j in pars_pen){
                cc <- cc + 1
                #print(update.pars[j])
                #print(pen_diff[j])
                #print(soft(pen_diff[j],lambda,type="lasso",step=alpha1,e_alpha,gamma))
                update.pars[j] <- update.pars[j] -
                  pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=1,e_alpha,gamma)
              }
            }else if(type=="alasso" & lambda > 0){
              for(j in pars_pen){
                #print(update.pars[j])
                update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=1,e_alpha,gamma)
              }
            }



            v <- update.pars - new.pars[count,]

          # http://www.stat.cmu.edu/~ryantibs/convexopt-S15/lectures/24-prox-newton.pdf


            h <- function(pars){
              if(type=="enet"){
                (1-e_alpha)*sum(abs(pars[pars_pen])) + (e_alpha)*sqrt(sum(pars[pars_pen]**2))
              }else if(type=="ridge"){
                1*sqrt(sum(pars[pars_pen]**2))
              }else if(type=="lasso"){
                1*sum(abs(pars[pars_pen]))
              }else if(type=="alasso"){
                sum(abs(pars[pars_pen])*abs(1/pen_vec_ml))
              }else if(type=="diff_lasso"){
                sum(abs(pars[pars_pen] - pen_diff))
              }else{
                stop("quasi is currently not supported for that type of penalty")
              }
            }

            if(line.search==FALSE){
              alpha=step
            }else{
              vv <- t(grad.vec[count,])%*%v
              fmin.old <-  func(new.pars[count,])
              soft.old <- h(new.pars[count,])
              alpha= alpha
              c = 0.5 # previously 0.001 worked well; removed 0.5 0.8.1 set at 0.05
              # p. 102 of "sparsity" recommend 0.5, and then 0.8*alpha below


              # not changing alpha
              bool=FALSE
              while(bool==FALSE){

                if(alpha < 0.1){
                  alpha = 0.1; break
                }

                if(is.na(func(new.pars[count,]+alpha*v))){
                  alpha=.1
                  break
                }else if(is.na(fmin.old+c*alpha*(vv))){
                  alpha = .1
                  break
                }else if(is.na(c*((h(new.pars[count,]+alpha*v)-soft.old)))){
                  alpha = .1
                  break
                }else{
                  if(func(new.pars[count,]+alpha*v) > fmin.old+c*alpha*(vv) + c*((h(new.pars[count,]+alpha*v)-soft.old))){
                    alpha = 0.8*alpha
                    bool=FALSE
                  }else{
                    bool=TRUE
                  }

                }


              }
            }

#alpha = .2

          update.pars <- new.pars[count,] + alpha*(update.pars-new.pars[count,])


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
              update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha1,e_alpha,gamma)
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
    }else if(hessFun!="none" & solver==TRUE & prerun==FALSE){
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
         update.pars[j] <- update.pars[j] -
           pen_diff[cc] - soft(pen_diff[cc],lambda,type="lasso",step=alpha,e_alpha,gamma)
       }
     }else if(type=="alasso" & lambda > 0){
       for(j in pars_pen){
         #print(update.pars[j])
         update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha,e_alpha,gamma)
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


    }else if(solver==TRUE & hessFun=="none" & prerun==FALSE){

      out <- nlminb(new.pars[count,],func,grad,control=list(iter.max=1,step.min=alpha,step.max=alpha))
      #out <- lbfgs::lbfgs(func,grad,new.pars[count,],invisible=1)
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
          update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha,e_alpha,gamma)
        }
      }

      # S



    }else if(prerun==TRUE){

#print(out.solution)
     # calc2 = function(start){
      #  10-calc(start)
     # }
      #out = GA::ga("real-valued", fitness = func,min=0,max=1000, nBits = length(start),maxiter=30)


      gg <- try(grad(out.solution),silent=TRUE)
      #print(round(gg,2))
      if(inherits(gg, "try-error")) {
        gg <- rnorm(length(new.pars[count,]),0,.01)
      }else{
        gg <- gg
      }





      update.pars <- out.solution - alpha*gg

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
          update.pars[j] <- soft(update.pars[j]*(1/pen_vec_ml[j]),lambda,type,step=alpha,e_alpha,gamma)
        }
      }



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
     # print(9999)
    }else if(is.na(st.crit)==TRUE){
      convergence=99
      break
      #print(8888)

    }else if(any(new.pars[count+1,] > par.lim[2]) | any(new.pars[count+1,] < par.lim[1])){
      break
      convergence=99
    }else if(abs(vals[count+1])>100){
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
