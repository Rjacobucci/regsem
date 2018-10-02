

soft <- function(par,lambda,type,step,e_alpha,gamma){
  if(type=="lasso"){
    lambda <- lambda*step
   # print(par)
      ret.val <- sign(par)*max(abs(par)-lambda,0)
     # print(ret.val)
  }else if(type=="enet"){
    #http://www.stat.washington.edu/courses/stat527/s13/readings/zouhastie05.pdf ****p. 305****
      lambda <- lambda*step
      lambda2 <- e_alpha*(lambda)
      lambda1 <- (1-e_alpha)*lambda
      ret.val <- sign(par)*(max(abs(par)-(lambda1),0)/(1+lambda2))

    #  ret.val <- (sign(par)*max(abs(par)-step*lambda,0))/(1+lambda2)

    #  if(abs(par) < e_alpha*lambda){
     #   ret.val <- 0
     # }else{
        # might be missing max(0,lambda)
    #    ret.val <- (sign(par)*(abs(par)-e_alpha*lambda))/(1+(1-e_alpha)*lambda)
    #    ret.val <- sign(par)*(max(abs(par)-(step*lambda)/2,0)/(1+step*lambda))
    #  }

  }else if(type=="alasso"){
    # ftp://ftp.stat.math.ethz.ch/Teaching/buhlmann/advanced-comput-statist/notes1.pdf
    ret.val <- sign(par)*max(abs(par)-(step*lambda)/(2*abs(par)),0)


  }else if(type=="scad"){
    lambda <- lambda*step
    gamma <- gamma*step
    #stop("currently not supported")



    if(abs(par) <= 2*lambda){
      ret.val <- sign(par)*max(abs(par)-lambda,0)

    #  ret.val <- sign(par)*max(abs(par)-lambda*gamma,0)
    }else if(abs(par) > 2*lambda & abs(par) <= lambda*gamma){
      ret.val <- ((gamma - 1)/(gamma - 2)) * sign(par)*max(abs(par)-((lambda*gamma)/(gamma-1)),0)

     # ret.val <- ((gamma-1)*par - sign(par)*gamma*lambda)/(gamma-2)
    }else if(abs(par) > lambda*gamma){
      ret.val <- par

    }

  }else if(type=="mcp"){

    lambda <- lambda*step
    gamma <- gamma*step
    #print(lambda*gamma)
    #stop("currently not supported")

    if(abs(par) <= lambda * gamma){
      ret.val <- (gamma/(gamma-1)) * sign(par)*max(abs(par)-lambda,0)
    }else if(abs(par) > lambda*gamma){
      ret.val <- par
    }

  }
  ret.val
}
