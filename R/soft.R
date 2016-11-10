

soft <- function(par,lambda,type,step,e_alpha=0.5){
  if(type=="lasso"){

      ret.val <- sign(par)*max(abs(par)-step*lambda,0)

  }else if(type=="enet"){

    stop("currently not supported")
    lambda1 = e_alpha*lambda
    lambda2 = (1-e_alpha)*lambda

    ret.val <- (sign(par)*(max(abs(par)-step*lambda/2,0)))/(1+step*lambda2)

  }else if(type=="alasso"){
    # ftp://ftp.stat.math.ethz.ch/Teaching/buhlmann/advanced-comput-statist/notes1.pdf
    ret.val <- sign(par)*max(abs(par)-(step*lambda)/(2*abs(par)),0)

  }else if(type=="scad"){

    stop("currently not supported")

    #if (abs(par) <= gm * (1+wq)) {
   #   par <- sign(par) * max(abs(par) - lambda, 0)
   # } else if (abs(par) >= alpha*lambda & abs(par) <= 2*lambda) {
   #   par <- sign(par) * max(abs(par) - gm * wq * dt / (dt - 1), 0) / (1 - (wq / (dt - 1)))
   # } else {}

  }else if(type=="mcp"){

    stop("currently not supported")

  #  if (abs(par) <= gm * dt) {
  #    ttq <- sign(par) * max(abs(par) - gm * wq, 0) / (1 - (wq / dt))
  #  } else {}

  }
  ret.val
}