#' Plot function for cv_regsem
#'
#' @param object An object from cv_regsem.
#' @param pars Which parameters to plot
#' @param color Vector of colors for each parameter
#' @param ... Other arguments.
#' @export



plot_cv <- function(object,pars,color=NULL,...){


  coef.mat <- object[[1]][,pars]

  if(is.null(color)){
    set.seed(1)
    colls <- sample(colorspace::rainbow_hcl(length(pars)))
  }else{
    colls <- color
  }



  plot(object[[2]][,"lambda"],coef.mat[,1],ylim=c(min(coef.mat)-.1,max(coef.mat)+.1),
       ylab="Loading",xlab="Penalty",type="l",lty=1,col=colls[1],lwd=3)
  for(i in 2:(ncol(coef.mat))){
    lines(object[[2]][,"lambda"],coef.mat[,i],lty=1,col=colls[i+1],lwd=3)
    #points(mat[,(nload)])
  }

  abline(a=0,b=0)
}
