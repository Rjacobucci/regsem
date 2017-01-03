#' Plot function for cv_regsem
#'
#' @param object An object from cv_regsem.
#' @param pars Which parameters to plot
#' @param color Whether to plot in color
#' @param ... Other arguments.
#' @export



plot_cv <- function(object,pars,color=TRUE,...){


  coef.mat <- object[[1]][,pars]

  if(color==TRUE){

    colls <- sample(colorspace::rainbow_hcl(length(pars)))
  }else{
    colls <- rep(1,length(pars))
  }



  plot(object[[2]][,"lambda"],coef.mat[,1],ylim=c(min(coef.mat)-.1,max(coef.mat)+.1),
       ylab="Loading",xlab="Penalty",type="l",lty=1,col=colls[1],lwd=3)
  for(i in 2:(ncol(coef.mat))){
    lines(object[[2]][,"lambda"],coef.mat[,i],lty=1,col=colls[i+1],lwd=3)
    #points(mat[,(nload)])
  }

  abline(a=0,b=0)
}
