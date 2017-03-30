#' Plot function for cv_regsem
#'
#' @param x An x from cv_regsem.
#' @param ... Other arguments.
#' @param pars Which parameters to plot
#' @param col A specification for the default plotting color.
#' @param type what type of plot should be drawn. Possible types are "p" for points, "l" for lines, or "b" for both
#' @param lwd line width
#' @param lty line type
#' @param xlab X axis label
#' @param ylab Y axis label
#' @method plot cvregsem
#' @export



plot.cvregsem <- function(x,..., pars=NULL,col=NULL,type="l",lwd=3,lty=1,xlab=NULL,ylab=NULL){


  if(is.null(pars)) pars<-x$pars_pen



  # check user input for validity
  if (!(type %in% c("p","b","l"))) stop("Unknown plot type given.")
  if (!class(x)=="cvregsem") stop("Specified x is not a x from cv_regsem(.")

  if (is.null(xlab)) xlab <- "Penalty"
  if (is.null(ylab)) ylab <- "Loading"
  coef.mat <- x[[1]][,pars]
  # determine colors either from rainbow or by repeating given colors
  if(is.null(col)){
    colls <- colorspace::rainbow_hcl(length(pars))
  }else{
    if (length(col) < length(pars)) col <- rep(col, ceiling(length(pars)/length(col)) )
    colls <- col

  }

  # empty plot
  plot(x[[2]][,"lambda"],coef.mat[,1],ylim=c(min(coef.mat)-.1,max(coef.mat)+.1),
       ylab=ylab,xlab=xlab,type="n")

  # add lines/points
  for(i in 1:(ncol(coef.mat))){
    if (type == "l" || type == "b")
      lines(x[[2]][,"lambda"],coef.mat[,i],lty=lty,col=colls[i],lwd=lwd)
    if (type == "p" || type == "b")
      points(x[[2]][,"lambda"],coef.mat[,i])
  }

  abline(a=0,b=0)
}
