#' Plot function for cv_regsem
#'
#' @param x An x from cv_regsem.
#' @param ... Other arguments.
#' @param pars Which parameters to plot
#' @param show.minimum What fit index to use
#' @param col A specification for the default plotting color.
#' @param type what type of plot should be drawn. Possible types are "p" for points, "l" for lines, or "b" for both
#' @param lwd line width
#' @param h_line Where to draw horizontal line
#' @param lty line type
#' @param xlab X axis label
#' @param ylab Y axis label
#' @param legend.x x-coordinate of legend. See ?legend
#' @param legend.y y-coordinate of legend. See ?legend
#' @param legend.cex cex of legend. See ?legend
#' @param legend.bg legend background color. See ?legend
#' @param grey.out Add grey to background
#' @method plot cvregsem
#' @export



plot.cvregsem <- function (x, ..., pars = NULL, show.minimum="BIC",
                              col = NULL, type = "l", lwd = 3,h_line=0,
                              lty = 1, xlab = NULL, ylab = NULL,
                              legend.x = NULL, legend.y = NULL,
                              legend.cex = 1, legend.bg=par("bg"), grey.out=FALSE)
{
  if (is.null(pars))
    pars <- x$pars_pen
  if (!(type %in% c("p", "b", "l")))
    stop("Unknown plot type given.")
  if (!class(x) == "cvregsem")
    stop("Specified x is not a x from cv_regsem(.")
  if (is.null(xlab))
    xlab <- "Penalty"
  if (is.null(ylab))
    ylab <- "Estimate"
  coef.mat <- x$parameters[, pars]
  if (is.null(col)) {
    colls <- colorspace::rainbow_hcl(length(pars))
  }
  else {
    if (length(col) < length(pars))
      col <- rep(col, ceiling(length(pars)/length(col)))
    colls <- col
  }

  # filter NA values in fit function
 if(is.null(dim(coef.mat))){
   ydat <- coef.mat
 }else{
   ydat <- coef.mat[, 1]
 }
  xdat <- x$fits[, "lambda"]
  rm.ids <- which(x$fits[,"conv"] != 0)
  if (length(rm.ids)>0) {
    xdat <- xdat[-rm.ids]
    ydat <- ydat[-rm.ids]
    coef.mat <- coef.mat[-rm.ids, ]
  }

  # grey-out colors based on threshold
  if (grey.out!=FALSE) {
    if (isTRUE(grey.out)) grey.out<-0.001
    min.id <- which.min(abs(x$fits[ x$fits[,"conv"]!=0 ,show.minimum]))
    min.pars <- coef.mat[min.id, ]
    min.pars.filt <- abs(min.pars)<grey.out
    colls[min.pars.filt] <- "grey"
  }

  # adjust plot limits relative to scale not by absolute increment
  plot(xdat, ydat, ylim = c(min(coef.mat) * 0.95, max(coef.mat) * 1.05),
       ylab = ylab, xlab = xlab,
       type = "n")

  if(is.null(dim(coef.mat))){

      if (type == "l" || type == "b")
        lines(xdat, coef.mat, lty = lty,
              col = colls, lwd = lwd)
      if (type == "p" || type == "b")
        points(xdat, coef.mat)

  }else{
    for (i in 1:(ncol(coef.mat))) {
      if (type == "l" || type == "b")
        lines(xdat, coef.mat[, i], lty = lty,
              col = colls[i], lwd = lwd)
      if (type == "p" || type == "b")
        points(xdat, coef.mat[, i])
    }
  }

  # draw horizontal line
  abline(a=h_line,b=0)

  # add minimum
  if (!is.null(show.minimum)) {
    min.id <- which.min(abs(x$fits[,show.minimum]))
    lambda <- x$fits[min.id,1]

    abline(v=lambda,lty=2)

    pnts <- x$parameters[min.id,pars]
    points(rep(lambda,length(pnts)),pnts, col=colls,cex=2, pch=19)
    points(rep(lambda,length(pnts)),pnts, col="black",cex=1)
  }
  # --

  # add legend
  if (!is.null(legend.x)) {
    graphics::legend(x=legend.x, y=legend.y,colnames(x$parameters[,pars]),
                     #lty=1,lwd=legend.lwd,
                     fill=colls,
                     cex=legend.cex,
                     y.intersp=0.75,x.intersp=0.5,
                     bg=legend.bg
                     )
  }
}
