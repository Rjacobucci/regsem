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
#' @param show.minimum What fit index to use
#' @method plot cvregsem
#' @export



plot.cvregsem <- function (x, ..., pars = NULL, col = NULL, type = "l", lwd = 3,
                              lty = 1, xlab = NULL, ylab = NULL, show.minimum=NULL)
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
  coef.mat <- x[[1]][, pars]
  if (is.null(col)) {
    colls <- colorspace::rainbow_hcl(length(pars))
  }
  else {
    if (length(col) < length(pars))
      col <- rep(col, ceiling(length(pars)/length(col)))
    colls <- col
  }

  # filter NA values in fit function
  ydat <- coef.mat[, 1]
  xdat <- x[[2]][, "lambda"]
  rm.ids <- which(x[[2]][,"conv"] != 0)
  if (length(rm.ids)>0) {
    xdat <- xdat[-rm.ids]
    ydat <- ydat[-rm.ids]
    coef.mat <- coef.mat[-rm.ids, ]
  }

  # adjust plot limits relative to scale not by absolute increment
  plot(xdat, ydat, ylim = c(min(coef.mat) * 0.95, max(coef.mat) * 1.05), ylab = ylab, xlab = xlab,
       type = "n")
  for (i in 1:(ncol(coef.mat))) {
    if (type == "l" || type == "b")
      lines(xdat, coef.mat[, i], lty = lty,
            col = colls[i], lwd = lwd)
    if (type == "p" || type == "b")
      points(xdat, coef.mat[, i])
  }
  abline(a = 0, b = 0)

  # add minimum
  if (!is.null(show.minimum)) {
    min.id <- which.min(x[[2]][,show.minimum])
    lambda <- x[[2]][min.id,1]

    abline(v=lambda,lty=2)

    pnts <- x[[1]][min.id,pars]
    points(rep(lambda,length(pnts)),pnts, col=colls,cex=2, pch=19)
    points(rep(lambda,length(pnts)),pnts, col="black",cex=1)
  }
  # --
}
