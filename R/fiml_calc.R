#'
#'
#' Calculates the FIML objective function
#' @param ImpCov expected covariance matrix.
#' @param data required dataset.
#' @param Areg A matrix with current parameter estimates.
#' @param lambda penalty value.
#' @param alpha mixture for elastic net.
#' @param type penalty type.
#' @param pen_vec vector of penalized parameters.
#' @param nvar number of variables.
#' @keywords fit maximum likelihood regularization
#'
#' @examples
#' \dontrun{
#' fiml_calc()
#' }



fiml_calc = function(ImpCov,data,Areg,lambda,alpha,type,pen_vec,nvar){


  m = dim(ImpCov)[1]
  IntCol = which(colnames(Areg) == "1")
  IntCol2 = colnames(Areg[,1:nvar])
  fit = 0

  if(type=="none"){

    for(i in 1:nrow(data)){
      person1 = data[i,IntCol2]

      ind1 = which(is.na(person1)==T)
      misVar = colnames(data)[ind1]
      K = (nvar - sum(unique(ind1))) * log(2*pi)


      if(sum(is.na(person1)==T) >0){
        meanvec = Areg[(rownames(Areg)!="1"),IntCol]
        meanvec2 = meanvec[names(meanvec) != misVar]
        sub1 = c(person1[is.na(person1)==FALSE],0) - meanvec2
        indFit = K - log(det(ImpCov[-ind1,-ind1])) + t(sub1) %*% solve(ImpCov[-ind1,-ind1]) %*% sub1
      }else{
        meanvec = Areg[(rownames(Areg)!="1"),IntCol]
        sub1 = as.numeric(cbind(person1,0) - meanvec)
        indFit = K - log(det(ImpCov))  + t(sub1) %*% solve(ImpCov) %*% sub1
      }
      fit = fit + indFit
    }


  }else{
    stop("Only type==none is currently supported")
  }


  #-2*fit
  -2 * fit
}
