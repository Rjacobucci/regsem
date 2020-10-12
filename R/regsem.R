#'
#'
#' Regularized Structural Equation Modeling. Tests a single penalty. For
#' testing multiple penalties, see cv_regsem().
#'
#' @param model Lavaan output object. This is a model that was previously
#'        run with any of the lavaan main functions: cfa(), lavaan(), sem(),
#'        or growth(). It also can be from the efaUnrotate() function from
#'        the semTools package. Currently, the parts of the model which cannot
#'        be handled in regsem is the use of multiple group models, missing
#'        other than listwise, thresholds from categorical variable models,
#'        the use of additional estimators other than
#'        ML, most notably WLSMV for categorical variables. Note: the model
#'        does not have to actually run (use do.fit=FALSE), converge etc...
#'        regsem() uses the lavaan object as more of a parser and to get
#'        sample covariance matrix.
#' @param lambda Penalty value. Note: higher values will result in additional
#'        convergence issues. If using values > 0.1, it is recommended to use
#'        mutli_optim() instead. See \code{\link{multi_optim}} for more detail.
#' @param alpha Mixture for elastic net. 1 = ridge, 0 = lasso
#' @param gamma Additional penalty for MCP and SCAD
#' @param type Penalty type. Options include "none", "lasso",
#'        "enet" for the elastic net,
#'        "alasso" for the adaptive lasso
#'        and "diff_lasso". If ridge penalties are desired, use type="enet" and
#'        alpha=1. diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par. Two methods for
#'        sparser results than lasso are the smooth clipped absolute deviation,
#'        "scad", and the minimum concave penalty, "mcp". Last option is "rlasso"
#'        which is the randomised lasso to be used for stability selection.
#' @param dual_pen Two penalties to be used for type="dual", first is lasso, second ridge
#' @param random.alpha Alpha parameter for randomised lasso. Has to be between
#'        0 and 1, with a default of 0.5. Note this is only used for
#'        "rlasso", which pairs with stability selection.
#' @param data Optional dataframe. Only required for missing="fiml" which
#'        is not currently working.
#' @param optMethod Solver to use. Two main options for use: rsoolnp and coord_desc.
#'        Although slightly slower, rsolnp works much better for complex models.
#'        coord_desc uses gradient descent with soft thresholding for the type of
#'        of penalty. Rsolnp is a nonlinear solver that doesn't rely on gradient
#'        information. There is a similar type of solver also available for use,
#'        slsqp from the nloptr package. coord_desc can also be used with hessian
#'        information, either through the use of quasi=TRUE, or specifying a hess_fun.
#'        However, this option is not recommended at this time.
#' @param estimator Whether to use maximum likelihood (ML) or unweighted least squares
#'        (ULS) as a base estimator.
#' @param gradFun Gradient function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        Only for use with optMethod="coord_desc".
#' @param hessFun Hessian function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        This is currently not recommended.
#' @param prerun Logical. Use rsolnp to first optimize before passing to
#'        gradient descent? Only for use with coord_desc.
#' @param parallel Logical. Whether to parallelize the processes?
#' @param Start type of starting values to use. Only recommended to use
#'        "default". This sets factor loadings and variances to 0.5.
#'        Start = "lavaan" uses the parameter estimates from the lavaan
#'        model object. This is not recommended as it can increase the
#'        chances in getting stuck at the previous parameter estimates.
#' @param subOpt Type of optimization to use in the optimx package.
#' @param longMod If TRUE, the model is using longitudinal data? This changes
#'        the sample covariance used.
#' @param pars_pen Parameter indicators to penalize. There are multiple ways to specify.
#'        The default is to penalize all regression parameters ("regressions"). Additionally,
#'        one can specify all loadings ("loadings"), or both c("regressions","loadings").
#'        Next, parameter labels can be assigned in the lavaan syntax and passed to pars_pen.
#'        See the example.Finally, one can take the parameter numbers from the A or S matrices and pass these
#'        directly. See extractMatrices(lav.object)$A.
#' @param diff_par Parameter values to deviate from. Only used when
#'        type="diff_lasso".
#' @param LB lower bound vector. Note: This is very important to specify
#'        when using regularization. It greatly increases the chances of
#'        converging.
#' @param UB Upper bound vector
#' @param par.lim Vector of minimum and maximum parameter estimates. Used to
#'        stop optimization and move to new starting values if violated.
#' @param block Whether to use block coordinate descent
#' @param full Whether to do full gradient descent or block
#' @param calc Type of calc function to use with means or not. Not recommended
#'        for use.
#' @param nlminb.control list of control values to pass to nlminb
#' @param max.iter Number of iterations for coordinate descent
#' @param tol Tolerance for coordinate descent
#' @param round Number of digits to round results to
#' @param solver Whether to use solver for coord_desc
#' @param quasi Whether to use quasi-Newton
#' @param solver.maxit Max iterations for solver in coord_desc
#' @param alpha.inc Whether alpha should increase for coord_desc
#' @param line.search Use line search for optimization. Default is no, use fixed step size
#' @param step Step size
#' @param momentum Momentum for step sizes
#' @param step.ratio Ratio of step size between A and S. Logical
#' @param missing How to handle missing data. Current options are "listwise"
#'        and "fiml". "fiml" is not currently working well.
#' @return out List of return values from optimization program
#' @return convergence Convergence status. 0 = converged, 1 or 99 means the model did not converge.
#' @return par.ret Final parameter estimates
#' @return Imp_Cov Final implied covariance matrix
#' @return grad Final gradient.
#' @return KKT1 Were final gradient values close enough to 0.
#' @return KKT2 Was the final Hessian positive definite.
#' @return df Final degrees of freedom. Note that df changes with lasso
#'         penalties.
#' @return npar Final number of free parameters. Note that this can change
#'         with lasso penalties.
#' @return SampCov Sample covariance matrix.
#' @return fit Final F_ml fit. Note this is the final parameter estimates
#'         evaluated with the F_ml fit function.
#' @return coefficients Final parameter estimates
#' @return nvar Number of variables.
#' @return N sample size.
#' @return nfac Number of factors
#' @return baseline.chisq Baseline chi-square.
#' @return baseline.df Baseline degrees of freedom.
#' @keywords optim calc
#' @useDynLib regsem, .registration=TRUE
#' @import Rcpp
#' @import parallel
#' @import lavaan
#' @import Rsolnp
#' @importFrom stats cov na.omit nlminb pchisq rnorm runif sd uniroot var weighted.mean cov2cor quantile
#' @importFrom graphics abline lines plot points par
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' @examples
#' # Note that this is not currently recommended. Use cv_regsem() instead
#' library(lavaan)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ 1*x1 + l1*x2 + l2*x3 + l3*x4 + l4*x5 + l5*x6 + l6*x7 + l7*x8 + l8*x9
#' '
#' # Recommended to specify meanstructure in lavaan
#' outt = cfa(mod, HS, meanstructure=TRUE)
#'
#' fit1 <- regsem(outt, lambda=0.05, type="lasso",
#'   pars_pen=c("l1", "l2", "l6", "l7", "l8"))
#' #equivalent to pars_pen=c(1:2, 6:8)
#' #summary(fit1)





regsem = function(model,lambda=0,alpha=0.5,gamma=3.7, type="lasso",
                  dual_pen=NULL,
                  random.alpha=0.5,
                  data=NULL,optMethod="rsolnp",
                  estimator="ML",
                 gradFun="none",hessFun="none",prerun=FALSE,parallel="no",Start="lavaan",
                 subOpt="nlminb",longMod=F,
                 pars_pen="regressions",
                 diff_par=NULL,
                 LB=-Inf,
                 UB=Inf,
                 par.lim=c(-Inf,Inf),
                 block=TRUE,
                 full=TRUE,
                 calc="normal",
                 max.iter=500,
                 tol=1e-5,
                 round=3,
                 solver=FALSE,
                 quasi=FALSE,
                 solver.maxit=5,
                 alpha.inc=FALSE,
                 line.search=FALSE,
                 step=.1,
                 momentum=FALSE,
                 step.ratio=FALSE,
                 nlminb.control=list(),
                 missing="listwise"){

  e_alpha=alpha

  if(quasi==TRUE){
    warnings("The quasi-Newton method is currently not recommended")
  }




  if (class(model)!="lavaan") stop("Input is not a 'lavaan' object")

  match.arg(type,c("lasso","none","ridge","scad","alasso","mcp","diff_lasso","enet","rlasso","rlasso2","dual"))


 # if(type == "mcp" & optMethod!= "coord_desc"){
  #  stop("For both scad and mcp must use optMethod=coord_desc")
 # }

 # if(type == "scad" & optMethod != "coord_desc"){
 #   stop("For both scad and mcp must use optMethod=coord_desc")
 # }



  if(type=="ridge"){
    type="enet";alpha=1
  }


  mats = extractMatrices(model)

  if(estimator=="ML"){
    estimator2 = 1
  }else if(estimator=="ULS"){
    estimator2 = 2
  }


  # ------------------------------

  # create pars_pen

  # ------------------------------

  # turn parameter labels into numbers

  pars_pen2 = NULL
  if(any(pars_pen=="regressions") & is.null(mats$regressions)){
    stop("No regression parameters to regularize")
  }

  if(type=="dual"){
    pars_pen=pars_pen
  }else if(any(pars_pen == "loadings")){
    # inds = mats$A[,mats$name.factors]
    # pars_pen2 = c(inds[inds != 0],pars_pen2)
    pars_pen2 = mats$loadings
  }else if(any(pars_pen == "regressions") | is.null(pars_pen)){
    pars_pen2 = c(pars_pen2,mats$regressions)
    # if(is.na(mats$name.factors)==TRUE){
    #   if(any(colnames(mats$A) == "1")){
    #    IntCol = which(colnames(mats$A) == "1")
    #    A_minusInt = mats$A[,-IntCol]
    #    A_pen = A_minusInt != 0
    #    pars_pen2 = c(A_minusInt[A_pen],pars_pen2)
    #  }else{
    #    A_pen = mats$A != 0
    #    pars_pen2 = c(mats$A[A_pen],pars_pen2)
    #  }
    # }else{
    # remove factor loadings
    #   if(any(colnames(mats$A) == "1")){
    #     IntCol = which(colnames(A) == "1" | colnames(mats$A) != mats$name.factors)
    #     A_minusInt = mats$A[,-IntCol]
    #    A_pen = A_minusInt != 0
    #    pars_pen2 = c(A_minusInt[A_pen],pars_pen2)
    #  }else{
    #    inds2 = mats$A[,colnames(mats$A) != mats$name.factors]
    #
    #   pars_pen2 = c(inds2[inds2 != 0],pars_pen2)
    #  }
    # }
  }else if(is.null(pars_pen)==FALSE & is.numeric(pars_pen)==FALSE){
    #pars_pen2 <- parse_parameters(pars_pen,model)

    ids = which(mats$pars.align[,2] %in% pars_pen)

    if(length(ids) != length(pars_pen) & type != "dual"){
      stop("Have to specify parameter number in pars_pen for equality constrained labels")
    }

    pars_pen2 = as.numeric(mats$pars.align[ids,1])

  }else if(is.numeric(pars_pen)){
    pars_pen2 = pars_pen
  }#else if(is.null(pars_pen)==TRUE){
  #  if(any(colnames(mats$A) == "1")){
  #    IntCol = which(colnames(mats$A) == "1")
  #    A_minusInt = mats$A[,-IntCol]
  #    A_pen = A_minusInt != 0
  #    pars_pen2 = A_minusInt[A_pen]
  #  }else{
  #    A_pen = mats$A != 0
  #    pars_pen2 = mats$A[A_pen]
  #  }
  #}


if(type != "dual"){
  pars_pen = as.numeric(pars_pen2)
}

if(type=="rlasso"){
  ralpha <- runif(length(pars_pen),random.alpha,1) # can alter and add argument
}
if(type=="rlasso2"){
  ralpha <- runif(length(pars_pen),0,random.alpha) # can alter and add argument
}

#  if(optMethod=="nlminb"& type !="ridge" | type != "none"){
#    stop("Only optMethod=coord_desc is recommended for use")
#  }

  if(length(nlminb.control)==0){
    nlminb.control <- list(abs.tol=1e-6,
                    iter.max=60000,
                    eval.max=60000,
                    rel.tol=1e-6,
                    x.tol=1e-6,
                    xf.tol=1e-6)
  }


  if(missing=="fiml" & is.null(model@SampleStats@missing[[1]])){
    stop("need to change missing=fiml in lavaan")
  }

#  if(gradFun != "none" & missing=="fiml"){
#    stop("only gradFun = none is supported with missing data")
#  }

#  if(model@Data@nobs[[1]] != model@Data@norig[[1]]){
#    warning("regsem is currently not working well in the presence of missing data")
#  }


  #  if(gradFun=="norm"){
  #    stop("Only recommended grad function is ram or none at this time")
 #   }

#  if(type=="ridge" & gradFun != "none"){
#    warning("At this time, only gradFun=none recommended with ridge penalties")
#  }

 # if(type=="ridge" & optMethod != "nlminb"){
 #   stop("For ridge, only use optMethod=nlminb and gradFun=none")
 # }

  if(type=="lasso"  & gradFun != "ram"){
    warning("At this time, only gradFun=ram recommended with lasso penalties")
  }


  if(type=="diff_lasso"  & gradFun != "ram"){
    warning("At this time, only gradFun=ram recommended with lasso penalties")
  }
 #   parL = parTable(model)[,"label"]
  #  if(sum(duplicated(parL[parL != ""])) > 0){
 #     stop("regsem currently does not allow equality constraints")
 #   }


  if(model@SampleStats@ngroups > 1){
    stop("regsem currently does not work with multiple group models")
  }


 # lav.fits <- fitmeasures(model)

#  if(model@Fit@converged == FALSE){
#    sat.lik = NA
#  }else{
#    sat.lik <- as.numeric(lav.fits["unrestricted.logl"])
#  }





    nvar = model@pta$nvar[[1]][1]
    nfac = model@pta$nfac[[1]][1]

    if(missing=="listwise"){
      calc_fit = "cov"
      nobs = model@SampleStats@nobs[[1]][1]



      # get defined parameters
      if(any(is.na(mats$defined_params))==FALSE){
        defined_params_vals <- mats$defined_params
      }else{
        defined_params_vals <- NA
      }

      if(missing=="listwise"){
        SampCov <- model@SampleStats@cov[][[1]]
      }else{
        SampCov <- model@implied$cov[[1]]
      }


      if(mats$mean == TRUE){
        mm = mats$A[,"1"]

        SampMean <- model@SampleStats@mean[][[1]]
        ss = match(names(mm[mm > 0]),model@Data@ov$name)
        ss <- ss[!is.na(ss)] # NAs may arise from fixed means
        SampMean[-c(ss)] = 0


        SampCov2 <- SampCov + SampMean%*%t(SampMean)
        # try changing size of SampCov
        SampCov3 = cbind(SampCov2,SampMean)
        SampCov = rbind(SampCov3,append(SampMean,1))
      }else if(mats$mean == FALSE){
       # SampCov <- model@SampleStats@cov[][[1]]
        SampMean = NULL
      }

      #for grad ram with mean
      SampCov22 <- model@SampleStats@cov[][[1]]# + SampMean %*% t(SampMean)

    }else if(missing=="fiml"){
      #stop("FIML is currently not supported at this time")
      calc_fit = "ind"
      SampCov <- model@SampleStats@cov[][[1]]
      defined_params_vals <- NA

      nobs = model@SampleStats@nobs[[1]][1]
     # if(is.null(data)==TRUE){
     #   stop("Dataset needs to be provided for missing==fiml")
     # }


     # if(length(model@ParTable$op[model@ParTable$op == "~1"]) == 0){
     #   stop("meanstructure needs to be equal to TRUE for FIML")
     # }

    }
    #SampCov <- fitted(model)$cov
    #SampMean <- rep(0,nvar)


  #  if(estimator=="ULS"){
   #   poly_vec = SampCov[lower.tri(SampCov)]
  #  }


    type2 = 0
    if(type=="lasso"){
      type2 = 1
    }else if(type=="ridge"){
      type2=2
    }else if(type=="diff_lasso"){
      type2=3
    }else if(type=="enet"){
      type2=4
    }else if(type=="alasso"){ ## try just creating new pen_vec
      type2=1
    }else if(type=="rlasso"){ ## try just creating new pen_vec
        type2=1
    }else if(type=="rlasso2"){ ## try just creating new pen_vec
      type2=5
    }else if(type=="scad"){
      type2=6
    }else if(type=="mcp"){
      type2=7
    }else if(type=="dual"){
      type2=8
    }





    #nUniq = nvar
    #nFacCov
    df = model@test[[1]]$df
    npar = model@Fit@npar
    nload = length(model@ParTable$op[model@ParTable$op == "=~"])




  A <- mats$A
  A_est <- mats$A_est
  A_fixed <- mats$A_fixed
  S <- mats$S
  S_est <- mats$S_est
  S_fixed <- mats$S_fixed
  F <- mats$F
  I <- diag(nrow(A))

 # pars_pen <- parse_parameters(pars_pen, model)





   if(class(Start)=="numeric"){
      start=Start
   }else if(class(Start) != "numeric"){
     if(Start=="lavaan"){
       # get starting values
       start <- mats$parameters

     } else if(Start == "default"){
       nstart <- max(max(A),max(S))
       start <- rep(0.5,nstart)

     }
   }


  pen_vec1 = 0;pen_vec2=0
  dual_pen1=0;dual_pen2=0

  if(calc == "normal"){
    calc = function(start){
         mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
         #print(mult)

         #mult2 = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
         #print(mult2)
         pen_vec = c(mult$A_est22[match(pars_pen,A,nomatch=0)],mult$S_est22[match(pars_pen,S,nomatch=0)])
         if(type=="diff_lasso"){
           pen_diff = pen_vec - diff_par
         }else{
           pen_diff=0
         }



         # two penalties
         if(type=="dual"){
           pen_vec1 = c(mult$A_est22[match(pars_pen[[1]],A,nomatch=0)],mult$S_est22[match(pars_pen[[1]],S,nomatch=0)])
           pen_vec2 = c(mult$A_est22[match(pars_pen[[2]],A,nomatch=0)],mult$S_est22[match(pars_pen[[2]],S,nomatch=0)])
           dual_pen1 = dual_pen[1];dual_pen2=dual_pen[2]
         }

         #### for alasso - weight the parameters ####
         #### overwrite pen_vec ########
         if(type=="alasso"){
           pen_vec_ml = c(mats$A_est[match(pars_pen,A,nomatch=0)],mats$S_est[match(pars_pen,S,nomatch=0)])
           pen_vec = abs(pen_vec)*(1/(abs(pen_vec_ml)))
         }
         if(type=="rlasso"){
           pen_vec = abs(pen_vec)/ralpha
         }
         if(type=="rlasso2"){
           rlasso_pen = sum((ralpha + lambda) * abs(pen_vec))
         }else{
           rlasso_pen=0
         }

         if(calc_fit=="cov"){

           #fit = fit_fun(ImpCov=mult$ImpCov,SampCov,Areg=mult$A_est22,lambda,alpha,type,pen_vec)
           if(estimator=="ULS"){
             imp_vec = mult$ImpCov[lower.tri(mult$ImpCov)]
           }

           fit = rcpp_fit_fun(ImpCov=mult$ImpCov,SampCov,type2,lambda,gamma,
                              pen_vec,pen_diff,e_alpha,rlasso_pen,pen_vec1,pen_vec2,
                              dual_pen1,dual_pen2)#,estimator2,poly_vec,imp_vec)
          # print(fit)
           #print(fit)
          # print(type2)
           #print(round(fit,3))#;print(pen_diff)
           #print(fit)
           fit
         }else if(calc_fit=="ind"){
           stop("Not currently supported")
           #print(mult$ImpCov)

           #fit = fiml_calc(ImpCov=mult$ImpCov,mu.hat=model@SampleStats@missing.h1[[1]]$mu,
           #                h1=model@SampleStats@missing.h1[[1]]$h1,
           #                Areg=mult$A_est22,lambda,alpha,type,pen_vec,nvar,
           #                lav.miss=model@SampleStats@missing[[1]])
          # fit = fiml_calc2(ImpCov=mult$ImpCov,F,mats2=mult,
          #                 type=type,lambda=lambda,
          #                 model=model,sat.lik=sat.lik,
          #                 pen_vec=pen_vec)
         }

    }
  }else if(calc == "calc2"){
    calc = function(start){
      #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
      fit = ram_calc(par=start,SampCov22,A,S,F,SampMean)
      like = fit$lik
      like
    }
  }else if(calc == "boot"){
    calc = function(start){
      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
      pen_vec = c(mult$A_est22[match(pars_pen,A,nomatch=0)],mult$S_est22[match(pars_pen,S,nomatch=0)])
      if(type=="diff_lasso"){
        pen_diff = pen_vec - diff_par
      }else{
        pen_diff=0
      }
      n.boot=100
      fits = rep(NA,n.boot)
      # boot part
      for(i in 1:n.boot){
        dat1 = model@Data@X[[1]]
        ids <- sample(nrow(dat1),nrow(dat1),replace=TRUE)
        new.dat <- dat1[ids,]
        SampCov.boot <- cov(new.dat)
        #fits[i] = rcpp_fit_fun(ImpCov=mult$ImpCov,SampCov.boot,type2,lambda,pen_vec,pen_diff)
        fits[i] = fit_fun(ImpCov=mult$ImpCov,SampCov.boot,lambda,alpha,type,pen_vec)
      }
      mean(fits)
    }
  }
#------------------------ only works for CFA models --------------------------------
#  grad = function(start){
#    mult = RAM_multSimp(A,A.fixed,S,S.fixed,F,start,nfac,nvar)
#    ret = gradient(ExpCov=mult$ExpCov,cov,A=mult$A,S=mult$S,
#                   start,lambda,alpha,type,nvar,nload,nfac,nUniq,nFacCov)
#    ret
# }

  if(gradFun=="norm"){
    grad = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
                #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
                ret = grad_fun(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                               Sreg=mult$S_est22,A,A_fixed,A_est,S,S_fixed,S_est,
                               F,lambda,alpha,type,pars_pen)
               ret
  }
} else if(gradFun=="ram"){
    grad = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)

      if(optMethod=="coord_desc"){
        if(type2==1 | type2==3 | type2 ==4 | type2==6 | type2==7) type2=0


        ret = 2*rcpp_grad_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                            Sreg=mult$S_est22,A,S,
                            F,lambda,type2=type2,pars_pen,diff_par=0)
      }else{

           ret = 2*rcpp_grad_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                           Sreg=mult$S_est22,A,S,
                             F,lambda,type2=type2,pars_pen,diff_par=0)
        }


      ret
    }

}else if(gradFun=="ram_mean"){
  grad = function(start){

    mult = ram_calc(par=start,SampCov22,A,S,F,SampMean)
    ret = grad_ram_wMean(par=start,ImpCov=mult$ImpCov,SampCov22,Areg = mult$A2,
                    Sreg=mult$S2,A=mult$A.pars,S=mult$S.pars,
                    F=mult$F,SampMean,lambda,type,m=mult$m,mu=mult$mu,m.pars=mult$m.pars)
    ret
  }
} else if(gradFun=="numDeriv"){
    grad = function(start){

    #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
    ret = numDeriv::grad(calc,x=start)
    ret
  }
} else if(gradFun=="auto"){
  grad = function(start){
stop("gradFun = auto is not supported at this time")

  #  ret = numderiv(calc,x=start)
  #  ret
  }
}


  if(hessFun=="norm"){



    hess = function(start){

        mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
        #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
        retH = hessian(par=start,ImpCov=mult$ImpCov,SampCov,A,A_fixed,A_est,
                 S,S_fixed,S_est,F)
        retH
    }

}  else if(hessFun=="ram"){


    hess = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
      ret = hess_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                      Sreg=mult$S_est22,A,S,F)
      ret
    }

} else if(hessFun=="numDeriv"){
  hess = function(start){

    #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
    ret = numDeriv::hessian(calc,x=start)
    ret
  }
}else{
  hess = NULL
}

    res <- list()
if(optMethod=="nlminb"){
    if(gradFun=="norm"){
      if(hessFun=="norm"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(abs(diag(S)-max(A)),rep(-10,max(S)-max(diag(S))))
        #UB = c(100,100,100,1)
        out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                 iter.max=max.iter))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
       #LB = c(rep(-6,max(A)),rep(1e-6,rep(-10,max(S)-max(diag(S))),max(diag(S))-max(A)))
        out <- nlminb(start,calc,grad,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                     iter.max=max.iter))
        res$out <- out

        res$convergence = out$convergence
        #res$optim_fit <- out$objective
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="ram"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
       out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=nlminb.control)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="norm"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=nlminb.control)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        res$iteration = out$iterations
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
       suppressWarnings(out <- nlminb(start,calc,grad,lower=LB,upper=UB,
                     control=nlminb.control)) #,x.tol=1.5e-6
        res$out <- out
        res$iteration = out$iterations
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,lower=LB,upper=UB,control=nlminb.control)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
    }else if(gradFun=="ram_mean"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,grad,hess,lower=LB,control=list(eval.max=max.iter,
                                                                 iter.max=max.iter))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc2,grad,lower=LB,upper=UB,eval.max=max.iter,
                      iter.max=max.iter)
        res$out <- out
        res$convergence = out$convergence
        #res$optim_fit <- out$objective
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="none"){
      #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
      out <- nlminb(start,calc,lower=LB,upper=UB,control=nlminb.control)
      res$out <- out
      #res$optim_fit <- out$objective
      res$convergence = out$convergence
      par.ret <- out$par
      #res$iterations <- out$iterations
    }else if(gradFun=="numDeriv"){
      if(hessFun=="numDeriv"){
        warning("numDeriv does not seem to be accurate at this time")
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=nlminb.control)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        warning("numDeriv does not seem to be accurate at this time")
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,lower=LB,upper=UB,gradient=grad,control=nlminb.control)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }
}else if(optMethod=="optimx"){
    if(gradFun=="norm"){
      if(hessFun=="norm"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals

        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,lower=LB,upper=UB,method=subOpt,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        res$convergence = out$convcode
        par.ret <- coef(out)
      }
    }else if(gradFun=="ram_mean"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,lower=LB,upper=UB,method=subOpt,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        res$convergence = out$convcode
        par.ret <- coef(out)
       }
    }else if(gradFun=="ram"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE,itnmax = max.iter))
        res$out <- out
        #res$optim_fit <- out$value
        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        #res$optim_fit <- out$value
        res$convergence = out$convcode
        par.ret <- coef(out)
      }
    }else if(gradFun=="none"){
      #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-S[1,1]+1),rep(-10,max(S)-max(diag(S))))
      out <- optimx::optimx(start,calc,method=subOpt,itnmax=1500,lower=LB,upper=UB,
                    hessian=T,control=list(maxit=1500,starttests=FALSE,all.methods=TRUE,abstol=1e-12))
      res$out <- out
      #res$iterations <- out$fevals
      res$convergence = out$convcode
      #pars <- coef(out)
      #res$pars <- pars
      par.ret <- coef(out)
    }
}else if(optMethod=="rsolnp"){
       # if(UB == Inf) UB=NULL
       # if(LB == -Inf) LB=NULL

        suppressWarnings(out <- Rsolnp::solnp(start,calc,#LB=LB,UB=UB,
                                              control=list(trace=0)))#tol=1e-16
        #out <- optim(par=start,fn=calc,gr=grad)
        res$out <- out
        #res$iterations <- out$nfunevals
        res$optim_fit <- out$values[length(out$values)]

        if(abs(res$optim_fit)>100){
          res$convergence=1
        }else{
          res$convergence = out$convergence
        }

        par.ret <- out$pars
       # print(par.ret)

}else if(optMethod=="slsqp"){

  out <- nloptr::slsqp(start,calc)

    par.ret <- out$par
   res$optim_fit <- out$value
   res$convergence = ifelse(out$convergence>1,0,1)

}else if(optMethod=="NlcOptim"){

  out <- NlcOptim::solnl(start,calc)

  par.ret <- out$par
  res$optim_fit <- out$fn
  res$convergence = 0#ifelse(out$convergence>1,0,1)

}else if(optMethod=="lbfgs"){

  suppressWarnings(out <- lbfgs::lbfgs(calc,grad,start,orthantwise_c =lambda,
                                       orthantwise_start=0,orthantwise_end=6))
  print(out)
  par.ret <- out$par
  res$optim_fit <- out$value
  res$convergence = out$convergence

}else if(optMethod=="GA"){
  calc2 = function(start){
    10-calc(start)
  }
  out = GA::ga("real-valued", fitness = calc2, nBits = length(start),
               min=LB,max=UB,monitor=FALSE,pcrossover=0.9,pmutation=0.3,
               maxiter=10000)
  res$out <- summary(out)
  res$optim_fit <- 10 - summary(out)$fitness
  res$convergence = 0
  res$par.ret <- summary(out)$solution
}else if(optMethod=="coord_desc"){

  if(type=="alasso"){
    mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
    pen_vec = c(mult$A_est22[match(pars_pen,A,nomatch=0)],mult$S_est22[match(pars_pen,S,nomatch=0)])
    pen_vec_ml = c(mats$A_est[match(pars_pen,A,nomatch=0)],mats$S_est[match(pars_pen,S,nomatch=0)])
    pen_vec = abs(pen_vec)*(1/(abs(pen_vec_ml)))
  }



  out = coord_desc(start=start,func=calc,type=type,grad=grad,
                   hess=hess,hessFun=hessFun,prerun=prerun,
                   pars_pen=pars_pen,model=model,max.iter=max.iter,
                   lambda=lambda,mats=mats,block=block,tol=tol,full=full,
                   solver=solver,solver.maxit=solver.maxit,
                   alpha.inc=alpha.inc,step=step,momentum=momentum,
                   e_alpha=e_alpha,gamma=gamma,
                   par.lim=par.lim,
                   step.ratio=step.ratio,diff_par=diff_par,pen_vec=pen_vec,quasi=quasi,
                   line.search=line.search,pen_vec_ml)
  res$out <- out
  res$optim_fit <- out$value
  #print(out$convergence)

    res$convergence <- out$convergence

  par.ret <- out$pars
  res$iterations <- out$iterations
}



   # imp_cov = RAMmult(res$out$par,A,S,F,A.fixed,A.est,S.fixed,S.est) [[1]]
   # res$imp_cov <- imp_cov
    #res$df <- df
    #res$nobs <- nobs
    #res$nload <- nload

    #hess <- ret_hess(res$out$par,A,S,F,A_fixed,A_est,S_fixed,S_est)
    #res$hess <- hess
    #res$info <- ginv((nobs/2)*hess)
    #res$A <- A
    #res$S <- S
    #res$A.est <- A_est
    #res$A.fixed <- A_fixed
    #res$S_est <- S_est
    #res$S.fixed <- S_fixed
   # res$F <- F


    pars.df <- data.frame(matrix(NA,1,max(max(A),max(S))))
    pars.df[1,] <- par.ret

#    if(any(pars.df[diag(S[diag(S) != 0])] < 0)){
#      warning("Some Variances are Negative!")
#      res$convergence <- 2
 #   }





    if(type=="ridge"){
      hess = function(start){

        mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
        #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
        ret = hess_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                       Sreg=mult$S_est22,A,S,F)
        ret
      }
      hess.mat = hess(as.numeric(pars.df))
      res$hess <- hess.mat
    }


    #res$ftt = rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      # get Implied Covariance matrix

    mult.out <- rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
    Imp_Cov1 <- mult.out$ImpCov
    #Imp_Cov <- RAMmult(par=as.numeric(pars.df),A,S,F,A_fixed,A_est,S_fixed,S_est)$ImpCov

    pen_vec = c(mult.out$A_est22[match(pars_pen,A,nomatch=0)],mult.out$S_est22[match(pars_pen,S,nomatch=0)])

    if(mats$mean==TRUE & missing=="listwise"){
      Imp_Cov = Imp_Cov1[1:(nrow(Imp_Cov1)-1),1:(ncol(Imp_Cov1)-1)] - SampMean %*% t(SampMean)
    }else{
      Imp_Cov = Imp_Cov1
    }

    res$Imp_Cov <- Imp_Cov

    res$Imp_Cov1 <- Imp_Cov1

     # N = nobs; p=nvar; SampCov00 <- model@SampleStats@cov[][[1]]
     # c <- N*p/2 * log(2 * pi)
      #res$logl_sat <- -c -(N/2) * log(det(SampCov00)) - (N/2)*p

#    if(model@Fit@converged == FALSE){
#      res$logl_sat= NA
#    }else{
 #     res$logl_sat <- as.numeric(lav.fits["unrestricted.logl"])
 #   }


    #res$grad <- grad(as.numeric(pars.df))
    #### KKT conditions #####
    if(gradFun=="none"){
      res$KKT1 = "grad not specified"
    }else{
      res$grad <- grad(as.numeric(pars.df))
      kk = try(all(grad(as.numeric(pars.df)) < 0.001))
      if(inherits(kk, "try-error")){
        res$KKT1 = "error"
      }else{
        if(kk == TRUE){
          res$KKT1 = TRUE
        }else if(kk < 0.001){
          res$KKT1 = FALSE
        }else{
          res$KKT1 = NA
        }
      }
    }


#    if(hessFun=="none"){
#      res$KKT2 = "hess not specified"
#    }else{
#      hess.mat = hess(as.numeric(pars.df))
#      eig = eigen(hess.mat)$values
#      hess.eigs = try(all(eig) > 1e-6)
#      if(inherits(hess.eigs, "try-error")){
#        res$KKT2 = "error"
 #     }else{
 #       if(hess.eigs == TRUE){
#          res$KKT2 = TRUE
#        }else if(hess.eigs == FALSE){
#          res$KKT2 = FALSE
#        }else{
#          res$KKT2 = NA
#        }
#      }
#    }




   # rettt = rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
   # A_new2 <- rettt$A_est22
   # S_new2 <- rettt$S_est22
    #A_new <- A
    #S_new <- S
    #I = diag(ncol(A))

    #for(i in 1:max(A)) A_new[A_new==i] = pars.df[i]
    #A_new2 = matrix(unlist(A_new),nrow(A),ncol(A))

    #for(i in min(S[S>0]):max(S)) S_new[S_new==i] = pars.df[i]
    #S_new2 = matrix(unlist(S_new),nrow(S),ncol(S))
    #res$retttt = calc(as.numeric(pars.df))
    #res$Imp_Cov = F %*% solve(I-A_new2) %*% S_new2 %*% t(solve(I-A_new2)) %*% t(F)
    #res$fitt =  calc(as.numeric(pars.df))

    for(i in 1:ncol(pars.df)){
      if(any(A == i)){
        pos = which(A == i,arr.ind=T)
        one = colnames(A)[pos[2]]
        two = rownames(A)[pos[1]]
        names(pars.df)[i] = paste(one,"->",two)
      }else if(any(S==i)){
        pos = which(S == i,arr.ind=T)
        one = colnames(S)[pos[2]]
        two = rownames(S)[pos[1]]
        names(pars.df)[i] = paste(one,"~~",two)
      }
    }


    if(type=="none" | lambda==0){
      res$df = df
      res$npar = npar
    }else if(type=="lasso" | type=="alasso" | type=="rlasso" | type=="rlasso2" | type=="enet" | type=="scad" | type=="mcp" & alpha < 1){
      #A_estim = A != 0
      #pars = A_est[A_estim]
      pars_sum = pars.df[pars_pen]
      pars_l2 = sqrt(pars_sum**2)
      res$df = df + sum(pars_l2 < 1/(10^round))
      res$npar = npar - sum(pars_l2 < 1/(10^round))

    }else if(type=="dual"){
      pars_sum = pars.df[pars_pen[[1]]]
      pars_l2 = sqrt(pars_sum**2)
      res$df = df + sum(pars_l2 < 1/(10^round))
      res$npar = npar - sum(pars_l2 < 1/(10^round))
    }else if(type=="ridge" | alpha == 1){
      #ratio1 <- sqrt(pars.df[pars_pen]**2)/sqrt(mats$parameters[pars_pen]**2)
      res$df = df #+ length(ratio1) - sum(ratio1)
      res$npar = npar #- sum(ratio1)
    }else if(type=="diff_lasso"){
      pars_sum = as.numeric(pars.df[pars_pen])
      #print(pars_sum);print(duplicated(round(pars_sum,3)))
      res$df = df + sum(duplicated(round(pars_sum,3)))
      res$npar = npar - sum(duplicated(round(pars_sum,3)))

    }


      if(optMethod=="nlminb"){
        optFit <- out$objective
      }else{
        optFit <- res$optim_fit
      }


    if(missing == "listwise"){
     # SampCov <- model@SampleStats@cov[][[1]]
    #  res$SampCov = SampCov
      #res$fit = 0.5*(log(det(Imp_Cov1)) + trace(SampCov %*% solve(Imp_Cov1)) -
       #              log(det(SampCov))  - nvar)
     # pen_vec = c(mult$A_est22[A %in% pars_pen],mult$S_est22[S %in% pars_pen])
      if(type=="diff_lasso"){
        pen_diff = pen_vec - diff_par
      }else{
        pen_diff=0
      }
      if(estimator=="ULS"){
        imp_vec = Imp_Cov1[lower.tri(Imp_Cov1)]
      }
      res$fit = rcpp_fit_fun(Imp_Cov1, SampCov,type2=0,lambda=0,pen_vec=0,
                             pen_diff=pen_diff,e_alpha=0,gamma=0,rlasso_pen=0,
                             pen_vec1,pen_vec2,
                             dual_pen1,dual_pen2)#,
                           #  estimator2,poly_vec,imp_vec)
    }else if(missing == "fiml" & type == "none"){
      #print(res$optim_fit)
      res$fit = (optFit/nobs)*.5
     # res$fit = rcpp_fit_fun(ImpCov=Imp_Cov,SampCov,
     #                        type2,lambda,pen_vec=0,pen_diff=0)
      #SampCov <- model@implied$cov[[i]]
      #res$fit = rcpp_fit_fun(Imp_Cov1, SampCov,type2=0,lambda=0,pen_vec=0,pen_diff=0)
    }else if(missing=="fiml" & type != "none"){
      if(estimator=="ULS"){
        imp_vec = Imp_Cov[lower.tri(Imp_Cov)]
      }

      res$fit = rcpp_fit_fun(ImpCov=Imp_Cov,SampCov,
                             type2,lambda,pen_vec=0,
                             pen_diff=0,e_alpha=0,gamma=0,rlasso_pen=0,
                             pen_vec1,pen_vec2,
                             dual_pen1,dual_pen2)
                            # estimator2,poly_vec,imp_vec)

    }
    SampCov2 <- SampCov
    SampCov <- model@SampleStats@cov[][[1]]
    res$SampCov = SampCov
    res$SampCov2 <- SampCov2

    res$data <- as.data.frame(model@Data@X)

    res$coefficients <- round(pars.df,round)
    res$nvar = nvar
    res$N = nobs
    res$nfac = nfac

#    if(model@Fit@converged == FALSE){
#      res$baseline.chisq = NA
#      res$baseline.df = NA
#    }else{
#      res$baseline.chisq = lav.fits["baseline.chisq"]
#      res$baseline.df = lav.fits["baseline.df"]
#    }

    #res$grad <- grad(res$par.ret)

   # res$hess <- hess(as.numeric(res$par.ret))
    if(any(is.na(defined_params_vals))==FALSE){

      rettt = rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
       A_est <- rettt$A_est22
       S_est <- rettt$S_est22
       #A_new <- A

       #S_new <- S

     ppars <- defined_params_vals$pars.mult

     #mediation_vals



    parT <- lavaan::parTable(model)

    par.labels <- parT$free > 0 & parT$label != ""
    val <- rep(NA,nrow(parT[par.labels,]))

    for(i in 1:nrow(parT[par.labels,])){
      pp = parT[par.labels,]
      par <- pp[i,"label"]
      par.num <- pp[i,"free"]


        if(any(A==par.num & par.num != 0)){
        val[i] <- round(A_est[A==par.num],round)
        }else if(any(S==par.num & par.num != 0)){
        val[i] <- round(S_est[S==par.num],round)
        }
    }

    labs <- parT$label[par.labels]

    for(i in length(labs):1){
      for(j in 1:length(ppars)){
        ppars[j] <- gsub(labs[i],val[i],ppars[j])
      }
    }

   return.vals <- rep(NA,length(ppars))
   for(i in 1:length(ppars)){
     return.vals[i] <- eval(parse(text=ppars[i]))
   }


   ppars <- defined_params_vals$pars.mult

   for(i in 1:length(ppars)){
     first <- paste("=",return.vals[i])
     ppars[i] <- paste(ppars[i],first)
   }

   res$defined_params <- ppars
   res$defined_params_vals <- return.vals


    }


    if(lambda > 0){
      res$pars_pen <- pars_pen
    }else{
      res$pars_pen <- NULL
    }

    res$mean <- mats$mean

    res$lav.model <- model

    if(res$convergence != 0){
      warning("WARNING: Model did not converge! It is recommended to try multi_optim()")
    }

    res$call <- match.call()
    class(res) <- "regsem"
    return(res)
}
