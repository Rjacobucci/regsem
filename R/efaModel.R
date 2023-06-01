#'
#'
#' Generates an EFA model to be used by lavaan and regsem
#' Function created by Florian Scharf for the paper
#' Should regularization replace simple structure rotation in
#' Exploratory Factor Analysis -- Scharf & Nestler (in press at SEM)
#' @param nFactors Number of latent factors to generate.
#' @param variables Names of variables to be used as indicators
#' @return model Full EFA model parameters.
#' @keywords fa factor analysis
#' @export
#' @examples
#' \dontrun{
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#'# Note to find number of factors, recommended to use
#'# fa.parallel() from the psych package
#'# using the wrong number of factors can distort the results
#'mod = efaModel(3, colnames(HS))
#'
#'semFit = sem(mod, data = HS, int.ov.free = FALSE, int.lv.free = FALSE,
#'             std.lv = TRUE, std.ov = TRUE, auto.fix.single = FALSE, se = "none")
#'
#'# note it requires smaller penalties than other applications
#'reg.out2 = cv_regsem(model = semFit, pars_pen = "loadings",
#'                     mult.start = TRUE, multi.iter = 10,
#'                     n.lambda = 100, type = "lasso", jump = 10^-5, lambda.start = 0.001)
#'reg.out2
#'plot(reg.out2) # note that the solution jumps around -- make sure best fit makes sense
#' }



efaModel = function(nFactors, variables){
  factors = paste("f",1:nFactors,sep = "") #prepare factor labels
  loadings = matrix(paste("lambda_",expand.grid(1:length(variables), 1:nFactors)[,2],"_", expand.grid(1:length(variables), 1:nFactors)[,1], sep = "")
                    ,length(variables),nFactors) # prepare loading matrices

  # don't fix first loading to 1
  loadings[1,] = NA


  #set restricions
  # sum of loadings > 0 per factor (to prefer positive loadings)
  # Note the factors can be sign-inverted any time, it is just for convenience (e.g., plotting)
  # that we prefer positive loadings, this procedure imitates the behaviour of Mplus
  loadingConstraints = paste(apply(loadings[-1,], MARGIN = 2, function(x){
    paste(paste(x,collapse = " + "), " > 0", collapse = "")
  }), collapse = "\n")

  loadings = apply(loadings,MARGIN = 2,FUN = function(x){paste(x,variables,sep = "*")})

  # These lines are only an automatized way to write the model definitions
  # for a large number of variables.
  # In order to make this step more transparent, the model definition
  # is saved in a document called "model.txt". We recommend you refer to this file
  # if you want to see the model specification in detail.
  measurementModels = sapply(1:nFactors,FUN = function(x){paste(factors[x],"=~",paste(loadings[,x],collapse = "+"),collapse = "\n")})
  factorCombinations = expand.grid(factors,factors)
  test = expand.grid(1:nFactors,1:nFactors) #remove redundant covariances
  factorCombinations = factorCombinations[test[,1] <= test[,2],]
  factorCovariances = sapply(1:dim(factorCombinations)[1], function(x){ifelse(factorCombinations[x,1] == factorCombinations[x,2],paste(factorCombinations[x,1],"~~1*",factorCombinations[x,2],sep = ""),paste(factorCombinations[x,1],"~~NA*",factorCombinations[x,2],sep = ""))})

  model = c(measurementModels,loadingConstraints, factorCovariances)
 # writeLines(model, "model.txt")
  return(model)
}
