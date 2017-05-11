#'
#'
#' Function to performed exploratory mediation with categorical variables
#'
#' data Name of the dataset
#' iv Name of independent variable
#' mediators Name of mediators
#' dv Name of dependent variable
#' type What type of penalty. Options include lasso, ridge, and enet.
#'



# example

College1 = College[which(College$Private=="Yes"),]
Data = data.frame(scale(College1[c(3,4,9:12,15,17)]))
#lavaan model with all mediators
model1 <-
  ' # direct effect (c)
Enroll ~ c*Accept
# mediators
Outstate ~ a1*Accept
Room.Board ~ a2*Accept
Books ~ a3*Accept
Personal ~ a4*Accept
S.F.Ratio ~ a5*Accept
Expend ~ a6*Accept
Enroll ~ b1*Outstate + b2*Room.Board + b3*Books + b4*Personal + b5*S.F.Ratio + b6*Expend
# indirect effects (a*b)
a1b1 := a1*b1
a2b2 := a2*b2
a3b3 := a3*b3
a4b4 := a4*b4
a5b5 := a5*b5
a6b6 := a6*b6
# total effect (c_prime)
total := c + (a1*b1) + (a2*b2) + (a3*b3) + (a4*b4) + (a5*b5) + (a6*b6)
#Enroll~~0.5*Enroll
'
#p-value approach using delta method standard errors
fit.delta = sem(model1,data=Data,fixed.x=T)
summary(fit.delta)





library(ISLR)

iv <- "Accept"
dv <- "Enroll"
mediators <- c("Outstate","Room.Board","Books","Personal","S.F.Ratio","Expend")

xmed_cat <- function(data,iv,mediators,dv,type="lasso"){
  library(glmnet) # remove
  res <- list()

  data <- Data

  if(type=="lasso"){
    alpha=1
  }else if(type=="ridge"){
    alpha=0
  }else if(type=="enet"){
    alpha=0.5
  }


  iv.mat <- as.matrix(Data[,iv])
  mediators.mat <- as.matrix(data[,mediators])
  dv.mat <- as.matrix(data[,dv])

  if(sum(is.na(iv.mat)) > 0 |
     sum(is.na(mediators.mat)) > 0 |
     sum(is.na(dv.mat)) > 0){
    stop("Missing values are not allowed")
  }

  if(is.numeric(dv.mat)){
    dv.class="gaussian"
  }else if(is.integer(dv.mat)){
    dv.class = "gaussian"
  }else{
    dv.class == "binomial"
  }

  #b's
  b.cv.lasso = cv.glmnet(mediators.mat,dv.mat,alpha=alpha,family=dv.class,
                         penalty.factor=c(rep(1,ncol(data)-2),0))

  b.fit.lasso = glmnet(mediators.mat, dv.mat, alpha=alpha,
                       penalty.factor=c(rep(1,ncol(data)-2),0))

  b.coefs = coef(b.fit.lasso, s=b.cv.lasso$lambda.min)[-1,1]
  #b.coefs = b.coefs[-length(b.coefs)]  ??????? Why -- Need to include?
  res$b.coefs <- round(b.coefs,3)

  # a

  a.cv.lasso = a.fit.lasso = vector("list",ncol(mediators.mat))
  for(i in 1:ncol(mediators.mat)){
    if(is.numeric(mediators.mat[,i])){
      med.class="gaussian"
    }else if(is.integer(mediators.mat[,i])){
      med.class = "gaussian"
    }else{
      med.class == "binomial"
    }
    a.cv.lasso[[i]] = cv.glmnet(as.matrix(cbind(rnorm(nrow(data),1,0.0001),iv.mat)), mediators.mat[,i],
                                alpha=alpha, family=med.class,intercept=F,penalty.factor=c(0,1))
    a.fit.lasso[[i]] = glmnet(as.matrix(cbind(rnorm(nrow(data),1,0.0001),iv.mat)), mediators.mat[,i],
                              alpha=alpha, family=med.class,intercept=F,penalty.factor=c(0,1))
  }

  a.coefs = numeric(length(b.coefs))
  for(i in 1:length(a.coefs)){
    if(!is.null(a.cv.lasso[[i]])){
      a.coefs[i] = coef(a.fit.lasso[[i]], s=a.cv.lasso[[i]]$lambda.min)[-1,1][2]
    }
  }
  names(a.coefs) = mediators
  res$a.coefs <- round(a.coefs,3)

  # indirect

  res$indirect = round(a.coefs * b.coefs,3)


  res$important <- names(res$indirect[res$indirect != 0])

  # add class for summary function

  # return list
  res
}


xmed_cat(data,iv,mediators,dv)
