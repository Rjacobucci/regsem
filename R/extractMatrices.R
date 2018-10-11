#'
#'
#' This function extracts RAM matrices from a lavaan object.
#' @param model Lavaan model object.
#' @return The RAM matrices from \code{model}.
#' @keywords extract
#' @import lavaan
#' @export
#' @examples
#'
#' library(lavaan)
#' data(HolzingerSwineford1939)
#' HS.model <- ' visual =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed =~ x7 + x8 + x9 '
#' mod <- cfa(HS.model, data=HolzingerSwineford1939)
#' mats = extractMatrices(mod)


extractMatrices <- function(model){

matrices <-  list()


#library(Matrix)
#library(lavaan)
pars <- lavaan::parameterestimates(model)
parT <- lavaan::parTable(model)

#parT = parT[parT$exo != 1,]
#pars = pars[parT$exo != 1,]


nfac.hold1 <- pars[pars$op == "=~",]
nfac1 <- length(unique(nfac.hold1$lhs))

mean =FALSE

if(any(parT$op == "~1")){
  parTT = parT[parT$op == "~1",]
  if(any(parTT$free > 0)){
    mean = TRUE
  }
}



nfac2 = nfac1 + ifelse(mean==TRUE,1,0)



# check for groups
#model@pta$ngroups


nvar = model@pta$nvar[[1]][1]

A_init <- matrix(0, nrow = nvar + nfac2, ncol = nvar + nfac2)

#unique(nfac.hold1$lhs)
name <- unique(pars$lhs)

name.vars <- model@pta$vnames$ov[[1]]


name.factors <- model@pta$vnames$lv[[1]]

if(identical(name.factors, character(0))){
  name.factors2=NA
}else{
  name.factors2 = name.factors
}

if(length(name.factors)!=0){
  if(mean==TRUE){
    colnames(A_init) <- c(name.vars,"1",name.factors)
    rownames(A_init) <- c(name.vars,"1",name.factors)
  }else{
    colnames(A_init) <- c(name.vars,name.factors)
    rownames(A_init) <- c(name.vars,name.factors)
  }
}else{
  if(mean==TRUE){
    colnames(A_init) <- c(name.vars,"1")
    rownames(A_init) <- c(name.vars,"1")
  }else{
    colnames(A_init) <- c(name.vars)
    rownames(A_init) <- c(name.vars)
  }
}

A <- A_init

A.parT = parT[parT$op == "=~" | parT$op == "~1" | parT$op == "~",]
A.pars = pars[parT$op == "=~" | parT$op == "~1" | parT$op == "~",]


A.parFree <- A.parT[A.parT$free > 0,]


uniq <- unique(A.parFree[,"label"])

uniq2 <- uniq[table(A.parFree[,"label"]) == 1]
pars.align.A = matrix(NA,length(A.parFree$free),2)

A.parFree33 = A.parFree
if(length(uniq2)>0){
  for(i in 1:length(uniq2)){
    A.parFree[A.parFree[,"label"] == uniq2[i],"label"] <- ""
  }
}



# any equality?
if(any(duplicated(A.parFree$label[A.parFree$label != ""]) == T)){
  labels = unique(A.parFree$label[A.parFree$label != ""])
  for(i in 1:length(labels)){
    equals = A.parFree$free[A.parFree$label == labels[i]]
    min.equal = min(equals)
    max.equal = max(equals)
    A.parFree$free[A.parFree$label == labels[i]] <- min.equal

  #  dec = max.equal - min.equal
   # A.parFree$free[A.parFree$label == labels[i]] <- A.parFree$free[A.parFree$label == labels[i]] - dec
  }
  #equals = A.parFree$free[A.parFree$label == labels]

}


A.parFree[,"free"] <- rank(A.parFree[,"free"],ties.method="min")


A.parFree2 <- A.parFree

loadings <- NULL
regressions <- NULL

if(nrow(A.parFree2) > 0){
for(i in 1:nrow(A.parFree2)){
  if(A.parFree2$op[i] == "=~"){
  colNum <- which(A.parFree2$lhs[i] == colnames(A))
  rowNum <- which(A.parFree2$rhs[i] == rownames(A))
  A[rowNum,colNum] = A.parFree2[i,"free"]
  loadings = c(loadings,A.parFree2[i,"free"])
  }else if(A.parFree2$op[i] == "~1"){
    A[which(rownames(A)==A.parFree2$lhs[i]),which(colnames(A) == "1")] = A.parFree2[i,"free"]
  }else if(A.parFree2$op[i] == "~"){
    colNum <- which(A.parFree2$rhs[i] == colnames(A))
    rowNum <- which(A.parFree2$lhs[i] == rownames(A))
    A[rowNum,colNum] = A.parFree2[i,"free"]
    regressions = c(regressions,A.parFree2[i,"free"])
  }
  pars.align.A[i,] = c(i,A.parFree33[i,"label"])
}
}else{
  A = A
}






# A of free parameters
A_fixed <- A > 10000

parA_fixed = A.parT[A.parT$free == 0,]

if(nrow(parA_fixed) > 0){
  for(gg in 1:nrow(parA_fixed)){
    if(parA_fixed[gg,"op"]=="=~"){
      coll = which(parA_fixed$lhs[gg] == colnames(A_fixed))
      roww = which(parA_fixed$rhs[gg] == rownames(A_fixed))
    }else if(parA_fixed[gg,"op"]=="~"){
      coll = which(parA_fixed$rhs[gg] == colnames(A_fixed))
      roww = which(parA_fixed$lhs[gg] == rownames(A_fixed))
    }else if(parA_fixed[gg,"op"]=="~1"){
      coll = which("1" == colnames(A_fixed))
      roww = which(parA_fixed$lhs[gg] == rownames(A_fixed))
    }
    A_fixed[roww,coll] = T
  }
}




# A_est
A_est = A_init
if(nrow(A.pars) > 0){
for(i in 1:nrow(A.pars)){
  if(A.pars$op[i] == "=~"){
    colNum <- which(A.pars$lhs[i] == colnames(A_est))
    rowNum <- which(A.pars$rhs[i] == rownames(A_est))
    A_est[rowNum,colNum] = A.pars$est[i]
  }else if(A.pars$op[i] == "~1"){
    A_est[which(rownames(A_est)==A.pars$lhs[i]),which(colnames(A_est) == "1")] = A.pars$est[i]
  }else if(A.pars$op[i] == "~"){
    colNum <- which(A.pars$rhs[i] == colnames(A_est))
    rowNum <- which(A.pars$lhs[i] == rownames(A_est))
    A_est[rowNum,colNum] = A.pars$est[i]
}
}
}else{
  A_est = A_est
}
# create F
F <- A_init[1:(nvar + ifelse(mean==TRUE,1,0)),]

varNames <- model@pta$vnames$ov.model[[1]]
f.names <- colnames(A)[1:(length(varNames) + ifelse(mean==TRUE,1,0))]

for(ii in 1:length(f.names)){
  #Numm = which(colnames(F)[ii] == f.names[ii])
  F[ii,ii] = 1
}


# create S

covars <- pars[pars$op == "~~",]
covarT <- parT[parT$op == "~~",]
S_init <- matrix(0, nrow = nvar + nfac2, ncol = nvar + nfac2)

colnames(S_init) <- colnames(A); rownames(S_init) <- rownames(A)


# S_est -- equality doesn't effect
S_est <- as.matrix(S_init)
for(jj in 1:nrow(covars)){
  col1 = which(colnames(S_est) == covars$lhs[jj])
  row1 = which(rownames(S_est) == covars$rhs[jj])
  S_est[row1,col1] = covars$est[jj]
}



if(any(S_est[lower.tri(S_est)] != 0)){
  S_est <- S_est + t(S_est) - diag(diag(S_est))
}



S_est[rownames(S_est) == "1",colnames(S_est) == "1"] <- 1

# S
S <- S_init

covarT.free <- covarT[covarT$free > 0, ]
pars.align.S = matrix(NA,length(covarT$free),2)

# any equality?
if(any(duplicated(covarT.free$label[covarT.free$label != ""]) == T)){
  labels = unique(covarT.free$label[covarT.free$label != ""])
  for(i in 1:length(labels)){
    equals = covarT.free$free[covarT.free$label == labels[i]]
    min.equal = min(equals)
    max.equal = max(equals)
    covarT.free$free[covarT.free$label == labels[i]] <- min.equal

    dec = max.equal - min.equal
    #covarT.free$free[covarT.free$label != labels[i]] <- covarT.free$free[covarT.free$label != labels[i]] - dec
  }
}


for(jj in 1:nrow(covarT.free)){
  col1 = which(colnames(S) == covarT.free$lhs[jj])
  row1 = which(rownames(S) == covarT.free$rhs[jj])
  S[row1,col1] = covarT.free$free[jj]
}

if(any(S[lower.tri(S)] != 0)){
  S <- S + t(S) - diag(diag(S))
}

ss= S[S != 0]
for(i in 1:length(unique(ss))){
  val = ss[ss == min(ss)]
  S[S == val[1]] <- i
  ss = ss[ss != val[1]]
}


if(sum(S >0 ) > 0){
dec2 = max(A) - min(S[S != 0]) + 1
S[S != 0] = S[S != 0] + dec2
}else{
  S = S
}

for(i in 1:length(covarT.free$free)){
  pars.align.S[i,] = c(sort(unique(S[S>0]))[i],covarT.free$label[i])
}

# S_fixed
#S_fixed <- S_init
S_fixed <- S_init > 10000

covarT.fixed <- covarT[covarT$free == 0, ]

if(nrow(covarT.fixed) > 0){
for(jj in 1:nrow(covarT.fixed)){
  col1 = which(colnames(S_fixed) == covarT.fixed$lhs[jj])
  row1 = which(rownames(S_fixed) == covarT.fixed$rhs[jj])
  S_fixed[row1,col1] = TRUE
}
}


if(any(S_fixed - diag(S_fixed) == 0 )){
  S_fixed = S_fixed + t(S_fixed) -diag(diag(S_fixed))
  S_fixed = S_fixed == 1
}else{
  S_fixed = S_fixed
}

S_fixed[rownames(S_est) == "1",colnames(S_est) == "1"] <- TRUE




pars <- rep(NA,length(unique(c(A[A !=0],S[S !=0]))))
#pars <- rep(NA,max(max(A),max(S)))
count = 0
for(tt in 1:max(max(A),max(S))){
  if(any(A == tt)==TRUE){
   count = count+1
    pars[count] = A_est[A==tt][1]
  }
  else if(any(S == tt)==TRUE){
   count = count +1
    pars[count] = S_est[S==tt][1]
  }
}



count = 0
for(i in 1:max(max(A),max(S))){
  if(any(A == i)){
    count = count + 1
    pos = which(A == i,arr.ind=T)
    one = colnames(A)[pos[1,2]]
    two = rownames(A)[pos[1,1]]
    names(pars)[count] = paste(one,"->",two)
  }else if(any(S==i)){
   count = count + 1
    pos = which(S == i,arr.ind=T)

#    if(nrow(pos) == 1){
     one = colnames(S)[pos[1,2]]
     two = rownames(S)[pos[1,1]]
#    }else if(nrow(pos) > 1){
#     one = colnames(S)[pos[1,2]]
#     two = rownames(S)[pos[1,1]]
#    }
    names(pars)[count] = paste(one,"~~",two)
  }
}


# add parameter labels

#var.labs <- unique(parT[parT[,"label"] != "","label"])

#if(length(var.labs)>0){
#  for(i in 1:length(var.labs)){



#    pars
#  }
#}


for(i in 1:length(unique(A[A>0]))){
  ord <- sort(unique(A[A>0]))
  A[A==ord[i]] <- i
}

for(i in 1:length(unique(S[S>0]))){
  ord <- sort(unique(S[S>0]))
  S[S==ord[i]] <- max(A) + i
}


# get mediation parameters
# only record arguments

if(any(parT$op == ":=")){

med.pars <- parT[parT$op == ":=",]
args <- med.pars$rhs
args.labs <- med.pars$lhs

sub09 = (parT$label != "")
labels <- parT$label[sub09]
labs.inds <- labels %in% args.labs
pars.labels <- labels[labs.inds]

# now need parameter numbers

labs.single <- labels[labs.inds==FALSE]

par.extra <- rep(NA,length(labs.single))
for(i in 1:length(labs.single)){
  par.extra[i] <- parT[parT$label == labs.single[i],]$free
}

# now match parameter numbers with labels
# for each in args

#for(i in 1:length(args)){
#  for(j in 1:length(labs.single)){
#    args[i] <- gsub(labs.single[j],par.extra[j],args[i])
#  }
#}



mediation <- list()

mediation$pars.mult <- args
mediation$pars.labs <- args.labs

}else{
  mediation <- NA
}




# return Matrices
matrices$A <- A
matrices$A_est <- A_est
matrices$A_fixed <- A_fixed
matrices$S <- S
matrices$S_est <- S_est
matrices$S_fixed <- S_fixed
matrices$F <- F
matrices$parameters <- round(pars,3)
matrices$mean <- mean
matrices$mediation <- mediation
matrices$name.factors <- name.factors2
matrices$loadings <- loadings
matrices$regressions <- regressions
matrices$pars.align <- rbind(pars.align.A,pars.align.S)
matrices

}
