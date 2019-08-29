#'
#'
#' Penalized model syntax.
#'
#' This function create a lavaan model syntax with paths corresponding to paremeters penalized to 0 removed.
#' @param model lavaan output object.
#' @param nm names(regsemOutput$coefficients).
#' @param pars_pen a vector of numbers corresponding to paths to be removed (same sequence as regsemOutput$coefficients).
#' @return new.mod new model in lavaan syntax.
#' @export


pen_mod<-function(model,nm=NULL,pars_pen=NULL){
  #model: lavaan output object
  #nm: names(regsem$coefficients)
  #pars_pen: parameters to be set to 0 (same sequence as regsem$coefficients)

  parT<-parTable(model)
  id<-parT$id[parT$user==1]
  l<-parT$lhs;r<-parT$rhs
  parl<-l;parl[parT$op=="~"]<-r[parT$op=="~"]
  parr<-r;parr[parT$op=="~"]<-l[parT$op=="~"]

  if (length(pars_pen)==0){
    #warning("no path is removed")
    id.l<-length(id)
    LHS<-parT$lhs;OP<-parT$op;RHS<-parT$rhs
    for (i in 1:id.l){
      #if (parT$free[i]==0){
      RHS[i]<-paste(parT$ustart[i],"*",RHS[i])
      #}
    }
    lhs<-LHS[1];op<-OP[1];rhs<-RHS[1]
    new.list<-list()
    l=1
    for (i in 2:id.l){
      if (OP[i]==op){
        if(LHS[i]==lhs){
          rhs<-append(rhs,RHS[i])
        }else{
          new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
          l=l+1
          lhs<-LHS[i];rhs<-RHS[i]
        }
      }else{
        new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
        l=l+1
        op<-OP[i];lhs<-LHS[i];rhs<-RHS[i]
      }
    }
    new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
    new.mod<-paste(new.list,collapse="\n")
    new.mod
  }else{
    nm.pen<-nm[pars_pen]
    #library(stringr)

    ##need to add " ~~ "
    pen<-data.frame(matrix(unlist(stringr::str_split(nm.pen,pattern=" -> ")),ncol=2,byrow=T))
    colnames(pen)<-c('parl','parr')
    #library(plyr)
    #pars_pen2: parameters to be set to 0 (sequence matched with lavann model)
    pars_pen2<-as.numeric(rownames(plyr::match_df(data.frame(parl,parr),pen)))
    ###matched op "~=" and "~"; not sure about "~~"
    dif.id<-setdiff(id,pars_pen2)
    parT<-parT[dif.id,]
    dif.l<-length(dif.id)
    LHS<-parT$lhs;OP<-parT$op;RHS<-parT$rhs
    for (i in 1:dif.l){
      #if (parT$free[i]==0){
      RHS[i]<-paste(parT$ustart[i],"*",RHS[i])
      #}
    }
    lhs<-LHS[1];op<-OP[1];rhs<-RHS[1]
    new.list<-list()
    l=1
    for (i in 2:dif.l){
      if (OP[i]==op){
        if(LHS[i]==lhs){
          rhs<-append(rhs,RHS[i])
        }else{
          new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
          l=l+1
          lhs<-LHS[i];rhs<-RHS[i]
        }
      }else{
        new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
        l=l+1
        op<-OP[i];lhs<-LHS[i];rhs<-RHS[i]
      }
    }
    new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
    new.mod<-paste(new.list,collapse="\n")
    new.mod
  }
}
