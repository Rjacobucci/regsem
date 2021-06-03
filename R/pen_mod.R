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
  id<-parT$id[parT$free!=0]

  if (length(pars_pen)==0){#no path to be removed
    #warning("no path is removed")
    id.l<-length(id)
    LHS<-parT$lhs;OP<-parT$op;RHS<-parT$rhs
    for (i in 1:id.l){
      t<-id[i]
      if (OP[t]!="~1"){
        RHS[t]<-paste(parT$ustart[t],"*",RHS[t])
      }
    }
    lhs<-LHS[id[1]];op<-OP[id[1]];rhs<-RHS[id[1]]
    new.list<-list()
    l=1
    for (i in 2:id.l){
      t=id[i]
      if (OP[t]==op){
        if(LHS[t]==lhs){
          rhs<-append(rhs,RHS[t])
        }else{
          new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
          l=l+1
          lhs<-LHS[t];rhs<-RHS[t]
        }
      }else{
        new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
        l=l+1
        op<-OP[t];lhs<-LHS[t];rhs<-RHS[t]
      }
    }
  }else{#exist path to be removed
    nm.pen<-nm[pars_pen]
    var<-nm[grep("->",nm)]
    reg<-unlist(stringr::str_split(var,pattern=" -> "))
    reg<-data.frame(matrix(reg,ncol=2,byrow=T))
    reg<-reg[reg$X1!=1,]
    var<-unique(c(reg$X1,reg$X2))
    #library(stringr)
    #further drop covariance structure and mean structure if a variable is removed from model.
    nm.nonpen<-nm[-pars_pen]
    var.reg<-nm.nonpen[grep("->",nm.nonpen)]
    un.var.reg<-unlist(stringr::str_split(var.reg,pattern=" -> "))
    unpen.reg<-data.frame(matrix(un.var.reg,ncol=2,byrow=T))
    unpen.reg<-unpen.reg[unpen.reg$X1!=1,]
    var.inc<-unique(c(unpen.reg$X1,unpen.reg$X2))
    var.rm<-setdiff(var,var.inc)
    if (length(var.rm)!=0) {
      for (c in 1:length(var.rm)){
        rm<-var.rm[c]
        rm.nm<-nm[grep(rm,nm)]
        nm.pen<-union(nm.pen,rm.nm)
      }
    } #nm.pen updated

    reg.nm<-nm.pen[grep("->",nm.pen)]
    un.reg<-unlist(stringr::str_split(reg.nm,pattern=" -> "))
    n.reg<-length(un.reg)

    cor.nm<-nm.pen[grep("~~",nm.pen)]
    un.cor<-unlist(stringr::str_split(cor.nm,pattern=" ~~ "))
    n.cor<-length(un.cor)


    n.pen<-length(nm.pen)
    if (n.reg==0){
      pen<-data.frame(matrix(un.cor,ncol=2,byrow=T))
    }else if (n.cor==0){
      pen<-data.frame(matrix(un.reg,ncol=2,byrow=T))
    }else{
      pen.reg<-data.frame(matrix(un.reg,ncol=2,byrow=T))
      pen.cor<-data.frame(matrix(un.cor,ncol=2,byrow=T))
      pen<-rbind(pen.reg,pen.cor)
    }
    colnames(pen)<-c('parl','parr')

    l<-parT$lhs;r<-parT$rhs
    parl<-l;parl[parT$op=="~"]<-r[parT$op=="~"] #switch to same sequence as using "->"
    parr<-r;parr[parT$op=="~"]<-l[parT$op=="~"]
    parl[parT$op=="~1"]<-1;parr[parT$op=="~1"]<-l[parT$op=="~1"]
    #library(plyr)
    #pars_pen2: parameters to be set to 0 (sequence matched with lavann model)
    pars_pen2<-as.numeric(rownames(plyr::match_df(data.frame(parl,parr),pen)))

    dif.id<-setdiff(id,pars_pen2)#to be kept
    parT<-parT[dif.id,]
    dif.l<-length(dif.id)
    LHS<-parT$lhs;OP<-parT$op;RHS<-parT$rhs
    for (i in 1:dif.l){
      if (OP[i]!="~1"){
        RHS[i]<-paste(parT$ustart[i],"*",RHS[i])
      }
    }
    lhs<-LHS[1];op<-OP[1];rhs<-RHS[1]
    new.list<-list()
    l=1
    for (i in 2:dif.l){
      if (OP[i]==op & op!="~~"){
        if(LHS[i]==lhs){
          rhs<-append(rhs,RHS[i])
        }else{
          new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
          l=l+1
          lhs<-LHS[i];rhs<-RHS[i]
        }
      }else if (OP[i]!=op & op!="~~"){
        new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
        l=l+1
        op<-OP[i];lhs<-LHS[i];rhs<-RHS[i]
      }else if (op=="~~"){
        new.list[[l]]<-paste(lhs,op,rhs)
        l=l+1
        op<-OP[i];lhs<-LHS[i];rhs<-RHS[i]
      }
    }
  }
  new.list[[l]]<-paste(paste(lhs,op),paste(rhs,collapse="+"))
  new.mod<-paste(new.list,collapse="\n")
  new.mod
}
