# prediction with gam

predict.scores.gam=function(Y,discretspace,map
                            ,formula="~s(F1,k=2)+s(F2,k=2)"
                            ,pred.na=FALSE){
  ## map is a data.frame with F1 and F2 obtained after DR on the explained data Y
  notespr= nbconsos=matrix(0,nrow(discretspace),ncol(discretspace))
  regs=vector("list",ncol(Y))
  # preference=array(0,dim=c(nrow(discretspace),ncol(discretspace),ncol(X)))
  pred.conso=preference=matrix(0,nrow(discretspace),ncol(Y))
  ## Firts we preform all regressions
  nb.NA=vector("list",ncol(Y))
  pos.NA=vector("list",ncol(Y))
  nbconsos=c()
  for(j in 1:ncol(Y)){
    map.reg=cbind.data.frame(Y[,j],map)
    colnames(map.reg)[1]="Conso"
    modele=as.formula(paste("Conso",formula))
    regs[[j]]=gam(modele,data=map.reg)
    pred.conso[,j]=predict(regs[[j]],newdata=discretspace)

    if (pred.na==TRUE) {
      x=pred.conso[,j]
      x[x<0]=NA
      x[x>10]=NA
      pred.conso[,j]=x
      x=as.data.frame(x)
      nb.NA[[j]] <- apply(x,2,function(a) sum(is.na(a)))# nbre de NA pour chaque conso
      pos.NA[[j]]=which(is.na(x))# le point qui contient NA pour chaque consom
      occur.NA <- unlist(pos.NA)
      occur.NA=as.vector(occur.NA)
      occur.NA=as.data.frame(table(occur.NA))# nbre de NA en chaque point du plan pour tous les consos

    }
    else {
      nb.NA=0
      pos.NA=0
      occur.NA=0
    }
    preference[,j]=(pred.conso[,j]> mean(Y[,j]))
  }
  return(list(regression=regs,pred.conso=pred.conso,preference=preference,nb.NA=nb.NA,
              pos.NA=pos.NA, occur.NA=occur.NA))

}

