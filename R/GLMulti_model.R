# predict scores glmulti

predict.scores.glmulti=function(Y,discretspace,
                                map,formula="~I(F1*F1)+I(F2*F2)+F1*F2"
                                ,pred.na=FALSE){
  ## map is a data.frame with F1 and F2 obtained after DR on the explained data Y
  notespr= nbconsos=matrix(0,nrow(discretspace),ncol(discretspace))
  regs=datas=vector("list",ncol(Y))
  pred.conso=preference=matrix(0,nrow = nrow(discretspace),ncol=ncol(Y))
  nb.NA=vector("list",ncol(Y))
  pos.NA=vector("list",ncol(Y))
  dt=cbind.data.frame(rep(1,nrow(map)),map)
  colnames(dt)[1]="y"
  modele=as.formula(paste("y",formula))
  m0=lm(modele,data=dt)
  X0=model.matrix(m0)
  X0=X0[,-1]
  p=ncol(X0)
  cn=paste("x",1:p,sep="")
  colnames(X0)=cn

  dt=cbind.data.frame(rep(1,nrow(discretspace)),discretspace)
  colnames(dt)[1]="y"
  m.pred=lm(modele,data=dt)
  Xpred=model.matrix(m.pred)
  Xpred=Xpred[,-1]
  colnames(Xpred)=cn
  Xpred=data.frame(Xpred)
  ## Firts we preform all regressions
  for(j in 1:ncol(Y)){

    print(j)
    temp=cbind.data.frame(Y[,j],X0)

    colnames(temp)[1]="y"

    m1=c()
    m1=glmulti(y = "y",xr = cn,data=temp,level=1,method = "h",
               fitfunction=lm,plotty=F)

    regs[[j]]=m1@objects[[1]]
    rm(list=c("m1","temp"))
    print(regs[[j]])
    pred.conso[,j]=predict(regs[[j]],newdata=Xpred)

    if (pred.na==TRUE) {
      x=pred.conso[,j]
      x[x<0]=NA
      x[x>10]=NA
      pred.conso[,j]=x
      x=as.data.frame(x)
      nb.NA[[j]] <- apply(x,2,function(a) sum(is.na(a)))
      pos.NA[[j]]=which(is.na(x))
      occur.NA <- unlist(pos.NA)
      occur.NA=as.vector(occur.NA)
      occur.NA=as.data.frame(table(occur.NA))

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

