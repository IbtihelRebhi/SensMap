######  Distances before and after smoothing  when PCA

StabMap<-function(Y,X,S,n,axis=c(1,2),formula_lm,
                  formula_gam,dimredumethod=1,
                  pred.na=FALSE,
                  nbpoints=50) {

      hedotrain <- hedotest <- list()
      dprob_lm <- dprob_lm_loess <- 0
      dprob_gam <- dprob_gam_loess <- 0
      dprob_glmulti <- dprob_glmulti_loess <- 0
      dprob_bayes <- dprob_bayes_loess <- 0


      for (i in (1:n)) {

        y1=sample(1:ncol(Y),size=round(ncol(Y)/2))
        hedotrain[[i]]=as.data.frame(Y[,y1]) # sapmle 1
        hedotest[[i]]=as.data.frame(Y[,setdiff(1:ncol(Y),y1)]) # complementary sample

        if(dimredumethod==1) #PCA on Y
        {
          map<-map.with.pca(X = X,axis = axis)
          map<-cbind.data.frame(map$F1,map$F2)
          colnames(map)=c("F1","F2")
        }

        if(dimredumethod==2) # MFA
        {
          map<-map.with.mfa(X = X,Y = Y,axis = axis)
          map<-cbind.data.frame(map$F1,map$F2)
          colnames(map)=c("F1","F2")

        }

        if(dimredumethod==3) #Canonical Analysis
        {
          map<-map.with.ca(X=X,S=S,Y=Y)
          map<-cbind.data.frame(map$F1,map$F2)
          colnames(map)=c("F1","F2")

        }

        discretspace=discrete.function(map = map)

        ## LM before loess
        reg1<-predict.scores.lm(Y = hedotrain[[i]],formula = formula_lm,discretspace = discretspace,map = map)
        reg2<-predict.scores.lm(Y = hedotest[[i]],formula =  formula_lm,discretspace = discretspace,map = map)
        z.lm1=rowMeans(reg1$pred.conso)
        z.lm2=rowMeans(reg2$pred.conso)
        p.lm1=100*rowMeans(reg1$preference)
        p.lm2=100*rowMeans(reg2$preference)
        graph.surfconso1=as.image(Z=p.lm1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso2=as.image(Z=p.lm2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso1=as.image(Z=z.lm1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso2=as.image(Z=z.lm2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_lm[[i]]=sum((graph.surfconso1$z-graph.surfconso2$z)^2)/(nbpoints^2)


        ## LM after loess

        reg1.loess<-denoising.loess.global(Y = hedotrain[[i]],X,formula = formula_lm,dimredumethod=1,
                                           predmodel=1,discretspace = discretspace)
        reg2.loess<-denoising.loess.global(Y = hedotest[[i]],X,formula = formula_lm,dimredumethod=1,
                                           predmodel=1,discretspace = discretspace)
        z1.loess=reg1.loess$pred.conso
        p1.loess=reg1.loess$preference
        z2.loess=reg2.loess$pred.conso
        p2.loess=reg2.loess$preference
        graph.surfconso.l1=as.image(Z=p1.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso.l2=as.image(Z=p2.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l1=as.image(Z=z1.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l2=as.image(Z=z2.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_lm_loess[[i]]=sum((graph.surfconso.l1$z-graph.surfconso.l2$z)^2)/(nbpoints^2)




        ##GAM before loess

        reg1gam<-predict.scores.gam(Y = hedotrain[[i]],formula = formula_gam,discretspace = discretspace,map = map)
        reg2gam<-predict.scores.gam(Y = hedotest[[i]],formula = formula_gam,discretspace = discretspace,map = map)
        z.gam1=rowMeans(reg1gam$pred.conso)
        z.gam2=rowMeans(reg2gam$pred.conso)
        p.gam1=100*rowMeans(reg1gam$preference)
        p.gam2=100*rowMeans(reg2gam$preference)
        graph.surfconso1g=as.image(Z=p.gam1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso2g=as.image(Z=p.gam2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso1g=as.image(Z=z.gam1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso2g=as.image(Z=z.gam2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_gam[[i]]=sum((graph.surfconso1g$z-graph.surfconso2g$z)^2)/(nbpoints^2)


        ## GAM after loess
        reg1gam.loess<-denoising.loess.global(Y = hedotrain[[i]],Y,formula = formula_gam,dimredumethod=1,
                                              predmodel=2,discretspace = discretspace)
        reg2gam.loess<-denoising.loess.global(Y = hedotest[[i]],Y,formula = formula_gam,dimredumethod=1,
                                              predmodel=2,discretspace = discretspace)
        z1g.loess=reg1gam.loess$pred.conso
        p1g.loess=reg1gam.loess$preference
        z2g.loess=reg2gam.loess$pred.conso
        p2g.loess=reg2gam.loess$preference

        graph.surfconso.l.g1=as.image(Z=p1g.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso.l.g2=as.image(Z=p2g.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.g1=as.image(Z=z1g.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.g2=as.image(Z=z2g.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_gam_loess[[i]]=sum((graph.surfconso.l.g1$z-graph.surfconso.l.g2$z)^2)/(nbpoints^2)



        ##glmulti before loess

        reg1glm<-predict.scores.glmulti(Y = hedotrain[[i]],formula = formula_lm,discretspace = discretspace,map = map)
        reg2glm<-predict.scores.glmulti(Y = hedotest[[i]],formula = formula_lm,discretspace = discretspace,map = map)
        z.glm1=rowMeans(reg1glm$pred.conso)
        z.glm2=rowMeans(reg2glm$pred.conso)
        p.glm1=100*rowMeans(reg1glm$preference)
        p.glm2=100*rowMeans(reg2glm$preference)
        graph.surfconso1glm=as.image(Z=p.glm1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso2glm=as.image(Z=p.glm2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso1glm=as.image(Z=z.glm1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso2glm=as.image(Z=z.glm2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_glmulti[[i]]=sum((graph.surfconso1glm$z-graph.surfconso2glm$z)^2)/(nbpoints^2)


        ## glmulti after loess
        reg1glm.loess<-denoising.loess.global(Y = hedotrain[[i]],Y,formula = formula_lm,dimredumethod=1,
                                              predmodel=3,discretspace = discretspace)
        reg2glm.loess<-denoising.loess.global(Y = hedotest[[i]],Y,formula = formula_lm,dimredumethod=1,
                                              predmodel=3,discretspace = discretspace)
        z1glm.loess=reg1glm.loess$pred.conso
        p1glm.loess=reg1glm.loess$preference
        z2glm.loess=reg2glm.loess$pred.conso
        p2glm.loess=reg2glm.loess$preference

        graph.surfconso.l.glm1=as.image(Z=p1glm.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso.l.glm2=as.image(Z=p2glm.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.glm1=as.image(Z=z1glm.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.glm2=as.image(Z=z2glm.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_glmulti_loess[[i]]=sum((graph.surfconso.l.glm1$z-graph.surfconso.l.glm2$z)^2)/(nbpoints^2)



        ## Bayes before loess

        reg1bayes<-predict.scores.bayes(Y = hedotrain[[i]],discretspace= discretspace,map=map,formula=formula_lm)
        reg2bayes<-predict.scores.bayes(Y = hedotest[[i]],discretspace= discretspace,map=map,formula=formula_lm)
        z.bayes1=rowMeans(reg1bayes$pred.conso)
        z.bayes2=rowMeans(reg2bayes$pred.conso)
        p.bayes1=100*rowMeans(reg1bayes$preference)
        p.bayes2=100*rowMeans(reg2bayes$preference)
        graph.surfconso.b1=as.image(Z=p.bayes1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso.b2=as.image(Z=p.bayes2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.b1=as.image(Z=z.bayes1,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.b2=as.image(Z=z.bayes2,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_bayes[[i]]=sum((graph.surfconso.b1$z-graph.surfconso.b2$z)^2)/(nbpoints^2)


        ## Bayes after loess
        reg1bayes.loess<-denoising.loess.global(Y = hedotrain[[i]],X,dimredumethod=1,
                                                predmodel=4,discretspace = discretspace,
                                                formula=formula_lm)
        reg2bayes.loess<-denoising.loess.global(Y = hedotest[[i]],X,dimredumethod=1,
                                                predmodel=4,discretspace = discretspace,
                                                formula=formula_lm)
        z1b.loess=reg1bayes.loess$pred.conso
        p1b.loess=reg1bayes.loess$preference
        z2b.loess=reg2bayes.loess$pred.conso
        p2b.loess=reg2bayes.loess$preference

        graph.surfconso.l.b1=as.image(Z=p1b.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.surfconso.l.b2=as.image(Z=p2b.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.b1=as.image(Z=z1b.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso
        graph.predconso.l.b2=as.image(Z=z2b.loess,x=discretspace,nrow=nbpoints,ncol=nbpoints) # Surface d'un seul conso

        dprob_bayes_loess[[i]]=sum((graph.surfconso.l.b1$z-graph.surfconso.l.b2$z)^2)/(nbpoints^2)

      }


      result <- rbind(dprob_lm,dprob_lm_loess,
                      dprob_gam,dprob_gam_loess,
                      dprob_glmulti,dprob_glmulti_loess,
                      dprob_bayes,dprob_bayes_loess)

      rownames(result) <- c("dprob_lm","dprob_lm_loess",
                            "dprob_gam","dprob_gam_loess",
                            "dprob_glmulti","dprob_glmulti_loess",
                            "dprob_bayes","dprob_bayes_loess")


      return(result)


    }
