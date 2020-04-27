clock2 <- function(predictors, KDE_layer, data, nb_replicates) {
  
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("The function requires 'raster' package to work, please install it.", call. = FALSE)
  }
  if (!requireNamespace("dismo", quietly = TRUE)) {
    stop("The function requires 'dismo' package to work, please install it.", call. = FALSE)
  }
  
  # Preparation of the model, initialize vectors and matrices  
  #------------------------------------------------------------
  eval<-matrix(nrow=15,ncol=nb_replicates,data=NA); colnames(eval)<-seq(1,ncol(eval),1);rownames(eval)<-c("AUC_cv","AUC_self","AUC_all","maxSSS","COR","pCOR","TSS", "valid_test_pres","nb_NA_test","nb_test","moran_resi","pval_resi","sd_moran_resi","ntrees","sampled_longitude")
  stack.pred<-subset(predictors,1);values(stack.pred)<-NA
  testvaluesp<-rep(NA,nrow(unique(data))) ; testvaluesa<-testvaluesp
  
  # Define presence and background data, and associated environmental values 
  #--------------------------------------------------------------------------
  presvals.unique_glob <- raster::extract(predictors, unique(data)) 
  
  background_data <- xyFromCell(KDE_layer, sample(which(!is.na(values(KDE_layer))), dim(presvals.unique_glob)[1], prob=values(KDE_layer)[!is.na(values(KDE_layer))])) # define lat and long of the background data
  colnames(background_data) <- colnames(data)

  pseudoabs_glob <- extract(predictors,background_data) #extract the value of the environment of the background data
  
  for (j in 1:nb_replicates){
    
    # for each replicate, the background data is sampled  
    background_data <- xyFromCell(KDE_layer, sample(which(!is.na(values(KDE_layer))), dim(presvals.unique_glob)[1], prob=values(KDE_layer)[!is.na(values(KDE_layer))]))
    colnames(background_data) <- colnames(data)
    
    random_long <- seq(-180,180,1)
    random_long_opp <- rep(NA,361)
    
    # nombres opposés à chaque longitude pour délimiter la borne supérieure de l'échantillonnage
    for (i in 1:361){
      if (random_long[i] >= 0){
        random_long_opp[random_long > 0] <- random_long[random_long > 0]-180
      } else{
        random_long_opp[random_long <= 0] <- random_long[random_long <= 0]+180
      }
    }
    
    #tirer un nombre entre -180 et 180 pour définir la borne inférieure de l'échantillonnage
    random_long_tirage <- sample(random_long,1)
    random_long_tirage_t <- random_long_opp[random_long==random_long_tirage]
    
    ## échantillonner le groupe test et le groupe training
    presence_tot <- unique(data)
    
    if (random_long_tirage>0){
      presence_training <- subset(presence_tot,presence_tot[,1] <random_long_tirage & presence_tot[,1] >random_long_tirage_t)
      background_training <- subset(background_data,background_data[,1] <random_long_tirage & background_data[,1] >random_long_tirage_t)
      presence_test <- presence_tot[-which(presence_tot[,1] <random_long_tirage & presence_tot[,1] >random_long_tirage_t),]
      background_test <- background_data[-which(background_data[,1] <random_long_tirage & background_data[,1] >random_long_tirage_t),]
      borne_sup <- random_long_tirage
      borne_inf <- random_long_tirage_t
      
    } else {
      presence_training <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] <random_long_tirage_t)
      background_training <- subset(background_data,background_data[,1] >random_long_tirage & background_data[,1] <random_long_tirage_t)
      presence_test <- presence_tot[-which(presence_tot[,1] > random_long_tirage & presence_tot[,1] <random_long_tirage_t),]
      background_test <- background_data[-which(background_data[,1] >random_long_tirage & background_data[,1] <random_long_tirage_t),]
      borne_sup <- random_long_tirage_t
      borne_inf <- random_long_tirage
    }
    
    # plot(subset(predictors,1))
    # points(background_training, pch=20, col="grey")
    # points(presence_training, pch=20)
    # points(background_test, pch=20, col="pink")
    # points(presence_test, pch=20, col="darkred")
    
    training_presences.occ <- presence_training
    training_backg.occ <- background_training
    test_presences.occ <- presence_test
    test_back.occ <- background_test
  
    
    ####### LAUNCH THE MODEL ON THE TRAINING DATA
    #-----------------------------------------------------------
    # BUILD THE MATRIX TO CALIBRATE BRST
    #------------------------------------------------
    ### PRESENCE DATA
    
    presence_data <- training_presences.occ
    presvals.unique <- extract(predictors, presence_data) 
    
    
    ## BACKGROUND DATA
    #-----------------------
    pseudoabs <- extract(predictors,training_backg.occ)
    
    
    ## BRT CALIBRATION
    #-----------------------
    # MATRIX DATA
    id<-0;sdmdata.unique<-0;  id<-c(rep(1,nrow(presvals.unique)),rep(0,nrow(pseudoabs))) 
    FICHIER<-data.frame(cbind(id,rbind(presvals.unique,pseudoabs)))
    
    # Calibration parameters 
    tc=4   # tree complexity
    lr=0.007 # learning rate
    bf=0.75  # bag fraction
    
    # BRT Launch
    model.res<- gbm.step(data=FICHIER, 
                         gbm.x = 2:ncol(FICHIER),
                         gbm.y = 1,
                         family = "bernoulli",
                         tree.complexity = tc,
                         learning.rate = lr,
                         bag.fraction = bf)

    eval[14,j] <- model.res$n.trees
    
    #------------------------------------------------------------
    # BRT OUTPUTS  
    #------------------------------------------------------------
    # prediction stack
    p<-predict(predictors,model.res,n.trees=model.res$gbm.call$best.trees,type="response", na.rm=F)
    p <- mask(p,subset(pred,1))
    stack.pred<-stack(stack.pred,p)
    
    # valeurs des AUC, testp et testa
    testp<-dismo::predict(model.res,data.frame(presvals.unique_glob),n.trees=model.res$gbm.call$best.trees,type="response")
    testa<-dismo::predict(model.res,data.frame(pseudoabs_glob),n.trees=model.res$gbm.call$best.trees,type="response")
    
    eval3<-dismo::evaluate(p=testp,a=testa)
    # récupérer l'AUC 
    eval[3,j]<-eval3@auc
    eval[5,j]<-eval3@cor # point biserial correlation => correlation entre présence/background et valeur de proba
    eval[6,j]<-eval3@pcor
    # récupérer l'AUC 
    eval[1,j]<-model.res$cv.statistics$discrimination.mean
    eval[2,j]<-model.res$self.statistics$discrimination
    
    # mesure TSS
    specificity <- (eval3@TPR/(eval3@TPR+eval3@FPR))
    sensitivity <- (eval3@TNR/(eval3@TNR+eval3@FNR))
    tss <- specificity + sensitivity -1
    eval[7,j]<- mean(tss, na.rm=T)
    
    # mesure maxSSS
    tab<-cbind(eval3@t,eval3@TPR+eval3@TNR)
    eval[4,j]<-(subset(tab,tab[,2]==max(tab[,2])))[1,1]
    
    # récupération des données sur l'influence des param
    influ<-summary(model.res)
    #write.table(influ, "results/classic/influence_param.csv",append=T)
    
    # valeurs des prédictions sur lesquelles tombent les présences et les pseudoabsences 
    testvaluesp<-cbind(testvaluesp,testp)
    testvaluesa<-cbind(testvaluesa,testa)
    
    # evaluation des tests data
    values_test_pres <- extract(p,test_presences.occ)
    #hist(values_test_pres)
    
    # plot(p, col=palettecol)
    # points(map, type="l")
    # points(test_presences.occ, pch=20, col="pink")
    # points(test_back.occ,pch=20, col="grey")
    # 
    # evaluer par rapport à la valeur du maxSSS
    eval[8,j] <- 100*length(which(values_test_pres>eval[4,j]))/(length(values_test_pres)-length(which(is.na(values_test_pres))))
    eval[9, j] <- length(which(is.na(values_test_pres)))
    eval[10, j] <-length(values_test_pres)
    
    # partial dependance plot 
    y<- model.res$fitted
    x.partial <- model.res$gbm.call$dataframe
    dataframe.partial <- cbind (y,x.partial)
    names(dataframe.partial) <- c("y",names(x.partial))
    #write.table(dataframe.partial,paste("results/classic/partial/dataframe.partial",j,".csv"))
    
    #mesure de Moran 
    latlong_resi <- rbind(training_presences.occ,training_backg.occ)
    dist.mat_visit_resi <- as.matrix(geosphere::distm (latlong_resi))
    dist.mat_visit.inv_resi <- 1/dist.mat_visit_resi
    diag(dist.mat_visit.inv_resi ) <-0
    #dist.mat_visit.inv [1:5, 1:5]
    dist.mat_visit.inv_resi <- base::replace(dist.mat_visit.inv_resi,dist.mat_visit.inv_resi==Inf,0)
    
    # sur les résidus 
    eval_moran_resi <- ape::Moran.I(model.res$residuals,dist.mat_visit.inv_resi,na.rm=T) # numéro pixel
    eval[11,j] <- eval_moran_resi$observed
    eval[12,j] <- eval_moran_resi$p.value
    eval[13,j] <- eval_moran_resi$sd
    
    
  }
  
stack.pred <- dropLayer(stack.pred,1)

   return(list("stack"=stack.pred,"pdp"= dataframe.partial,"eval.file"= eval,"param.contri"= influ))
}