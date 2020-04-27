clock6 <- function(predictors, KDE_layer, data, nb_replicates) {
  
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
    
    # initialise the random sampling
    random_longA <- seq(0,179,1)
    random_longB <- seq(-180,-1,1)
    random_longC <- seq(0,180,1)
    random_longD <- seq(-179,0,1)
    random_long <- c(random_longA,random_longB,random_longC,random_longD)
    
    # sample a random number between -180 and 180 to define the inferior boundary of the sampling 
    tirage <- sample(seq(181,541,1),1)
    random_long_tirage <- random_long[tirage]
    random_long_tirage_A1 <- random_long[tirage+60]
    random_long_tirage_A2 <- random_long[tirage+120]
    random_long_tirage_A3 <- random_long[tirage+180]
    random_long_tirage_A4 <- random_long[tirage-120]
    random_long_tirage_A5 <- random_long[tirage-60]
    random_long_tirage_A6 <- random_long_tirage
    
    eval[15,j]<-tirage
    
    # define the test and training groups 
    presence_tot <- unique(data)
    
    training_presences.occ1 <- NA;training_presences.occ3 <-NA;training_presences.occ5 <-NA
    training_backgr.occ1 <- NA; training_backgr.occ3 <- NA;training_backgr.occ5 <-NA
    test_presences.occ2 <- NA;test_presences.occ4 <- NA;test_presences.occ6 <- NA
    test_backgr.occ2 <- NA;test_backgr.occ4 <- NA;test_backgr.occ6 <- NA
    
    
    #---------------------
    # TRAINING PRESENCES
    #---------------------
    
    ### ZONE 1 ####
    
    if (random_long_tirage <= 0 & random_long_tirage_A1 <=0){
      training_presences.occ1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < random_long_tirage_A1 )
      training_backgr.occ1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < random_long_tirage_A1)
    } 
    if (random_long_tirage >= 0 & random_long_tirage_A1 >=0){
      training_presences.occ1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < random_long_tirage_A1 )
      training_backgr.occ1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < random_long_tirage_A1)
    }
    
    if (random_long_tirage >= 0 & random_long_tirage_A1 <=0){
      training_presences.occ1_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < 180 )
      training_presences.occ1_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A1 )
      training_presences.occ1 <- rbind(training_presences.occ1_part1,training_presences.occ1_part2)
      
      training_backgr.occ1_part1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < 180 )
      training_backgr.occ1_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A1 )
      training_backgr.occ1 <- rbind(training_backgr.occ1_part1,training_backgr.occ1_part2)
    }
    
    if (random_long_tirage <= 0 & random_long_tirage_A1 >=0){
      training_presences.occ1_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < 0 )
      training_presences.occ1_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A1 )
      training_presences.occ1 <- rbind(training_presences.occ1_part1,training_presences.occ1_part2)
      
      training_backgr.occ1_part1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < 0 )
      training_backgr.occ1_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A1 )
      training_backgr.occ1 <- rbind(training_backgr.occ1_part1,training_backgr.occ1_part2)
    }
    
    ### ZONE 2 ####
    
    if (random_long_tirage_A1 <= 0 & random_long_tirage_A2 <=0){
      training_presences.occ2 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < random_long_tirage_A2 )
      training_backgr.occ2 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < random_long_tirage_A2)
    } 
    if (random_long_tirage_A1 >= 0 & random_long_tirage_A2 >=0){
      training_presences.occ2 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < random_long_tirage_A2 )
      training_backgr.occ2 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < random_long_tirage_A2)
    }
    
    if (random_long_tirage_A1 >= 0 & random_long_tirage_A2 <=0){
      training_presences.occ2_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < 180 )
      training_presences.occ2_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A2 )
      training_presences.occ2 <- rbind(training_presences.occ2_part1,training_presences.occ2_part2)
      
      training_backgr.occ2_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < 180 )
      training_backgr.occ2_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A2 )
      training_backgr.occ2 <- rbind(training_backgr.occ2_part1,training_backgr.occ2_part2)
    }
    
    if (random_long_tirage_A1 <= 0 & random_long_tirage_A2 >=0){
      training_presences.occ2_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < 0 )
      training_presences.occ2_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A2 )
      training_presences.occ2 <- rbind(training_presences.occ2_part1,training_presences.occ2_part2)
      
      training_backgr.occ2_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < 0 )
      training_backgr.occ2_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A2 )
      training_backgr.occ2 <- rbind(training_backgr.occ2_part1,training_backgr.occ2_part2)
    }
    
    ### ZONE 3 ####
    
    if (random_long_tirage_A2 <= 0 & random_long_tirage_A3 <=0){
      training_presences.occ3 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < random_long_tirage_A3 )
      training_backgr.occ3 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < random_long_tirage_A3)
    } 
    if (random_long_tirage_A2 >= 0 & random_long_tirage_A3 >=0){
      training_presences.occ3 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < random_long_tirage_A3 )
      training_backgr.occ3 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < random_long_tirage_A3)
    }
    
    if (random_long_tirage_A2 >= 0 & random_long_tirage_A3 <=0){
      training_presences.occ3_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < 180 )
      training_presences.occ3_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A3 )
      training_presences.occ3 <- rbind(training_presences.occ3_part1,training_presences.occ3_part2)
      
      training_backgr.occ3_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < 180 )
      training_backgr.occ3_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A3 )
      training_backgr.occ3 <- rbind(training_backgr.occ3_part1,training_backgr.occ3_part2)
    }
    
    if (random_long_tirage_A2 <= 0 & random_long_tirage_A3 >=0){
      training_presences.occ3_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < 0 )
      training_presences.occ3_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A3 )
      training_presences.occ3 <- rbind(training_presences.occ3_part1,training_presences.occ3_part2)
      
      training_backgr.occ3_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < 0 )
      training_backgr.occ3_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A3 )
      training_backgr.occ3 <- rbind(training_backgr.occ3_part1,training_backgr.occ3_part2)
    }
    
    ### ZONE 4 ####
    
    if (random_long_tirage_A3 <= 0 & random_long_tirage_A4 <=0){
      training_presences.occ4 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < random_long_tirage_A4 )
      training_backgr.occ4 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < random_long_tirage_A4)
    } 
    if (random_long_tirage_A3 >= 0 & random_long_tirage_A4 >=0){
      training_presences.occ4 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < random_long_tirage_A4 )
      training_backgr.occ4 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < random_long_tirage_A4)
    }
    
    if (random_long_tirage_A3 >= 0 & random_long_tirage_A4 <=0){
      training_presences.occ4_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < 180 )
      training_presences.occ4_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A4 )
      training_presences.occ4 <- rbind(training_presences.occ4_part1,training_presences.occ4_part2)
      
      training_backgr.occ4_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < 180 )
      training_backgr.occ4_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A4 )
      training_backgr.occ4 <- rbind(training_backgr.occ4_part1,training_backgr.occ4_part2)
    }
    
    if (random_long_tirage_A3 <= 0 & random_long_tirage_A4 >=0){
      training_presences.occ4_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < 0 )
      training_presences.occ4_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A4 )
      training_presences.occ4 <- rbind(training_presences.occ4_part1,training_presences.occ4_part2)
      
      training_backgr.occ4_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < 0 )
      training_backgr.occ4_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A4 )
      training_backgr.occ4 <- rbind(training_backgr.occ4_part1,training_backgr.occ4_part2)
    }
    
    ### ZONE 5 ####
    
    if (random_long_tirage_A4 <= 0 & random_long_tirage_A5 <=0){
      training_presences.occ5 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < random_long_tirage_A5 )
      training_backgr.occ5 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < random_long_tirage_A5)
    } 
    if (random_long_tirage_A4 >= 0 & random_long_tirage_A5 >=0){
      training_presences.occ5 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < random_long_tirage_A5 )
      training_backgr.occ5 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < random_long_tirage_A5)
    }
    
    if (random_long_tirage_A4 >= 0 & random_long_tirage_A5 <=0){
      training_presences.occ5_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < 180 )
      training_presences.occ5_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A5 )
      training_presences.occ5 <- rbind(training_presences.occ5_part1,training_presences.occ5_part2)
      
      training_backgr.occ5_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < 180 )
      training_backgr.occ5_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A5 )
      training_backgr.occ5 <- rbind(training_backgr.occ5_part1,training_backgr.occ5_part2)
    }
    
    if (random_long_tirage_A4 <= 0 & random_long_tirage_A5 >=0){
      training_presences.occ5_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < 0 )
      training_presences.occ5_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A5 )
      training_presences.occ5 <- rbind(training_presences.occ5_part1,training_presences.occ5_part2)
      
      training_backgr.occ5_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < 0 )
      training_backgr.occ5_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A5 )
      training_backgr.occ5 <- rbind(training_backgr.occ5_part1,training_backgr.occ5_part2)
    }
    
    ### ZONE 6 ####
    
    if (random_long_tirage_A5 <= 0 & random_long_tirage_A6 <=0){
      training_presences.occ6 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < random_long_tirage_A6 )
      training_backgr.occ6 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < random_long_tirage_A6)
    } 
    if (random_long_tirage_A5 >= 0 & random_long_tirage_A6 >=0){
      training_presences.occ6 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < random_long_tirage_A6 )
      training_backgr.occ6 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < random_long_tirage_A6)
    }
    
    if (random_long_tirage_A5 >= 0 & random_long_tirage_A6 <=0){
      training_presences.occ6_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < 180 )
      training_presences.occ6_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A6 )
      training_presences.occ6 <- rbind(training_presences.occ6_part1,training_presences.occ6_part2)
      
      training_backgr.occ6_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < 180 )
      training_backgr.occ6_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A6 )
      training_backgr.occ6 <- rbind(training_backgr.occ6_part1,training_backgr.occ6_part2)
    }
    
    if (random_long_tirage_A5 <= 0 & random_long_tirage_A6 >=0){
      training_presences.occ6_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < 0 )
      training_presences.occ6_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A6 )
      training_presences.occ6 <- rbind(training_presences.occ6_part1,training_presences.occ6_part2)
      
      training_backgr.occ6_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < 0 )
      training_backgr.occ6_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A6 )
      training_backgr.occ6 <- rbind(training_backgr.occ6_part1,training_backgr.occ6_part2)
    }
    
    # palette.bleue<-colorRampPalette(c("blue4","blue","dodgerblue", "deepskyblue","lightskyblue"))(800)
    # plot(subset(predictors,1), col=palette.bleue)
    # points(bathybiaster.occ, col="black",pch=20)
    # points(training_backgr.occ1, col="pink", pch=16)
    # points(training_presences.occ1, col="red", pch=20)
    # points(training_backgr.occ2, col="palegreen1", pch=16)
    # points(training_presences.occ2, col="green4", pch=20)
    # points(training_backgr.occ3, col="grey", pch=16)
    # points(training_presences.occ3, col="black", pch=20)
    # points(training_backgr.occ4, col="yellow", pch=16)
    # points(training_presences.occ4, col="orange", pch=20)
    # points(training_backgr.occ5, col="cadetblue1", pch=16)
    # points(training_presences.occ5, col="blue", pch=20)
    # points(training_backgr.occ6, col="mediumorchid1", pch=16)
    # points(training_presences.occ6, col="blueviolet", pch=20)
    # points(bathybiaster.occ, col="black",pch=20)
    # 
    #---------------------
    # TEST PRESENCE
    #---------------------
    test_presences.occ <- rbind(training_presences.occ1,training_presences.occ3,training_presences.occ5)
    test_back.occ <- rbind(training_backgr.occ1,training_backgr.occ3,training_backgr.occ5)
    
    training_presences.occ <- rbind(training_presences.occ2,training_presences.occ4,training_presences.occ6)
    training_backgr.occ <- rbind(training_backgr.occ2,training_backgr.occ4,training_backgr.occ6)
    
    # plot(subset(predictors,1), col=palette.bleue)
    # points(test_back.occ, col="palegreen1", pch=20)
    # points(test_presences.occ, col="green4", pch=20)
    # points(training_backgr.occ, col="pink", pch=20)
    # points(training_presences.occ, col="red", pch=20)
    
    
    ####### LAUNCH THE MODEL ON THE TRAINING DATA
    #-----------------------------------------------------------
    # BUILD THE MATRIX TO CALIBRATE BRST
    #------------------------------------------------
    ### PRESENCE DATA
    
    presence_data <- training_presences.occ
    presvals.unique <- extract(predictors, presence_data) 
    
    
    ## BACKGROUND DATA
    #-----------------------
    pseudoabs <- extract(predictors,training_backgr.occ)
    
    
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
    latlong_resi <- rbind(training_presences.occ,training_backgr.occ)
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