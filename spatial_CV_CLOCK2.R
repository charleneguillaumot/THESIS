#############################################################################################################
## Paper "Species distribution models in a data-poor and broad scale context", Progress in Oceanography
## 11/2018
## Guillaumot Charlène < charleneguillaumot21@gmail.com>
## SCRIPT for application of 2-fold CLOCK spatial cross-validation procedure
## use of Boosted Regression Trees
#############################################################################################################
library(raster)
library(dismo)
library(geosphere)
library(ape)

# Load the occurrence dataset
#--------------------------------
data <- read.csv("O_validus_occ_modif.csv", sep=";",dec=".", header=T)
odontaster.occ <- data[,c(1,2)]

#======================================================================================#
# Load the environmental descriptors
#-----------------------------------------
# RasterStack 
predictors_2005_2012 <- raster::stack("predictors_2005_2012_ANT_final_withcurrent.grd")
predictors_2005_2012 <- dropLayer(predictors_2005_2012,16) # drop a layer for this analysis
profGEBCO <- raster("profGEBCO_final.asc")

prof_max_species <- -1500 # define the maximal depth according to the deepest record, to restrain model extrapolation in deep areas where the species has neven been observed 

source("delim.area.R") # from SDMPlay package
extent_map <- extent(predictors_2005_2012)
predictors_2005_2012_crop_species <- delim.area (predictors = predictors_2005_2012,longmin = extent_map[1],longmax = extent_map[2], latmin =extent_map[3],latmax = extent_map[4],interval = c(0,prof_max_species))
predictors <- predictors_2005_2012_crop_species # crop the RasterStack to the extent of the wished projection (lat, long, depth)

#==============================================================================================#
# KDE layer of sampling effort, on which the background data will be sampled (by weighting) 
#----------------------------------------------------------------------------------------------
KDE <- raster("bias.grd")
#plot(KDE)
prof <- subset(predictors,1)
#plot(prof)
prof2 <- crop(prof, extent(KDE))
#plot(prof2)
extent(KDE) <- extent(prof2)
#plot(prof2)
KDE_mask <- mask(KDE, prof2)
#plot(KDE_mask)
KDE <- KDE_mask

# Préparation du modèle 
#----------------------------------
#boot <- 1 # nombre de réplicats de tirages de background
cv.boot <- 100 # nombre de CV

eval<-matrix(nrow=20,ncol=cv.boot,data=NA); colnames(eval)<-seq(1,ncol(eval),1);rownames(eval)<-c("AUC_cv","AUC_self","AUC_all","maxSSS","COR","pCOR","TSS","test_gp", "valid_test_pres","nb_NA_test","nb_test","moran_resi","pval_resi","sd_moran_resi","moran_extract","pval_extract","sd_moran_extract","borne_sup","borne_inf","ntrees") 
stack.pred<-subset(predictors,1);values(stack.pred)<-NA
testvaluesp<-rep(NA,nrow(unique(odontaster.occ))) ; testvaluesa<-testvaluesp


## récupérer les valeurs de l'environnement pour ce tirage 
#-----------------------------------------------------------
presvals_glob <- extract(predictors, unique(odontaster.occ)) # on extrait les données environnementales là où tombent les occurrences 
# on a la couche de profondeur qui est vraiment plus précise, donc on va extraire les données sur cette couche et remplacer la colonne de données pour la profondeur
extract_prof_presences_glob <- raster::extract(profGEBCO,unique(odontaster.occ))
presvals_glob_bis <- cbind(depth=extract_prof_presences_glob,presvals_glob[,-1])  # on raboute l'ancien tableau en changeant la colonne de profondeur
presvals.unique_glob <- replace(presvals_glob_bis,presvals_glob_bis[,1] <= prof_max_species, prof_max_species) # on remplace les profondeurs abbérantes par la profondeur maximale choisie pour l'espèce 

# on tire aléatoirement les données de background dans la zone pour définir le MESS 
background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
colnames(background_data) <- colnames(odontaster.occ)
# background
pseudoabs_glob1 <- extract(predictors,background_data)
pseudoabs_prof_glob <- extract(profGEBCO,background_data)
pseudoabs_glob <- cbind (depth=pseudoabs_prof_glob,pseudoabs_glob1[,-1])


for (j in 1:cv.boot){
  
  # on tire aléatoirement les données de background dans la zone pour définir le MESS 
  background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
  colnames(background_data) <- colnames(odontaster.occ)

  
  #library(ENMeval) # marche pas à cause d'un défaut dans rJava....
  # du coup on source directement la fonction qui nous intéresse 
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
  presence_tot <- unique(odontaster.occ)
  
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
  
  eval[18,j] <- borne_sup
  eval[19,j] <- borne_inf
  
  # plot(subset(predictors,1))
  # points(background_training, pch=20, col="grey")
  # points(presence_training, pch=20)
  # points(background_test, pch=20, col="pink")
  # points(presence_test, pch=20, col="darkred")
  
  training_presences.occ <- presence_training
  training_backg.occ <- background_training
  test_presences.occ <- presence_test
  test_back.occ <- background_test
  
  ####### LANCER LE MODELE SUR LES DONNEES DE TRAINING
  #-----------------------------------------------------------
  # Construire la matrice à implémenter dans BRT 
  #------------------------------------------------
  ### DONNEES DE PRESENCE
  
  presence_data <- training_presences.occ
  presvals <- extract(predictors, presence_data) # on extrait les données environnementales là où tombent les occurrences 
  #head(presvals) # 
  # on a la couche de profondeur qui est vraiment plus précise, donc on va extraire les données sur cette couche et remplacer la colonne de données pour la profondeur
  #profGEBCO   # à charger dans fichier #1
  extract_prof_presences <- raster::extract(profGEBCO,presence_data)
  #extract_prof_presences
  
  presvals_bis <- cbind(depth=extract_prof_presences,presvals[,-1])  # on raboute l'ancien tableau en changeant la colonne de profondeur
  presvals.unique <- replace(presvals_bis,presvals_bis[,1] <= prof_max_species, prof_max_species) # on remplace les profondeurs abbérantes par la profondeur maximale choisie pour l'espèce 
  
  
  ## DONNEES DE BACKGROUND
  #-----------------------
  pseudoabs1 <- extract(predictors,training_backg.occ)
  pseudoabs_prof <- extract(profGEBCO,training_backg.occ)
  pseudoabs <- cbind (depth=pseudoabs_prof,pseudoabs1[,-1])
  
  
  
  ## PREPARATION DE BRT 
  #-----------------------
  # Construction de la matrice combinant les renseignements sur les présences et les pseudo-absences (détail de la fonction SDMdata de SDMPlay)
  id<-0;sdmdata.unique<-0;  id<-c(rep(1,nrow(presvals.unique)),rep(0,nrow(pseudoabs))) 
  FICHIER<-data.frame(cbind(id,rbind(presvals.unique,pseudoabs)))
  
  # Boucle pour lancer BRT
  #-----------------------------------------------------------------------------
  # renseigner les paramètres de calibration de BRT
  tc=4   # tree complexity
  lr=0.007 # learning rate
  bf=0.75  # bag fraction
  
  # renseigner le dossier et le no du fichier où les résultats de contribution des variables envi seront stockés 
  #output <- "results2/influence_param.csv"
  
  #----------------------
  # Lancement de BRT
  #----------------------
  model.res<- gbm.step(data=FICHIER, 
                       gbm.x = 2:ncol(FICHIER),
                       gbm.y = 1,
                       family = "bernoulli",
                       tree.complexity = tc,
                       learning.rate = lr,
                       bag.fraction = bf)
  #----------------------
  
  #------------------------------------------------------------
  # Récupération des données de sortie du modèle 
  #------------------------------------------------------------
  # stack de prédiction 
  p<-predict(predictors,model.res,n.trees=model.res$gbm.call$best.trees,type="response", na.rm=F)
  p <- mask(p,subset(predictors,1))
  stack.pred<-stack(stack.pred,p)
  #plot(p, col=palettecol)
  #points(map, type="l")
  
  # get number of trees
  eval[20,j] <- model.res$n.trees
  
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
  #eval[11, j] <- eval3@auc
  
  # mesurer le TSS
  specificity <- (eval3@TPR/(eval3@TPR+eval3@FPR))
  sensitivity <- (eval3@TNR/(eval3@TNR+eval3@FNR))
  tss <- specificity + sensitivity -1
  eval[7,j]<- mean(tss, na.rm=T)
  
  # valeurs des maxSSS
  tab<-cbind(eval3@t,eval3@TPR+eval3@TNR)
  eval[4,j]<-(subset(tab,tab[,2]==max(tab[,2])))[1,1]
  
  # récupération des données sur l'influence des param
  influ<-summary(model.res)
  write.table(influ, "results/O_validus/CLOCK/influence_param.csv",append=T)
  
  # valeurs des prédictions sur lesquelles tombent les présences et les pseudoabsences 
  testvaluesp<-cbind(testvaluesp,testp)
  testvaluesa<-cbind(testvaluesa,testa)
  
  # evaluation des tests data
  #test_back.occ
  #test_presences.occ
  values_test_pres <- extract(p,test_presences.occ)
  #hist(values_test_pres)
  
  #plot(p, col=palettecol)
  #points(map, type="l")
  #points(test_presences.occ, pch=20, col="pink")
  #points(test_back.occ,pch=20, col="grey")
  
  # evaluer par rapport à la valeur du maxSSS
  eval[9,j] <- 100*length(which(values_test_pres>eval[4,j]))/(length(values_test_pres)-length(which(is.na(values_test_pres))))
  eval[10, j] <- length(which(is.na(values_test_pres)))
  eval[11, j] <-length(values_test_pres)
  
  # partial dependance plot 
  y<- model.res$fitted
  x.partial <- model.res$gbm.call$dataframe
  dataframe.partial <- cbind (y,x.partial)
  names(dataframe.partial) <- c("y",names(x.partial))
  write.table(dataframe.partial,paste("results/O_validus/CLOCK/partial/dataframe.partial",j,".csv"))
  
  #mesure de Moran 
  latlong_resi <- rbind(training_presences.occ,training_backg.occ)
  dist.mat_visit_resi <- as.matrix(geosphere::distm (latlong_resi))
  dist.mat_visit.inv_resi <- 1/dist.mat_visit_resi
  diag(dist.mat_visit.inv_resi ) <-0
  #dist.mat_visit.inv [1:5, 1:5]
  dist.mat_visit.inv_resi <- base::replace(dist.mat_visit.inv_resi,dist.mat_visit.inv_resi==Inf,0)
  
  # sur les résidus 
  eval_moran_resi <- ape::Moran.I(model.res$residuals,dist.mat_visit.inv_resi,na.rm=T) # numéro pixel
  eval[12,j] <- eval_moran_resi$observed
  eval[13,j] <- eval_moran_resi$p.value
  eval[14,j] <- eval_moran_resi$sd
  
  # sur les prédictions
  # # on tire 10000 points de background dans la zone
  # background_data_eval <- randomPoints(subset(predictors,1), n=10000)
  # 
  # latlong_extract <- background_data_eval
  # dist.mat_visit_extract <- as.matrix(geosphere::distm (latlong_extract))
  # dist.mat_visit.inv_extract <- 1/dist.mat_visit_extract
  # diag(dist.mat_visit.inv_extract ) <-0
  # #dist.mat_visit.inv [1:5, 1:5]
  # dist.mat_visit.inv_extract <- base::replace(dist.mat_visit.inv_extract,dist.mat_visit.inv_extract==Inf,0)
  # 
  # 
  # eval_moran_extract <- ape::Moran.I(extract(p,background_data_eval),dist.mat_visit.inv_extract,na.rm=T) # numéro pixel
  # eval[15,j] <- eval_moran_extract$observed
  # eval[16,j] <- eval_moran_extract$p.value
  # eval[17,j] <- eval_moran_extract$sd
  # 
}

write.csv(eval, "results/O_validus/CLOCK/eval.csv",append=T)

#mean_stack <- raster::calc(stack.pred, mean, na.rm=T)
#sd_stack <- raster::calc(stack.pred,sd, na.rm=T)

#raster::writeRaster(mean_stack, "results4/results2/stack.pred_mean1.asc")
#raster::writeRaster(sd_stack, "results4/results2/stack.pred_sd1.asc")
raster::writeRaster(stack.pred, "results/O_validus/CLOCK/stack_pred.grd")

