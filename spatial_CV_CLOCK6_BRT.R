#############################################################################################################
## Paper "Species distribution models in a data-poor and broad scale context", Progress in Oceanography
## 12/2018
## Guillaumot Charlène < charleneguillaumot21@gmail.com>
## SCRIPT for application of 6-fold CLOCK spatial cross-validation procedure
## use of Boosted Regression Trees
#############################################################################################################
library(raster)
library(dismo)
library(geosphere)
library(ape)
library(dplyr)

#-----------------------------
# Load the occurrence dataset
#-----------------------------
data <- read.csv("O_validus_occ_modif.csv", sep=";",dec=".", header=T)
odontaster.occ <- data[,c(1,2)]

#------------------------------------
# Load the environmental descriptors
#------------------------------------
# RasterStack 
predictors_2005_2012 <- raster::stack("predictors_2005_2012_ANT_final_withcurrent.grd")
predictors_2005_2012 <- dropLayer(predictors_2005_2012,16) # drop a layer for this analysis
profGEBCO <- raster("profGEBCO_final.asc")

prof_max_species <- -1500 # define the maximal depth according to the deepest record, to restrain model extrapolation in deep areas where the species has neven been observed 

source("delim.area.R") # from SDMPlay package
extent_map <- extent(predictors_2005_2012)
predictors_2005_2012_crop_species <- delim.area (predictors = predictors_2005_2012,longmin = extent_map[1],longmax = extent_map[2], latmin =extent_map[3],latmax = extent_map[4],interval = c(0,prof_max_species))
predictors <- predictors_2005_2012_crop_species # crop the RasterStack to the extent of the desired projection (lat, long, depth)

#--------------------------------------------------------------------------------------------------------------
# KDE layer of sampling effort, on which the background data will be sampled (by weighting) 
# The KDE (Kernel Density Estimation) is a statistical tool that helps to measure the probability of finding
# an occurrencce on each pixel, according to the set of benthic records sampled in the entire Southern Ocean 
# (update from the Biogeographic Atlas of the Southern Ocean, supp.mat. of this study)
#--------------------------------------------------------------------------------------------------------------
KDE <- raster("bias.grd")
prof2 <- crop(prof, extent(KDE))
extent(KDE) <- extent(prof2)
KDE_mask <- mask(KDE, prof2)
KDE <- KDE_mask

#-------------------
# Model preparation 
#-------------------
cv.boot <- 100 # number of model replicates 

# create empty matrices, vectors and rasters that will be filled in the following loop 
eval<-matrix(nrow=20,ncol=cv.boot,data=NA); colnames(eval)<-seq(1,ncol(eval),1);rownames(eval)<-c("AUC_cv","AUC_self","AUC_all","maxSSS","COR","pCOR","TSS","test_gp", "valid_test_pres","nb_NA_test","nb_test","moran_resi","pval_resi","sd_moran_resi","moran_extract","pval_extract","sd_moran_extract","borne_sup","borne_inf","ntrees") 
stack.pred<-subset(predictors,1);values(stack.pred)<-NA
testvaluesp<-rep(NA,nrow(unique(odontaster.occ))) ; testvaluesa<-testvaluesp

## Initialise the matrix that contains occurrence records (presence and background data + extract of the environmental values contained in the pixels where occurrences are found)
#--------------------------------------------------------------------------------------------------------------------
presvals_glob <- extract(predictors, unique(odontaster.occ)) # extract environmental data where presences are found. The unique function restricts the extraction of a single occurrence in each pixel. Duplicates are removed. 
extract_prof_presences_glob <- raster::extract(profGEBCO,unique(odontaster.occ))
presvals_glob_bis <- cbind(depth=extract_prof_presences_glob,presvals_glob[,-1])
# the GEBCO layer is a bathymetry layer that is far more precise that the other one contained in the RasterStack. So, we extract the bathymetry value to calibrate the model from this layer but predict on the other bathymetry layer (due to resolution)
presvals_glob_bis <- cbind(depth=extract_prof_presences_glob,presvals_glob[,-1])  
presvals.unique_glob <- replace(presvals_glob_bis,presvals_glob_bis[,1] <= prof_max_species, prof_max_species) # removal of weird georefenced data (replacement of far too deep records)

# sampling of background data 
#-----------------------------
# 1000 background data are randomly sampled in the environment, according to the weighting scheme of the KDE layer 
background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
colnames(background_data) <- colnames(odontaster.occ)
# extract environmental conditions where the background data are sampled 
pseudoabs_glob1 <- extract(predictors,background_data)
pseudoabs_prof_glob <- extract(profGEBCO,background_data)
pseudoabs_glob <- cbind (depth=pseudoabs_prof_glob,pseudoabs_glob1[,-1])

for (j in 1:cv.boot){
  # sampling of new background data for each loop step 
  background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
  colnames(background_data) <- colnames(odontaster.occ)
  
  # the spatial sampling aims at defining areas containing either test or training data 
  # for the CLOCK-6 method, 6 triangle areas are defined, cutting the Southern Ocean circle into 6 equal areas. Three of these 6 areas contain the data that will be used to train the model, the other remaining three ones, the data that will be used to test the model after predictions 
  # the areas used to train and test the model are randomly defined at each loop step 
  
  #initialise vectors that are used to generate the spatial sampling structure
  random_longA <- seq(0,179,1)
  random_longB <- seq(-180,-1,1)
  random_longC <- seq(0,180,1)
  random_longD <- seq(-179,0,1)
  random_long <- c(random_longA,random_longB,random_longC,random_longD)
  
  # sample a number between -180 and 180 to define the random sampling transect 
  tirage <- sample(seq(181,541,1),1)
  random_long_tirage <- random_long[tirage]
  random_long_tirage_A1 <- random_long[tirage+60]
  random_long_tirage_A2 <- random_long[tirage+120]
  random_long_tirage_A3 <- random_long[tirage+180]
  random_long_tirage_A4 <- random_long[tirage-120]
  random_long_tirage_A5 <- random_long[tirage-60]
  random_long_tirage_A6 <- random_long_tirage
  
  
  ## define training and test groups (composed of both presence and background data)
  presence_tot <- unique(odontaster.occ)
  training_presences.occ1 <- NA;training_presences.occ3 <-NA;training_presences.occ5 <-NA
  training_backgr.occ1 <- NA; training_backgr.occ3 <- NA;training_backgr.occ5 <-NA
  test_presences.occ2 <- NA;test_presences.occ4 <- NA;test_presences.occ6 <- NA
  test_backgr.occ2 <- NA;test_backgr.occ4 <- NA;test_backgr.occ6 <- NA
  

  #---------------------
  # TRAINING PRESENCES
  #---------------------
  
  ### ZONE 1 ####
  
  if (random_long_tirage < 0 & random_long_tirage_A1 <0){
    training_presences.occ1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < random_long_tirage_A1 )
    training_backgr.occ1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < random_long_tirage_A1)
  } 
  if (random_long_tirage > 0 & random_long_tirage_A1 >0){
    training_presences.occ1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < random_long_tirage_A1 )
    training_backgr.occ1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < random_long_tirage_A1)
  }
  
  if (random_long_tirage > 0 & random_long_tirage_A1 <0){
    training_presences.occ1_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < 180 )
    training_presences.occ1_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A1 )
    training_presences.occ1 <- rbind(training_presences.occ1_part1,training_presences.occ1_part2)
    
    training_backgr.occ1_part1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < 180 )
    training_backgr.occ1_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A1 )
    training_backgr.occ1 <- rbind(training_backgr.occ1_part1,training_backgr.occ1_part2)
  }
  
  if (random_long_tirage < 0 & random_long_tirage_A1 >0){
    training_presences.occ1_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage & presence_tot[,1] < 0 )
    training_presences.occ1_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A1 )
    training_presences.occ1 <- rbind(training_presences.occ1_part1,training_presences.occ1_part2)
    
    training_backgr.occ1_part1 <- subset(background_data,background_data[,1] > random_long_tirage & background_data[,1] < 0 )
    training_backgr.occ1_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A1 )
    training_backgr.occ1 <- rbind(training_backgr.occ1_part1,training_backgr.occ1_part2)
  }
  
  ### ZONE 2 ####
  
  if (random_long_tirage_A1 < 0 & random_long_tirage_A2 <0){
    training_presences.occ2 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < random_long_tirage_A2 )
    training_backgr.occ2 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < random_long_tirage_A2)
  } 
  if (random_long_tirage_A1 > 0 & random_long_tirage_A2 >0){
    training_presences.occ2 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < random_long_tirage_A2 )
    training_backgr.occ2 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < random_long_tirage_A2)
  }
  
  if (random_long_tirage_A1 > 0 & random_long_tirage_A2 <0){
    training_presences.occ2_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < 180 )
    training_presences.occ2_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A2 )
    training_presences.occ2 <- rbind(training_presences.occ2_part1,training_presences.occ2_part2)
    
    training_backgr.occ2_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < 180 )
    training_backgr.occ2_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A2 )
    training_backgr.occ2 <- rbind(training_backgr.occ2_part1,training_backgr.occ2_part2)
  }
  
  if (random_long_tirage_A1 < 0 & random_long_tirage_A2 >0){
    training_presences.occ2_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A1 & presence_tot[,1] < 0 )
    training_presences.occ2_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A2 )
    training_presences.occ2 <- rbind(training_presences.occ2_part1,training_presences.occ2_part2)
    
    training_backgr.occ2_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A1 & background_data[,1] < 0 )
    training_backgr.occ2_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A2 )
    training_backgr.occ2 <- rbind(training_backgr.occ2_part1,training_backgr.occ2_part2)
  }
  
  ### ZONE 3 ####
  
  if (random_long_tirage_A2 < 0 & random_long_tirage_A3 <0){
    training_presences.occ3 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < random_long_tirage_A3 )
    training_backgr.occ3 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < random_long_tirage_A3)
  } 
  if (random_long_tirage_A2 > 0 & random_long_tirage_A3 >0){
    training_presences.occ3 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < random_long_tirage_A3 )
    training_backgr.occ3 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < random_long_tirage_A3)
  }
  
  if (random_long_tirage_A2 > 0 & random_long_tirage_A3 <0){
    training_presences.occ3_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < 180 )
    training_presences.occ3_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A3 )
    training_presences.occ3 <- rbind(training_presences.occ3_part1,training_presences.occ3_part2)
    
    training_backgr.occ3_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < 180 )
    training_backgr.occ3_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A3 )
    training_backgr.occ3 <- rbind(training_backgr.occ3_part1,training_backgr.occ3_part2)
  }
  
  if (random_long_tirage_A2 < 0 & random_long_tirage_A3 >0){
    training_presences.occ3_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A2 & presence_tot[,1] < 0 )
    training_presences.occ3_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A3 )
    training_presences.occ3 <- rbind(training_presences.occ3_part1,training_presences.occ3_part2)
    
    training_backgr.occ3_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A2 & background_data[,1] < 0 )
    training_backgr.occ3_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A3 )
    training_backgr.occ3 <- rbind(training_backgr.occ3_part1,training_backgr.occ3_part2)
  }
  
  ### ZONE 4 ####
  
  if (random_long_tirage_A3 < 0 & random_long_tirage_A4 <0){
    training_presences.occ4 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < random_long_tirage_A4 )
    training_backgr.occ4 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < random_long_tirage_A4)
  } 
  if (random_long_tirage_A3 > 0 & random_long_tirage_A4 >0){
    training_presences.occ4 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < random_long_tirage_A4 )
    training_backgr.occ4 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < random_long_tirage_A4)
  }
  
  if (random_long_tirage_A3 > 0 & random_long_tirage_A4 <0){
    training_presences.occ4_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < 180 )
    training_presences.occ4_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A4 )
    training_presences.occ4 <- rbind(training_presences.occ4_part1,training_presences.occ4_part2)
    
    training_backgr.occ4_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < 180 )
    training_backgr.occ4_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A4 )
    training_backgr.occ4 <- rbind(training_backgr.occ4_part1,training_backgr.occ4_part2)
  }
  
  if (random_long_tirage_A3 < 0 & random_long_tirage_A4 >0){
    training_presences.occ4_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A3 & presence_tot[,1] < 0 )
    training_presences.occ4_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A4 )
    training_presences.occ4 <- rbind(training_presences.occ4_part1,training_presences.occ4_part2)
    
    training_backgr.occ4_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A3 & background_data[,1] < 0 )
    training_backgr.occ4_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A4 )
    training_backgr.occ4 <- rbind(training_backgr.occ4_part1,training_backgr.occ4_part2)
  }
  
  ### ZONE 5 ####
  
  if (random_long_tirage_A4 < 0 & random_long_tirage_A5 <0){
    training_presences.occ5 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < random_long_tirage_A5 )
    training_backgr.occ5 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < random_long_tirage_A5)
  } 
  if (random_long_tirage_A4 > 0 & random_long_tirage_A5 >0){
    training_presences.occ5 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < random_long_tirage_A5 )
    training_backgr.occ5 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < random_long_tirage_A5)
  }
  
  if (random_long_tirage_A4 > 0 & random_long_tirage_A5 <0){
    training_presences.occ5_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < 180 )
    training_presences.occ5_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A5 )
    training_presences.occ5 <- rbind(training_presences.occ5_part1,training_presences.occ5_part2)
    
    training_backgr.occ5_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < 180 )
    training_backgr.occ5_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A5 )
    training_backgr.occ5 <- rbind(training_backgr.occ5_part1,training_backgr.occ5_part2)
  }
  
  if (random_long_tirage_A4 < 0 & random_long_tirage_A5 >0){
    training_presences.occ5_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A4 & presence_tot[,1] < 0 )
    training_presences.occ5_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A5 )
    training_presences.occ5 <- rbind(training_presences.occ5_part1,training_presences.occ5_part2)
    
    training_backgr.occ5_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A4 & background_data[,1] < 0 )
    training_backgr.occ5_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A5 )
    training_backgr.occ5 <- rbind(training_backgr.occ5_part1,training_backgr.occ5_part2)
  }
  
  ### ZONE 6 ####
  
  if (random_long_tirage_A5 < 0 & random_long_tirage_A6 <0){
    training_presences.occ6 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < random_long_tirage_A6 )
    training_backgr.occ6 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < random_long_tirage_A6)
  } 
  if (random_long_tirage_A5 > 0 & random_long_tirage_A6 >0){
    training_presences.occ6 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < random_long_tirage_A6 )
    training_backgr.occ6 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < random_long_tirage_A6)
  }
  
  if (random_long_tirage_A5 > 0 & random_long_tirage_A6 <0){
    training_presences.occ6_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < 180 )
    training_presences.occ6_part2 <- subset(presence_tot,presence_tot[,1] > -180 & presence_tot[,1] < random_long_tirage_A6 )
    training_presences.occ6 <- rbind(training_presences.occ6_part1,training_presences.occ6_part2)
    
    training_backgr.occ6_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < 180 )
    training_backgr.occ6_part2 <- subset(background_data,background_data[,1] > -180 & background_data[,1] < random_long_tirage_A6 )
    training_backgr.occ6 <- rbind(training_backgr.occ6_part1,training_backgr.occ6_part2)
  }
  
  if (random_long_tirage_A5 < 0 & random_long_tirage_A6 >0){
    training_presences.occ6_part1 <- subset(presence_tot,presence_tot[,1] > random_long_tirage_A5 & presence_tot[,1] < 0 )
    training_presences.occ6_part2 <- subset(presence_tot,presence_tot[,1] > 0 & presence_tot[,1] < random_long_tirage_A6 )
    training_presences.occ6 <- rbind(training_presences.occ6_part1,training_presences.occ6_part2)
    
    training_backgr.occ6_part1 <- subset(background_data,background_data[,1] > random_long_tirage_A5 & background_data[,1] < 0 )
    training_backgr.occ6_part2 <- subset(background_data,background_data[,1] > 0 & background_data[,1] < random_long_tirage_A6 )
    training_backgr.occ6 <- rbind(training_backgr.occ6_part1,training_backgr.occ6_part2)
  }
  
  #---------------------
  # TEST PRESENCE
  #---------------------
  test_presences.occ <- rbind(training_presences.occ1,training_presences.occ3,training_presences.occ5)
  test_back.occ <- rbind(training_backgr.occ1,training_backgr.occ3,training_backgr.occ5)
  
  training_presences.occ <- rbind(training_presences.occ2,training_presences.occ4,training_presences.occ6)
  training_backgr.occ <- rbind(training_backgr.occ2,training_backgr.occ4,training_backgr.occ6)
  
  # record the positions of the transect defined for sampling 
  borne <- random_long_tirage
  eval[18,j] <- borne

  ####### LAUNCH THE SDM ON TRAINING DATA #######
  #---------------------------------------------------------------------------------------------------
  # Build the matrix containing lat long data and the environmental values associated
  #---------------------------------------------------------------------------------------------------
  # presence data
  presence_data <- training_presences.occ
  presvals <- extract(predictors, presence_data) 
  extract_prof_presences <- raster::extract(profGEBCO,presence_data)
  presvals_bis <- cbind(depth=extract_prof_presences,presvals[,-1])  
  presvals.unique <- replace(presvals_bis,presvals_bis[,1] <= prof_max_species, prof_max_species) 
  
  ## background data
  #-----------------------
  pseudoabs1 <- extract(predictors,training_backg.occ)
  pseudoabs_prof <- extract(profGEBCO,training_backg.occ)
  pseudoabs <- cbind (depth=pseudoabs_prof,pseudoabs1[,-1])
  
  #-----------------------
  ## BRT calibration 
  #-----------------------
  # Initialise the matrix containing presence, background data and the environmental values associated 
  id<-0;sdmdata.unique<-0;  id<-c(rep(1,nrow(presvals.unique)),rep(0,nrow(pseudoabs))) 
  FICHIER<-data.frame(cbind(id,rbind(presvals.unique,pseudoabs)))
  
  # BRT parameters
  tc=4   # tree complexity
  lr=0.007 # learning rate
  bf=0.75  # bag fraction
  
  #-----------------------
  # BRT launch
  #----------------------
  model.res<- gbm.step(data=FICHIER, 
                       gbm.x = 2:ncol(FICHIER),
                       gbm.y = 1,
                       family = "bernoulli",
                       tree.complexity = tc,
                       learning.rate = lr,
                       bag.fraction = bf)
  
  #------------------------------------------------------------
  # Get model outputs 
  #------------------------------------------------------------
  # Predicted map  
  p<-predict(predictors,model.res,n.trees=model.res$gbm.call$best.trees,type="response", na.rm=F)
  p <- mask(p,subset(predictors,1))
  stack.pred<-stack(stack.pred,p) # stack all the maps replicates 
  
  # get number of trees
  eval[20,j] <- model.res$n.trees
  
  # Area Under the Curve, Point Biserial Correlation
  testp<-dismo::predict(model.res,data.frame(presvals.unique_glob),n.trees=model.res$gbm.call$best.trees,type="response")
  testa<-dismo::predict(model.res,data.frame(pseudoabs_glob),n.trees=model.res$gbm.call$best.trees,type="response")
  eval3<-dismo::evaluate(p=testp,a=testa)
  
  eval[3,j]<-eval3@auc # AUC
  eval[5,j]<-eval3@cor # point biserial correlation (COR)
  eval[6,j]<-eval3@pcor # probability of the COR
  eval[1,j]<-model.res$cv.statistics$discrimination.mean # AUC
  eval[2,j]<-model.res$self.statistics$discrimination # internal AUC 
  
  # True Skill Statistics (TSS)
  specificity <- (eval3@TPR/(eval3@TPR+eval3@FPR))
  sensitivity <- (eval3@TNR/(eval3@TNR+eval3@FNR))
  tss <- specificity + sensitivity -1
  eval[7,j]<- mean(tss, na.rm=T)
  
  # maxSSS: maximum sensitivity plus specificity threshold
  tab<-cbind(eval3@t,eval3@TPR+eval3@TNR)
  eval[4,j]<-(subset(tab,tab[,2]==max(tab[,2])))[1,1]
  
  # Contribution of environmental descriptors
  influ<-summary(model.res)
  write.table(influ, "results/O_validus/CLOCK/influence_param.csv",append=T)
  
  # Percentage of test data correcty predicted 
  testvaluesp<-cbind(testvaluesp,testp)
  testvaluesa<-cbind(testvaluesa,testa)
  
  values_test_pres <- extract(p,test_presences.occ)
  eval[9,j] <- 100*length(which(values_test_pres>eval[4,j]))/(length(values_test_pres)-length(which(is.na(values_test_pres))))
  eval[10, j] <- length(which(is.na(values_test_pres)))
  eval[11, j] <-length(values_test_pres)
  
  # partial dependance plots 
  y<- model.res$fitted
  x.partial <- model.res$gbm.call$dataframe
  dataframe.partial <- cbind (y,x.partial)
  names(dataframe.partial) <- c("y",names(x.partial))
  write.table(dataframe.partial,paste("results/O_validus/CLOCK/partial/dataframe.partial",j,".csv"))
  
  # spatial autocorrelation: Moran's I calculated on models residuals
  latlong_resi <- rbind(training_presences.occ,training_backg.occ)
  dist.mat_visit_resi <- as.matrix(geosphere::distm (latlong_resi))
  dist.mat_visit.inv_resi <- 1/dist.mat_visit_resi
  diag(dist.mat_visit.inv_resi ) <-0
  dist.mat_visit.inv_resi <- base::replace(dist.mat_visit.inv_resi,dist.mat_visit.inv_resi==Inf,0)
  
  eval_moran_resi <- ape::Moran.I(model.res$residuals,dist.mat_visit.inv_resi,na.rm=T) # numéro pixel
  eval[12,j] <- eval_moran_resi$observed
  eval[13,j] <- eval_moran_resi$p.value
  eval[14,j] <- eval_moran_resi$sd
  
}

#------------------------------------------------------------
# Save results 
#------------------------------------------------------------
write.csv(eval, "results/O_validus/CLOCK/eval.csv",append=T)

mean_stack <- raster::calc(stack.pred, mean, na.rm=T)
sd_stack <- raster::calc(stack.pred,sd, na.rm=T)
raster::writeRaster(mean_stack, "results4/results2/stack.pred_mean1.asc")
raster::writeRaster(sd_stack, "results4/results2/stack.pred_sd1.asc")
raster::writeRaster(stack.pred, "results/O_validus/CLOCK/stack_pred.grd")

