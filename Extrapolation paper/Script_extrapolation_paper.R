#--------------------------------------------------------------------------------------------------------------------------

# Script GUILLAUMOT C. et al. (2020) for paper "Extrapolation in species distribution modelling. Application to Southern Ocean marine species", Progress in Oceanography

# Example for Acodontaster hodgsoni, analysis #1 (model ran with delimiting the projection area to limit in depth). This code also includes comments to describe the parts to be modified to apply Analysis #0 (no limit in depth), or the increase/decrease of presence-only records (historical analysis).

#--------------------------------------------------------------------------------------------------------------------------
library(dismo)
library(raster)

#-----------------------------------
# open occurrences data 
#-----------------------------------
data <-read.csv("DATA/database2.csv", header=T, sep=";")
head(data)
species <- subset(data, data$ScientificName_accepted=="Acodontaster hodgsoni")
#species <- subset(species,species$Year <= 1980) # when running the model following the historical sampling (change in the number of presence-only records available according to the time of their sampling)
species.occ <- species[,c(6,5)]

#-----------------------------------
# open environmental database
#-----------------------------------
pred <- stack("DATA/stack_intermediate.grd")
pred <- dropLayer(pred,c(9:13,15:21,25,27,28,30:58))
depth <- subset(pred,1)

min_depth <- -1500  # indicate the maximal bathymetry where presence records can be found
source("DATA/delim.area.R")
pred <- delim.area(pred,longmin =extent(depth)[1], longmax=extent(depth)[2],latmin=extent(depth)[3], latmax=extent(depth)[4], interval=c(min_depth,0) )
# when running Analysis #0 (without cropping to depth), these 4 above lines are deleted 

pred <- mask(pred,subset(pred,1)) # delete the terrestrial areas (replace by the NA values of the "depth" layer)

#-----------------------------------------
# Load the KDE layer
# The KDE was created with the kde2d function of the MASS R package. It uses the ensemble of visited areas in Antarctica (see De Broyer et al. 2014, Guillaumot et al. 2019 supplementay material)
KDE <- raster("DATA/bias_OK.asc")
prof <- subset(pred,1)
prof2 <- crop(prof, extent(KDE))
extent(KDE) <- extent(prof2)
KDE_mask <- mask(KDE, prof2)
KDE <- KDE_mask

## Initialize empty tables that will be filled
prop <- 100 # Will be modified when less that 100% of the presence data are used to run the model (analysis, Table S6) 
replicates <- 100 # run several model replicates

lat <- matrix(data=NA, nrow = round(dim(species.occ)[1]*(prop/100),0), ncol= replicates); colnames(lat) <- seq(1,replicates,1)
long <- matrix(data=NA, nrow = round(dim(species.occ)[1]*(prop/100),0), ncol= replicates); colnames(long) <- seq(1,replicates,1)
eval_MESS_model <- matrix(data=NA, nrow = replicates, ncol= 6); colnames(eval_MESS_model) <- c("strict_extrapo","MESS_suitable","MESS_unsuitable","suitable_area","unsuitable_area","suitable_area_extrapolated")
table_mess_amelio <- matrix(data=NA, nrow = replicates, ncol= nlayers(pred)); colnames(table_mess_amelio) <- names(pred)
stack.pred <- subset(pred,1); values(stack.pred) <- NA
eval <- matrix(data=NA, nrow = 15, ncol= replicates); colnames(eval) <- seq(1,replicates,1)

for (i in 1:replicates) {
  
      #----------------
      # DATASET
      #----------------
      # sample a part of the occurrence database (is activated when prop <100%, modified above)
      sequence <- seq (1,dim(species.occ)[1],1)
      ech <- sample(sequence,round(dim(species.occ)[1]*(prop/100),0), replace=F)
      species.occ.portion <- species.occ[ech,]
      head(species.occ.portion)

      # Record the position of sampled occurrences
      lat[,i] <- species.occ.portion[,2]
      long[,i] <- species.occ.portion[,1]
      #write.csv(lat,paste("results/Acodontaster/analyse1/lat_",prop,".csv"), append=T)
      #write.csv(long,paste("results/Acodontaster/analyse1/long_",prop,".csv"), append=T)
      
      #----------------
      # Measure of extrapolation score= MESS
      #-----------------
      envi_MESS <- extract(pred, unique(species.occ.portion)) 
      head(envi_MESS)
      x <- dismo::mess(pred,envi_MESS)
      
      y <- x; values(y)<- values(x)>0
      y <- reclassify(y,cbind(FALSE,0)) # MODEL EXTRAPOLATES
      y <- reclassify(y,cbind(TRUE,1)) # MODEL DOES NOT EXTRAPOLATE
      
      ybis <- mask(y,subset(pred,1))
      
      # Calculation of the proportion of the area where extrapolation occurs
      MESS<- reclassify(ybis,cbind(1,NA))
      eval_MESS_model[i,1] <- length(which(!is.na(values(MESS))))*100 /length(which(!is.na(values(subset(pred,1       )))))

      #write.csv(eval_MESS_model,paste("results/Acodontaster/analyse1/eval_MESS_model_",prop,".csv"), append=T)
      #writeRaster(ybis,paste("results/Acodontaster/analyse1/maps/MESS_Acodontaster",i,"_",prop,".asc")) # record MESS LAYER
      
      #----------------
      ## "Improved" MESS: define the contribution of each variable to extrapolation
      #----------------
      stack_MESS_amelio <- subset(pred,1); values(stack_MESS_amelio) <- NA # create an empty raster to initiate a Rasterstack

      # Loop to calculate the value of dissimilarity of each environmental descriptor
      # For each pixel, it will be determined if extrapolation occurs for each environmental descriptor
      for (k in 1:nlayers(pred)){
        presvals <- extract(subset(pred, k), unique(species.occ.portion))  
        x_amelio <- dismo::mess(subset(pred, k),presvals)
        stack_MESS_amelio <- stack(stack_MESS_amelio,x_amelio) 
        }
      
      stack_MESS_amelio <- dropLayer(stack_MESS_amelio,1) # delete the first layer of the stack that was empty (initialization)
      names(stack_MESS_amelio) <- names(pred)
      
      # Search for the environmental layer that is responsible for the lower MESS score (i.e. responsible for the extrapolation)
      MESS_amelio <- which.min(stack_MESS_amelio)
      MESS_amelio <- mask(MESS_amelio, MESS) # indicate "improved" MESS solely in areas where extrapolation occurs 
      #plot(MESS_amelio)
      
      # Calculate the contribution of each environmental descriptor in MESS
      for (k2 in 1:nlayers(pred)){
        table_mess_amelio[i,k2] <- length(which(values(MESS_amelio)==k2))*100 /length(which(!is.na(values(subset(pred,1)))))
        }
      #writeRaster(MESS_amelio,paste("results/Acodontaster/analyse1/maps/MESS_amelio_Acodontaster",i,"_",prop,".asc"))
      write.csv(table_mess_amelio,paste("results/Acodontaster/analyse1/table_mess_amelio_",prop,".csv"), append=T)
      
      ########################################################################################################
      # LAUNCH THE SDM
      #------------------
      source("scripts/clock6.R") #open spatial cross-validation script
      sdm <- clock6(predictors=pred, data=species.occ.portion, KDE_layer=KDE,nb_replicates=1)

      # Open SDM results, see
      #sdm$stack
      #sdm$pdp
      #sdm$eval.file
      #sdm$param.contri
      
      stack.pred <- stack(stack.pred, sdm$stack)
      write.table(sdm$param.contri, "results/Acodontaster/analyse1/influence_param.csv",append=T)
      write.table(sdm$pdp,paste("results/Acodontaster/analyse1/partial/dataframe.partial",i,".csv"))
      eval[,i] <- sdm$eval.file ; rownames(eval)<- rownames(sdm$eval.file)
      #writeRaster(sdm$stack,paste("results/Acodontaster/analyse1/maps/SDMmap_Acodontaster",i,".asc"))
      
      ########################################################################################################
      # overlap of the MESS extrapolation with suitable areas (average map)
      #------------------
      
      # Define suitable areas (MAXSSS > x)
      #-----------------------------
      maxSSS_value <- sdm$eval.file[4,1]
      sdm_pred_favo <- sdm$stack > maxSSS_value # 1: suitable (TRUE); 0: non suitable (FALSE)
      sdm_pred_favo <- reclassify(sdm_pred_favo,cbind(FALSE,0))
      sdm_pred_favo <- reclassify(sdm_pred_favo,cbind(TRUE,1))
      eval_MESS_model[i,4] <- 100*length(which(values(sdm_pred_favo)==1))/length(which(!is.na(values(subset(pred,1))))) # proportion of the projection area that is predicted suitable 
      eval_MESS_model[i,5] <- 100*length(which(values(sdm_pred_favo)==0))/length(which(!is.na(values(subset(pred,1))))) # proportion of the projection area that is predicted unsuitable 
      #eval_MESS_model[i,10] <- length(which(values(sdm_pred_favo)==1)) # number of suitable pixels 
      #eval_MESS_model[i,11] <- length(which(values(sdm_pred_favo)==0)) # number of unsuitable pixels 
      
      # Overlap with the MESS
      #----------------------
      # area where extrapolation occurs (MESS=0) and suitable area (sdm_pred_favo=1)
      suitable_extrapo_raster <- sdm_pred_favo==1 & ybis==0 # value 1: TRUE: area where extrapolation occurs and SDM is suitable
      unsuitable_extrapo_raster <- sdm_pred_favo==0 & ybis==0 # value 1: TRUE: area where extrapolation occurs and SDM is UNsuitable
      eval_MESS_model[i,2] <- 100*length(which(values(suitable_extrapo_raster)==1))/ length(na.omit(which(values(ybis)==0))) # proportion of the extrapolation area where SDM is predicted suitable
      eval_MESS_model[i,3] <- 100*length(which(values(unsuitable_extrapo_raster)==1))/ length(na.omit(which(values(ybis)==0))) # proportion of the extrapolation area where SDM is predicted UNsuitable
      eval_MESS_model[i,6] <- 100*length(which(values(suitable_extrapo_raster)==1))/ length(which(values(sdm_pred_favo)==1)) # proportion of the area where SDM is predicted suitable for which extrapolation occurs 
      
      write.csv(eval_MESS_model,paste("results/Acodontaster/analyse1/eval_MESS_model_",prop,".csv"), append=T)
      
} # end of the loop of replicates

#  Save final results
write.csv(lat,paste("results/Acodontaster/analyse1/lat_",prop,"ALL.csv"))
write.csv(long,paste("results/Acodontaster/analyse1/long_",prop,"ALL.csv"))
write.csv(eval_MESS_model,paste("results/Acodontaster/analyse1/eval_MESS_model_",prop,"ALL.csv"))
write.csv(table_mess_amelio,paste("results/Acodontaster/analyse1/table_mess_amelio_",prop,"ALL.csv"))
write.csv(eval,paste("results/Acodontaster/analyse1/eval.file.SDM_",prop,"ALL.csv"))
#writeRaster(ybis,paste("results/Acodontaster/complete_analysis100/MESS_Acodontaster",prop,".asc")) # save a layer of MESS 

# SDM outputs (maps)
stack.pred <- dropLayer(stack.pred,1)
writeRaster(stack.pred,paste("results/Acodontaster/analyse1/stack.pred",prop,".grd"))
mean_SDM <- calc(stack.pred, mean,na.rm=T)
writeRaster(mean_SDM,paste("results/Acodontaster/analyse1/mean_SDM",prop,".asc"))
sd_SDM <- calc(stack.pred, sd,na.rm=T)
writeRaster(sd_SDM,paste("results/Acodontaster/analyse1/sd_SDM",prop,".asc"))

#REFERENCES 
# De Broyer, C., Koubbi, P., Griffiths, H.J., Raymond, B., d’Udekem d’Acoz, C., Van de Putte, A.P., … Ropert-Coudert, Y. (2014). Biogeographic atlas of the Southern Ocean (p. 498). C. De Broyer, & P. Koubbi (Eds.). Cambridge: Scientific Committee on Antarctic Research.

# Guillaumot, C., Artois, J., Saucède, T., Demoustier, L., Moreau, C., Eléaume, M. ... & Danis, B. (2019). Broad-scale species distribution models applied to data-poor areas. Progress in oceanography, 175, 198-207.