#------------------------------------------------------
# R scripts associated with the article
# Simple or hybrid? Next generation ecological models to study the distribution of Southern Ocean marine species
# by Guillaumot Charlène, Buba Yehezkel, Belmaker Jonathan, Fourcy Damien, Dubois Philippe, Danis Bruno, Saucède Thomas
# updated January 2021
#
# File: Spatial projection of the DEB
#------------------------------------------------------

library(raster)
library(R.matlab)
source("scripts/get_DEB_pars.R")
source("scripts/get_tj2.R")
source("scripts/t_corr.R")
source("scripts/ode_functs.R")
source("scripts/get_powers_j.R")
library(raster)
library(rgdal)
library(RColorBrewer)
library(scales)

# Set ups 
#----------
#my.palette <- brewer.pal(n = 9, name = "Oranges")
#my.palette.blue <- rev(brewer.pal(n = 9, name = "Blues"))
#palettecol <-colorRampPalette(c("deepskyblue", "darkseagreen","lightgreen","green","yellow","gold","orange", "red","firebrick"))(100)
#ker_map <- readOGR("envi/TC_Kerguelen_WGS84_EPSG4326_polygone.shp")

# open environmental layers
#---------------------------------------------------------------
Bathy <- raster("envi/bathymetry_Morbihan_crop_finale.asc")

library(raster)
depth <- raster("envi/bathymetry_Morbihan_crop_finale.asc")
f_09022017 <- raster("envi/f_layer_09022017_OK2.asc"); f_09022017 <- mask(crop(f_09022017, depth), depth)
# for this layer, artefact due to clouds -> replace values = 0 to NA 
f_09022017[f_09022017==0]<-NA

f_20082017 <- raster("envi/f_layer_20082017_OK2.asc"); f_20082017 <- mask(crop(f_20082017, depth), depth)
temp_09022017 <- raster("envi/SST_MUR09022017_OK2.asc")
temp_20082017 <- raster("envi/SST_MUR20082017_OK2.asc")

flayer <- f_09022017
temp <- temp_09022017

predictors <- stack(depth, temp,flayer)
names(predictors) <- c("depth","temp","f")

########### open Matlab file
#-------------------------------------------------
result_abatus <- get_DEB_pars("scripts_DEB/results_Abatus_cordatus.mat")

# convert each layer in a global table, with lat long T° and f values 
#---------------------------------------------------------------------
temp2 <- reclassify(temp, cbind(NA, 9999))
temp_v <- rasterToPoints(temp2)
temp_NA <- temp_v[,3]
temp_NA <- replace(temp_NA,temp_NA==9999, 55)
temp_v2 <- cbind(temp_v[,c(1,2)],temp_NA)
head(temp_v2)
colnames(temp_v2) <- c("longitude", "latitude","T")
head(temp_v2)
tail(temp_v2)

f.layer2 <- reclassify(flayer, cbind(NA, 1)) # all NA values (land) are replaced by 1 because the function cannot handle NA. No consequences on results beaucoup then all maps are masked by the f.layer raster layer (with replacement of the 1 by NA again)
f.layer_v <- rasterToPoints(f.layer2)
f.layer_NA <- f.layer_v[,3]
f.layer_NA <- replace(f.layer_NA,f.layer_NA==0, 0.01)
f.layer_v2 <- cbind(f.layer_v[,c(1,2)],f.layer_NA)
colnames(f.layer_v2) <- c("longitude", "latitude","f")
head(f.layer_v2)

tableTf <- cbind(temp_v2,f.layer_v2[,3])
colnames(tableTf)<- c("longitude", "latitude","T","f")
head(tableTf) 

# Add to the table the Arrhenis correction that varies according to the temperature contained in each pixel
#----------------------------------------------------------------------------------------------------------
# LOOP
TC_table <- matrix(nrow=nrow(temp_v2),ncol=1,data=NA)

temp_v3 <- as.data.frame(temp_v2)
temp_v3$T [is.na(temp_v3$T)] <- 0

for (i in 1:nrow(temp_v3)){
  TC_table[i,1] <- t_corr(p=result_abatus, Ti=(273.15+temp_v3[i,3]))
}
head(TC_table)

tableTf <- cbind(tableTf,TC_table)
colnames(tableTf)<- c("longitude", "latitude","T","f", "TC")
head(tableTf)

#----------------------------------------------------------
# PROJECT the physiological performance in space 
#----------------------------------------------------------
                                            #........................#
                                            # Project Ultimate size  #
                                            #........................#
library(deSolve)
library(pracma)

get_tj.li <- matrix(nrow=nrow(tableTf),ncol=1,data=NA)

 for (i in 1:nrow(tableTf)){
   get_tj.li[i,1] <- get_tj(p=result_abatus,f=tableTf[i,4])$li
 }

 L_i <- result_abatus$L_m * get_tj.li
 Lw_i <- L_i* result_abatus$del_M

# Layer ultimate size
rastervide <- temp; values(rastervide) <- 0
layer.ultimate.size <- rasterize(tableTf[,c(1,2)], rastervide,Lw_i)
layer.ultimate.size <- mask(layer.ultimate.size, flayer)
layer.ultimate.size <- mask(layer.ultimate.size, temp)

plot(layer.ultimate.size, main="Ultimate size (cm)")
plot(ker_map, add=T)
writeRaster(layer.ultimate.size, "layer.ultimate.size.asc")

                      #...............................................................#
                      # Project the ability to reproduce and Reproduction performance #
                      #...............................................................#
# If lp> li, the organism is not able to reproduce !
# If tp> ti, the organism die without reproducing!

get_tj.lp <- matrix(nrow=nrow(tableTf),ncol=1,data=NA)
get_tj.li <- matrix(nrow=nrow(tableTf),ncol=1,data=NA)
#get_tj.tp <- matrix(nrow=nrow(tableTf),ncol=1,data=NA)
#get_tj.ti <- matrix(nrow=nrow(tableTf),ncol=1,data=NA)

#Be careful, loops take a very long time to run (the more f values are small, the more computing time for the get.tj function)
 for (i in 1:nrow(tableTf)) {
    getx <- get_tj(p=result_abatus,f=tableTf[i,4])
    get_tj.lp[i,1] <- getx$lp
    get_tj.li[i,1] <- getx$li
    #get_tj.tp[i,1] <- getx$tp
    #get_tj.ti[i,1] <- getx$ti
 }

# Layers lp li tp ti
#---------------------
rastervide <- temp; values(rastervide) <- 0

layer.lp <- rasterize(tableTf[,c(1,2)], rastervide,get_tj.lp)
layer.li <- rasterize(tableTf[,c(1,2)], rastervide,get_tj.li)
#layer.tp <- rasterize(tableTf[,c(1,2)], rastervide,get_tj.tp)
#layer.ti <- rasterize(tableTf[,c(1,2)], rastervide,get_tj.ti)

# No reproduction plots
# If lp> li, the organism is not able to reproduce !
# If tp> ti, the organism die without reproducing!

#plot(layer.lp>layer.li, col=c("green","black")) # Black: TRUE: 1 lp>li no reproduction
#plot(layer.tp>layer.ti, col=c("black","green"))
#plot(layer.tp)

#writeRaster(layer.lp, "layer.lp.asc")
#writeRaster(layer.li, "layer.li.asc")
writeRaster(layer.lp>layer.li, "repro_DEB_aug.asc") # Black: TRUE: 1 lp>li no reproduction; Green= 0, False, reproduction
#writeRaster(layer.tp, "results/layer.tp.asc")
#writeRaster(layer.ti, "results/layer.ti.asc")

                                    #----------------------------
                                    # Evaluate if PM>PC adults  
                                    #----------------------------
Lw <- 3
L <- Lw*result_abatus$del_M # structural length (cm), based on a standard individual of 4 cm
pa <- matrix(nrow=nrow(tableTf),ncol=1,data=NA) # assimilation
pc <- matrix(nrow=nrow(tableTf),ncol=1,data=NA) # into reserves
pm <- matrix(nrow=nrow(tableTf),ncol=1,data=NA) # somatic maintenance
pj <- matrix(nrow=nrow(tableTf),ncol=1,data=NA) # maturity maintenance
pg <- matrix(nrow=nrow(tableTf),ncol=1,data=NA) # growth

# if the energy required for somatic maintenance (pM) is higher than the energy contained into the reserve compartment (pc), the individual is supposed to die 

# Calculate these 'powers'
for (i in 1:nrow(tableTf)){
  fi <- tableTf[i,4]
  s_M <- get_tj(p=result_abatus,f=tableTf[i,4])$s_M
  powers <- get_powers_j(p=result_abatus,L,e=fi, s_M)
  pa[i,1] <- powers[1]  / result_abatus$k_M / tableTf[i,5]
  pc[i,1] <- powers[2]  / result_abatus$k_M / tableTf[i,5]
  pm[i,1] <- powers[3]  / result_abatus$k_M / tableTf[i,5]
  pj[i,1] <- powers[4]  / result_abatus$k_M / tableTf[i,5]
  pg[i,1] <- powers[5]  / result_abatus$k_M / tableTf[i,5]
}

# layers
#---------------------
rastervide <- temp; values(rastervide) <- 0
layer.pa <- rasterize(tableTf[,c(1,2)], rastervide,pa)
#layer.pa <- mask(layer.pa, chla)
#plot(layer.pa)
writeRaster(layer.pa, "layer.pa_feb.asc")

layer.pc <- rasterize(tableTf[,c(1,2)], rastervide,pc)
#layer.pc <- mask(layer.pc, chla)
#plot(layer.pc)
writeRaster(layer.pc, "layer.pc_feb.asc")

layer.pm <- rasterize(tableTf[,c(1,2)], rastervide,pm)
#layer.pm <- mask(layer.pm, chla)
#plot(layer.pm)
writeRaster(layer.pm, "layer.pm_feb.asc")

layer.pj <- rasterize(tableTf[,c(1,2)], rastervide,pj)
#layer.pj <- mask(layer.pj, chla)
#plot(layer.pj)
writeRaster(layer.pj, "layer.pj_feb.asc")

layer.pg <- rasterize(tableTf[,c(1,2)], rastervide,pg)
#layer.pg <- mask(layer.pg, chla)
#plot(layer.pg)
writeRaster(layer.pg, "layer.pg_feb.asc")

## survival?
#plot(layer.pm>layer.pc, col=c("black","green")) # green: survival; black: death

