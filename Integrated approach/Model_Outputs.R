#------------------------------------------------------
# R scripts associated with the article
# Simple or hybrid? Next generation ecological models to study the distribution of Southern Ocean marine species
# by Guillaumot Charlène, Buba Yehezkel, Belmaker Jonathan, Fourcy Damien, Dubois Philippe, Danis Bruno, Saucède Thomas
# updated January 2021
#
# File: All model outputs that have been produced for the article

#--------------------------------------------------------------------------------------------------------

                                                          #------------------#
                                                          # POSTERIOR PRIORS #
                                                          #------------------#
# For all codes, when reaching the step
i=1
coef_b[min(which(is.na(coef_b))):(i*4000),] <- as.matrix(m1Results)

# the histograms for the posterior coefficients can be saved/plotted with : 
save_path <- "results/post_coeff_sdm_simple/"
col_name <- n
write.table(coef_b, paste(save_path,col_name,"_coef_b.csv",sep = ""))

save_path <- "results/"
for (i in 1:ncol(coef_b)){
  col_name <- colnames(coef_b)[i]
  png(filename = paste(save_path,col_name,"_hist_regular.png",sep = ""))
  hist(coef_b[,i], xlab = col_name, col = "blue", breaks = 50, main = paste("Histogram for coefficient",col_name))
  dev.off()
}

                                                #--------------------------#
                                                # PARTIAL DEPENDENCE PLOTS #
                                                #--------------------------#
library(rgdal)
library(RColorBrewer)
ker_map <- readOGR("envi/TC_Kerguelen_WGS84_EPSG4326_polygone.shp") # precise shapefile of Kerguelen
my.palette <- brewer.pal(n = 9, name = "Oranges")
#my.palette.blue <- rev(brewer.pal(n = 9, name = "Blues"))
palettecol <-colorRampPalette(c("deepskyblue", "darkseagreen","lightgreen","green","yellow","gold","orange", "red","firebrick"))(100)

# OPEN ENVIRONMENTAL LAYERS
#--------------------------
library(raster)
depth <- raster("envi/bathymetry_Morbihan_crop_finale.asc")
f_09022017 <- raster("envi/f_layer_09022017_OK2.asc"); f_09022017 <- mask(crop(f_09022017, depth), depth)
# for this layer, artefact due to clouds -> replace values = 0 to NA 
f_09022017[f_09022017==0]<-NA

f_20082017 <- raster("envi/f_layer_20082017_OK2.asc"); f_20082017 <- mask(crop(f_20082017, depth), depth)
temp_09022017 <- raster("envi/SST_MUR09022017_OK2.asc")
temp_20082017 <- raster("envi/SST_MUR20082017_OK2.asc")

f_layer_feb <- f_09022017
f_layer_aug <- f_20082017

temp_layer_feb <- temp_09022017
temp_layer_aug <- temp_20082017


# OPEN RESULTS 
#----------------
# Simple SDM
#-------------
sdm_feb_mean
sdm_aug_mean

# DEB within SDM
#---------------
DEBwSDM_feb_mean
DEBwSDM_aug_mean

# Bayesian integrated model 
#--------------------------
mean_combined_feb <- BAYESI_feb_mean
mean_combined_aug <- BAYESI_aug_mean

# Properly adjust the layers
temp_layer_feb <- crop(temp_layer_feb,extent(mean_combined_feb))
temp_layer_feb <- mask(temp_layer_feb,(mean_combined_feb))

temp_layer_aug <- crop(temp_layer_aug,extent(mean_combined_aug))
temp_layer_aug <- mask(temp_layer_aug,(mean_combined_aug))

f_layer_feb <- crop(f_layer_feb,extent(mean_combined_feb))
f_layer_feb <- mask(f_layer_feb,(mean_combined_feb))

f_layer_aug <- crop(f_layer_aug,extent(mean_combined_aug))
f_layer_aug <- mask(f_layer_aug,(mean_combined_aug))

depth1 <- crop(depth,extent(mean_combined_feb))
depth1 <- mask(depth1,(mean_combined_feb))
depth2 <- crop(depth,extent(mean_combined_aug))
depth2 <- mask(depth2,(mean_combined_aug))

# partial dependence plots 
#-------------------------
sdm_feb_mean_v <- rasterToPoints(sdm_feb_mean)
sdm_aug_mean_v <- rasterToPoints(sdm_aug_mean)
DEBwSDM_feb_mean_v <- rasterToPoints(DEBwSDM_feb_mean)
DEBwSDM_aug_mean_v <- rasterToPoints(DEBwSDM_aug_mean)
mean_combined_feb_v <- rasterToPoints(mean_combined_feb)
mean_combined_aug_v <- rasterToPoints(mean_combined_aug)

f_layer_feb_v <- rasterToPoints(f_layer_feb)
temp_layer_feb_v <- rasterToPoints(temp_layer_feb)
depth1_v <- rasterToPoints(depth1)
depth2_v <- rasterToPoints(depth2)
temp_layer_aug_v <- rasterToPoints(temp_layer_aug)
f_layer_aug_v <- rasterToPoints(f_layer_aug)


partial_temp_feb  <- data.frame(cbind(x_temp=temp_layer_feb_v[,3],  y_prob1 = sdm_feb_mean_v[,3], y_probC1 = DEBwSDM_feb_mean_v[,3],  y_probC2 = mean_combined_feb_v[,3]))
write.csv(partial_temp_feb,'partial_temp_feb.csv')

partial_temp_aug  <- data.frame(cbind(x_temp=temp_layer_aug_v[,3],  y_prob2 = sdm_aug_mean_v[,3], y_probC1 = DEBwSDM_aug_mean_v[,3],   y_probC2 = mean_combined_aug_v[,3]))


partial_f_feb  <- data.frame(cbind(x_f=f_layer_feb_v[,3],  y_prob1 = sdm_feb_mean_v[,3],  y_probC1 = DEBwSDM_feb_mean_v[,3], y_probC2 = mean_combined_feb_v[,3]))
write.csv(partial_f_feb,'partial_f_feb.csv')

partial_f_aug  <- data.frame(cbind(x_f=f_layer_aug_v[,3],  y_prob2 = sdm_aug_mean_v[,3],  y_probC1 = DEBwSDM_aug_mean_v[,3], y_probC2 = mean_combined_aug_v[,3]))


partial_depth1  <- data.frame(cbind(x_depth=depth1_v[,3],  y_prob1 = sdm_feb_mean_v[,3],  y_probC1 = DEBwSDM_feb_mean_v[,3], y_probC2 = mean_combined_feb_v[,3]))
write.csv(partial_depth1,'partial_depth1.csv')

partial_depth2  <- data.frame(cbind(x_depth=depth2_v[,3],    y_prob2 = sdm_aug_mean_v[,3], y_probC1 = DEBwSDM_aug_mean_v[,3],  y_probC2 = mean_combined_aug_v[,3]))


library(reshape)
partial_temp_feb_resh <- melt(partial_temp_feb, id.vars = "x_temp")
partial_temp_aug_resh <- melt(partial_temp_aug, id.vars = "x_temp")

partial_f_feb_resh <- melt(partial_f_feb, id.vars = "x_f")
partial_f_aug_resh <- melt(partial_f_aug, id.vars = "x_f")

partial_depth_resh1 <- melt(partial_depth1, id.vars = "x_depth")
partial_depth_resh2 <- melt(partial_depth2, id.vars = "x_depth")


library(ggplot2)
library(gridExtra)

p1 <- ggplot(partial_temp_feb_resh, aes(x=x_temp, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Temperature (°C)")+ 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ 
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)

ggplot(partial_temp_aug_resh, aes(x=x_temp, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Temperature (°C)")+ 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)

p2 <- ggplot(partial_f_feb_resh, aes(x=x_f, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Food availability")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)

ggplot(partial_f_aug_resh, aes(x=x_f, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Food availability")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)


p3 <- ggplot(partial_depth_resh1, aes(x=x_depth, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Depth (m)")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)

ggplot(partial_depth_resh2, aes(x=x_depth, y=value, colour=variable)) +
  #geom_point()+
  geom_smooth(lwd=1.5)+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Depth (m)")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text = element_text( size=14), axis.title=element_text(size=12))+
  ylim(0,1)

grid.arrange(p1, p2, p3, ncol=2, nrow = 2)


mean_august1_v <- rasterToPoints(mean_august1)
f_layer_v <- rasterToPoints(f_layer)
temp_layer_v <- rasterToPoints(temp_layer)
depth_v <- rasterToPoints(depth)

partial <- as.data.frame(cbind(x_temp=temp_layer_v[,3], x_f=f_layer_v[,3], x_depth=depth_v[,3],  y_prob = mean_august1_v[,3]))

library(ggplot2)
library(gridExtra)

p1 <- ggplot(partial, aes(x=x_temp, y=y_prob)) +
  geom_point(colour="lightgrey")+
  geom_smooth()+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Temperature (°C)")

p2 <- ggplot(partial, aes(x=x_f, y=y_prob)) +
  geom_point(colour="lightgrey")+
  geom_smooth()+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Food availability")

p3 <- ggplot(partial, aes(x=x_depth, y=y_prob)) +
  geom_point(colour="lightgrey")+
  geom_smooth()+
  theme_bw()+
  ylab("Density function of marginal effect") + xlab("Depth (m)")

grid.arrange(p1, p2, p3, ncol=2, nrow = 2)

                                                #--------------------------#
                                                # CORRECTLY CLASSIFIED DATA #
                                                #--------------------------#
# sdm_feb_mean; DEBwSDM_feb_mean; mean_combined_feb are the raster files (maps of predictions) for the simple, DEB within SDM and Bayesian integrated models. They are the average of the 50 model replicates

abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]
proba_pres_SDM_feb <- raster::extract(sdm_feb_mean, abatus.occ); mean(proba_pres_SDM_feb, na.rm=T)
proba_pres_DEBwSDM_feb <- raster::extract(DEBwSDM_feb_mean, abatus.occ); mean(proba_pres_DEBwSDM_feb, na.rm=T)
proba_pres_BAY_feb <- raster::extract(mean_combined_feb, abatus.occ); mean(proba_pres_BAY_feb, na.rm=T)

proba_pres_SDM_aug <- raster::extract(sdm_aug_mean, abatus.occ); mean(proba_pres_SDM_aug, na.rm=T)
proba_pres_DEBwSDM_aug <- raster::extract(DEBwSDM_aug_mean, abatus.occ); mean(proba_pres_DEBwSDM_aug, na.rm=T)
proba_pres_BAY_aug <- raster::extract(mean_combined_aug, abatus.occ); mean(proba_pres_BAY_aug, na.rm=T)

                                                #--------------------------#
                                                # CALCULATE EXTRAPOLATION #
                                                #--------------------------#
# Multivariate Environmental Similarity Surface (Elith et al. 2010)
library(raster)
abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]
depth <- raster("envi/bathymetry_Morbihan_crop_finale.asc")
f_09022017 <- raster("envi/f_layer_09022017_OK2.asc"); f_09022017 <- mask(crop(f_09022017, depth), depth)
# for this layer, artefact due to clouds -> replace values = 0 to NA 
f_09022017[f_09022017==0]<-NA

f_20082017 <- raster("envi/f_layer_20082017_OK2.asc"); f_20082017 <- mask(crop(f_20082017, depth), depth)
temp_09022017 <- raster("envi/SST_MUR09022017_OK2.asc")
temp_20082017 <- raster("envi/SST_MUR20082017_OK2.asc")

f_layer <- f_09022017
temp_layer <- temp_09022017
DEB_out <- raster("results/layer.pm_feb.asc"); DEB_out <- mask(DEB_out, temp_09022017)

predictors <- stack(depth, temp_layer,f_layer, DEB_out)


envi.presences <- unique(extract(predictors,abatus.occ))
x <- dismo::mess(predictors, na.omit(envi.presences))

y <- x; values(y)<- values(x)>0  # refers to Elith et al. (2010): when the calculated MESS values are negative, it means that it is extrapolating (outside of boundaries)
y <- reclassify(y,cbind(FALSE,0)) # extrapolation area (black)
y <- reclassify(y,cbind(TRUE,1))  # non extrapolation, inside the boundaries of calibration

plot(y)
writeRaster(y,"results/MESS_Abatus.asc")
MESS_Abatus <- raster("results/MESS_Abatus.asc")

# MESS WITH DETAILS (which layer responsible for extrapolation?)
stack_MESS_dissim <- subset(predictors,1); values(stack_MESS_dissim) <- NA

# calculer la valeur de la dissimilarité pour chacune des couches
for (i in 1:nlayers(predictors)){
  presvals <- extract(subset(predictors, i), abatus.occ)
  x <- dismo::mess(subset(predictors, i),presvals)
  stack_MESS_dissim <- stack(stack_MESS_dissim,x)
}

stack_MESS_dissim <- dropLayer(stack_MESS_dissim,c(1))
names(stack_MESS_dissim) <- c("depth","temp","f", "DEB layer")

# which layer contains the smallest MESS value?
MESS_amelio <- which.min(stack_MESS_dissim)
MESS_amelio

MESS_global_mask<- y
MESS_global_mask <- reclassify(MESS_global_mask,cbind(1,NA))
MESS_amelio_masked <- mask(MESS_amelio,MESS_global_mask)
plot(MESS_amelio_masked)
writeRaster(MESS_amelio_masked,"results/MESS_ameliore_Abatus_DEBpc_pmpJ_aug.asc")

library(rgdal)
library(RColorBrewer)
ker_map <- readOGR("envi/TC_Kerguelen_WGS84_EPSG4326_polygone.shp") # nouvelle couche shapefile précise de Ker

plot(MESS_amelio_masked, col=c("black","grey70","orange","green"))
plot(ker_map, add=T)

                                                    #---------------#
                                                    # CALCULATE AUC #
                                                    #---------------#

library(raster)
library(ROCR)

#open one output layer (e.g. replicate #50)
p <- raster("results/SDM_combined/50 replicats/result_combined_sdm_mean_august50.asc")
background_data <- read.csv("results/SDM_combined/50 replicats/background_data_aug_50.csv", header=T, sep=",")[,c(2,3)]

p <- raster("results/SDM_combined/50 replicats/result_combined_sdm_mean_feb30.asc")
background_data <- read.csv("results/SDM_combined/50 replicats/background_data_feb_30.csv", header=T, sep=",")[,c(2,3)]
abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]

predictions <- c(extract(p,background_data),extract(p, abatus.occ))
labels <- c(rep(0,nrow(background_data)),rep(1, nrow(abatus.occ)))
excluded <- c(-which(is.na(predictions)))
predictions <- predictions[-which(is.na(predictions))]
labels <-labels[excluded]

pred_ROCR <- prediction(predictions, labels)
roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
#plot(roc_ROCR, main = "ROC curve", colorize = T)
#abline(a = 0, b = 1)
auc_ROCR <- performance(pred_ROCR, measure = "auc")
auc_ROCR@y.values[[1]]

write.csv(auc_ROCR@y.values[[1]],"result_auc_aug_8.csv")
#apply(auc_ROCR_mat, mean, MARGIN=2)
#apply(auc_ROCR_mat, sd, MARGIN=2)

                                                    #------------------#
                                                    # CALCULATE MAXSSS #
                                                    #------------------#
library(ENMeval)
library(dismo)

abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]
backgr<- read.csv("results/SDM_combined/50_rep/background_data_aug_50.csv", header=T, sep=",")[,-1]
pred <- raster("results/SDM_combined/50_rep/result_combined_sdm_mean_august50.asc")
plot(pred)
testp <- extract(pred,abatus.occ)
testa <- extract(pred,backgr)

eval.data <- dismo::evaluate(p = testp, a = testa)
tab <- base::cbind(eval.data@t, eval.data@TPR + eval.data@TNR)
maxSSS <- (base::subset(tab, tab[, 2] == max(tab[, 2])))[1, 1]
maxSSS

100*length(which(extract(pred,abatus.occ)>maxSSS))/nrow(abatus.occ)

