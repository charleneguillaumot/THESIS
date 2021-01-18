#------------------------------------------------------
# R scripts associated with the article
# Simple or hybrid? Next generation ecological models to study the distribution of Southern Ocean marine species
# by Guillaumot Charlène, Buba Yehezkel, Belmaker Jonathan, Fourcy Damien, Dubois Philippe, Danis Bruno, Saucède Thomas
# updated January 2021
#
# File: procedure for the 'DEB-within-SDM' approach
#------------------------------------------------------
#--------------------------
# Upload librairies
#--------------------------
library(nlme)
library(minpack.lm)
library(AICcmodavg)
library(sp)
library(raster)
#library(BIOMOD)
library(stats4)
library(caTools)
library(dismo)
library(spaa)
#library(mopa)
#library(SDMTools)
source("scripts/ecospat.R")
source("scripts/ex1_globals.r")

#--------------------------
# OPEN OCCURRENCE RECORDS
#--------------------------
abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]
head(abatus.occ)

sp.occr_loc=data.frame(lon=abatus.occ$decimalLongitude,lat=abatus.occ$decimalLatitude)
sp_name=data.frame(Occurrence=rep(1,nrow(sp.occr_loc)))
species_occ=SpatialPointsDataFrame(coords=sp.occr_loc,data=sp_name)

#--------------------------
# OPEN ENVIRONMENTAL LAYERS
#--------------------------
depth <- raster("envi/bathymetry_Morbihan_crop_finale.asc")
f_09022017 <- raster("envi/f_layer_09022017_OK2.asc"); f_09022017 <- mask(crop(f_09022017, depth), depth)
# for this layer, artefact due to clouds -> replace values = 0 to NA 
f_09022017[f_09022017==0]<-NA

f_20082017 <- raster("envi/f_layer_20082017_OK2.asc"); f_20082017 <- mask(crop(f_20082017, depth), depth)
temp_09022017 <- raster("envi/SST_MUR09022017_OK2.asc")
temp_20082017 <- raster("envi/SST_MUR20082017_OK2.asc")

DEB_out <- raster("result_Topen/pc_pmpJ_aug.asc"); DEB_out <- mask(DEB_out, temp_09022017) # layer created by spatially projecting the DEB model (see other R file -> Spatial_DEB_XX.R)

f_layer <- f_20082017  # to be changed according to the studied period (August/February)
temp_layer <- temp_20082017

predictors <- stack(depth, temp_layer,f_layer, DEB_out)
names(predictors) <- c("depth","temp","f","DEB_out")

#----------------------------------------------------------------------------------------------
### LAUNCH MODEL ###
#----------------------------------------------------------------------------------------------

# Background data sampling 
#--------------------------
background_data <- xyFromCell(predictors, sample(which(!is.na(values(depth))), 200))

# different calculations on environmental descriptors
#----------------------------------------------------
predictors$temp_sq <- predictors$temp^2
predictors$f_sq <- predictors$f^2

### EXTRACT data
#---------------
pred_pres <- extract(predictors,species_occ)
pred_pres <-cbind(pred_pres, presence=rep(1,nrow(pred_pres)))

pred_abs <- extract(predictors,background_data)
pred_abs <-cbind(pred_abs, presence=rep(0,nrow(pred_abs)))

pred <- as.data.frame(rbind(pred_pres,pred_abs))
pred <- na.omit(pred)
head(pred)

# test for different GLM models and choose the best one 
#------------------------------------------------------
beta_poly_twoa <- glm(presence ~ depth + temp + f + temp_sq + f_sq + DEB_out, family = binomial, data=pred)
aicc_poly_twoa <- AICc(beta_poly_twoa, second.ord = TRUE)

beta_poly_onea <- glm(presence ~ depth + temp + f + DEB_out, family = binomial, data=pred)
aicc_poly_onea <- AICc(beta_poly_onea, second.ord = TRUE)

aicc_poly_p <- data.frame(aicc_poly_twoa,aicc_poly_onea)
aicc_poly_p <- sort(aicc_poly_p,decreasing = FALSE)

bPrior <- matrix(c(
  0.0, 1.0E-4,	# mean and tau for a0
  0.0, 1.0E-4,	#  a1 
  0.0, 1.0E-4,
  0.0, 1.0E-4,
  0.0, 1.0E-4,
  0.0, 1.0E-4,
  0.0, 1.0E-4,
  0.0, 1.0E-4,
  0.0, 1.0E-4),	
  byrow=TRUE, ncol=2)

if (aicc_poly_p[1]==aicc_poly_twoa){
  m1Data <- list(depth = pred$depth, temp = pred$temp, f = pred$f,DEB_out=pred$DEB_out,  presence = pred$presence, bPrior=bPrior)
  m1Pars <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5','b6')
}  else if (aicc_poly_p[1]==aicc_poly_onea){
  m1Data <- list(depth = pred$depth, temp = pred$temp, f = pred$f,DEB_out=pred$DEB_out, presence = pred$presence, bPrior=bPrior,b4=0,b5=0)
  m1Pars <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5','b6')
}

# Empty matrix to fill
coef_b <- matrix(nrow = 4000, ncol = length(m1Pars)); colnames(coef_b) <- m1Pars
colnames(coef_b) <- m1Pars
data_list <- vector(mode = "list")

# one model, one replicate (this script has to be run several times for several replicates)
model1 <- jags.model('scripts/sdm_model_poly_DEB.jags', data = m1Data, n.chains=settings$chains, n.adapt=settings$tuning)
update(model1, settings$burnin)

m1Results <- coda.samples(model1, m1Pars, settings$samples*settings$thin, thin=settings$thin)

#----------------------------------------------------------------------------------------------
### PREDICTIONS ###
#----------------------------------------------------------------------------------------------
# functions 
formula <- function(x,b){
  boot::inv.logit(b[1] + b[2]*x[1] + b[3]*x[2]+ b[4]*x[3] + b[5]*x[2]^2 + b[6]*x[4])}

formula_apply <- function(x,b){
  apply(x,1,function(y)formula(y,b)) 
}

predictions <- function(x,b){
  return(t(pbapply(b,1,function(y)formula_apply(x,y))))
}

i=1
coef_b[min(which(is.na(coef_b))):(i*4000),] <- as.matrix(m1Results)
#coef_b <- coef_b[sample(nrow(coef_b),1000,replace = F),] # random sampling of a part of the results -> for fast trials 

# MAP OF PREDICTIONS 
# #-----------------------------
print("calculate prob_sp")
predictors_v <- as.data.frame(na.omit(rasterToPoints(predictors)))
predictors_data <- predictors_v[,-c(1,2)]
head(predictors_data)

prob_sp <- predictions(predictors_data,coef_b)

prob_sp_mean <- colMeans(prob_sp)
max_prob_sp_mean <- max(prob_sp_mean)
SD_prob_sp <- apply(prob_sp, 2, sd)

p <- rasterFromXYZ(data.frame(x=predictors_v[,1],y=predictors_v[,2], z=prob_sp_mean))
writeRaster(p, "result_sdm_deb_mean_aug_1.asc")
p_sd <- rasterFromXYZ(data.frame(x=predictors_v[,1],y=predictors_v[,2], z=SD_prob_sp))
writeRaster(p_sd, "result_sdm_deb_sd_aug_1.asc")

#----------------------------------------------------------------------------------------------
### GET GLM PARAMETERS ###
#----------------------------------------------------------------------------------------------
mo_data <- as.matrix(m1Results)
head(mo_data)
b0 <- mo_data[,1]
b1 <- mo_data[,2]
b2 <- mo_data[,3]
b3 <- mo_data[,4]
b4 <- mo_data[,5]
b5 <- mo_data[,6]
b6 <- mo_data[,6]

depth_rng <- seq(-195,0, length=500)
temp_rng <- seq(6.62,7.36,length=500)
f_rng <- seq(0,1, length=500)
DEB_out_rng <- seq(0,2500, length=500)

prob <- data.frame(matrix(NA, nrow = length(b0), ncol = length(temp_rng))) 

library(foreach)
foreach (i= 1:length(temp_rng), j = 1:length(b0)) %do% {
  prob[j,i] <- b0[j] +  (b1[j] * depth_rng[i]) + (b2[j] * f_rng[i]) + (b3[j] * temp_rng[i]) + (b4[j] * temp_rng[i]^2) + (b5[j] * f_rng[i]^2) + (b6[j] * DEB_out_rng[i]) 
  prob[j,i] <- boot::inv.logit(prob[j,i])
}

# store the average values of the parameters 
b0_m <- mean(b0)
b1_m <- mean(b1)
b2_m <- mean(b2)
b3_m <- mean(b3)
b4_m <- mean(b4)
b5_m <- mean(b5)
b6_m <- mean(b6)

b0_sd <- sd(b0)
b1_sd <- sd(b1)
b2_sd <- sd(b2)
b3_sd <- sd(b3)
b4_sd <- sd(b4)
b5_sd <- sd(b5)
b6_sd <- sd(b6)

coef_data_sdm <- as.data.frame(rbind(b0_m,b1_m,b2_m,b3_m,b4_m,b5_m,b6_m))
coef_data_sdm$sd <- as.data.frame(rbind(b0_sd,b1_sd,b2_sd,b3_sd,b4_sd,b5_sd,b6_sd))
colnames (coef_data_sdm) <- c("mean","sd")

write.csv(coef_data_sdm, file = "parameters_sdm_deb_aug_1.csv")


