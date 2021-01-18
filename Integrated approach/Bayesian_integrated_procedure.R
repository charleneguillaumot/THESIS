#------------------------------------------------------
# R scripts associated with the article
# Simple or hybrid? Next generation ecological models to study the distribution of Southern Ocean marine species
# by Guillaumot Charlène, Buba Yehezkel, Belmaker Jonathan, Fourcy Damien, Dubois Philippe, Danis Bruno, Saucède Thomas
# updated January 2021
#
# File: procedure for the 'Bayesian integrated' approach
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

# use the physiological submodel and the simple SDM outputs to compute priors for the Bayesian integrated model 
# open results: 
#---------------
coef_data_physio <- read.csv(file = "results_examples/physio_data_aug.csv", sep=";") # from the Physiological submodel file
coef_data_physio$tau <- 1/(coef_data_physio$sd^2)

coef_data_simple_sdm <- read.csv(file = "results_examples/simple_sdm_data_aug.csv", sep=";") # from the simple SDM procedure file
coef_data_simple_sdm$tau <- 1/(coef_data_simple_sdm$sd^2)

# OPEN OCCURRENCE RECORDS
#--------------------------
abatus.occ <- read.csv("data/abatus.occ_modif.csv", header=T, sep=";")[,c(1,2)]
head(abatus.occ)

sp.occr_loc=data.frame(lon=abatus.occ$decimalLongitude,lat=abatus.occ$decimalLatitude)
sp_name=data.frame(Occurrence=rep(1,nrow(sp.occr_loc)))
species_occ=SpatialPointsDataFrame(coords=sp.occr_loc,data=sp_name)

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

f_layer <- f_20082017
temp_layer <- temp_20082017  # to be changed according to the studied period (August/February)

predictors <- stack(depth, temp_layer,f_layer)
names(predictors) <- c("depth","temp","f")


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

#---------------------------------------------------------------------------
#coef_data_physio -> coefficients of the physiological model
#coef_data_simple_sdm -> coefficients of the sdm
## EQUATION PHYSIO SUB-MODEL 
## a0 + a1* food_vector + a2* food_vector^2

## EQUATION SDM
## b0 + b1* depth + b2*food_vector + b3*temp + b4* temp^2 + b5* food_vector^2

##EQUATION COMBINED SDM
## a0 - b1 - a1 - b3 - b4 - a2
## c0 - c1 - c2 - c3 - c4- c5
# for the vague priors (coming from the simple SDM, we decrease the tau value (high uncertainty), to have a SD of around 10, so tau= 1/10^2= 0.01)
#---------------------------------------------------------------------------

bPrior <- matrix(c(
  coef_data_simple_sdm$mean[1], 0.01,	# c0 => a0 (intercept)-> physio model
  coef_data_simple_sdm$mean[2], 0.01,	# c1 => b1 (depth) -> simple sdm, vague prior
  coef_data_physio$mean[2], coef_data_physio$tau[2],	# c2 => a1 (food_vector)-> physio model
  coef_data_simple_sdm$mean[4], 0.01,	# c3 => b3 (temp) -> simple sdm, vague prior
  coef_data_simple_sdm$mean[5], 0.01,	# c4 => b4 (temp^2) -> simple sdm, vague prior
  coef_data_physio$mean[3], coef_data_physio$tau[3]),	# c5 => a2 (food_vector^2)-> physio model
  byrow=TRUE, ncol=2)

m1Data <- list(depth = pred$depth, temp = pred$temp, f = pred$f, presence = pred$presence, bPrior=bPrior)
m1Pars <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5')
#jags.inits <- list(b0=coef_data_physio$mean[1], b1=coef_data_simple_sdm$mean[2],b2=coef_data_physio$mean[2],b3=coef_data_simple_sdm$mean[4], b4=coef_data_simple_sdm$mean[5], b5=coef_data_physio$mean[3])
jags.inits <- list(b1=coef_data_simple_sdm$mean[2],b2=coef_data_physio$mean[2],b3=coef_data_simple_sdm$mean[4], b4=coef_data_simple_sdm$mean[5], b5=coef_data_physio$mean[3])

# Empty matrix to fill 
coef_b <- matrix(nrow = 4000, ncol = length(m1Pars)); colnames(coef_b) <- m1Pars
colnames(coef_b) <- m1Pars
data_list <- vector(mode = "list")

# create the model
#jags.inits[] = 0
#,inits=jags.inits
model1 <- jags.model("scripts/sdm_model_poly.jags", data = m1Data, n.chains=settings$chains, n.adapt=settings$tuning, inits=jags.inits  )
update(model1, settings$burnin)

m1Results <- coda.samples(model1, m1Pars, settings$samples*settings$thin, thin=settings$thin)

# PREDICTIONS 
##-------------
# fonctions 
formula <- function(x,b){
  boot::inv.logit(b[1] + b[2]*x[1] + b[3]*x[2]+ b[4]*x[3] + b[5]*x[2]^2 + b[6]*x[3]^2)}

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
predictors_v <- as.data.frame(na.omit(rasterToPoints(predictors)))
predictors_data <- predictors_v[,-c(1,2)]
head(predictors_data)

prob_sp <- predictions(predictors_data,coef_b)

prob_sp_mean <- colMeans(prob_sp)
max_prob_sp_mean <- max(prob_sp_mean)
SD_prob_sp <- apply(prob_sp, 2, sd)

p <- rasterFromXYZ(data.frame(x=predictors_v[,1],y=predictors_v[,2], z=prob_sp_mean))
#plot(p)
#points(abatus.occ, pch=20)
writeRaster(p, "result_combined_sdm_mean_august.asc")
p_sd <- rasterFromXYZ(data.frame(x=predictors_v[,1],y=predictors_v[,2], z=SD_prob_sp))
writeRaster(p_sd, "result_combined_sdm_sd_august.asc")
write.csv(background_data,"background_data_aug.csv")
