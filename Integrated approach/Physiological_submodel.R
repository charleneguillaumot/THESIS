#------------------------------------------------------
# R scripts associated with the article
# Simple or hybrid? Next generation ecological models to study the distribution of Southern Ocean marine species
# by Guillaumot Charlène, Buba Yehezkel, Belmaker Jonathan, Fourcy Damien, Dubois Philippe, Danis Bruno, Saucède Thomas
# updated January 2021
#
# File: Physiological submodel. The output file is the one that is used in the Bayesian integrated approach

#--------------------------------------------------------------------------------------------------------
library(rgdal)
library(raster)
library(R.matlab)

## Environmental files
#---------------------
depth <- raster("envi/bathymetry_Morbihan_crop_finale.asc")
f_09022017 <- raster("envi/f_layer_09022017_OK2.asc"); f_09022017 <- mask(crop(f_09022017, depth), depth)
# for this layer, artefact due to clouds -> replace values = 0 to NA 
f_09022017[f_09022017==0]<-NA

f_20082017 <- raster("envi/f_layer_20082017_OK2.asc"); f_20082017 <- mask(crop(f_20082017, depth), depth)
temp_09022017 <- raster("envi/SST_MUR09022017_OK2.asc")
temp_20082017 <- raster("envi/SST_MUR20082017_OK2.asc")

f_layer <- f_20082017
temp_layer <- temp_20082017 # to be changed according to the studied period (August/February)

# DEB model
#--------------
source("scripts_DEB/get_DEB_pars.R")
source("scripts_DEB/get_powers_j.R")
source("scripts_DEB/get_tj2.R")
source("scripts_DEB/ode_functs.R")
source("scripts_DEB/t_corr.R")

f_i <- na.omit(values(f_layer))
# random partition of f values 
f_i <- f_i[sample(length(f_i),round(length(f_i)/4,0),replace = F)] 

repli<-25 # number of replicates in environmental conditions and sea urchin lenght
pG_empty <- matrix(nrow=length(f_i),ncol=1,data=NA)
pG_stack <- matrix(nrow=length(f_i),ncol=repli,data=NA)
head(pG_stack)
result_abatus <- get_DEB_pars("scripts_DEB/results_Abatus_cordatus.mat") # parameters of the DEB model of Abatus cordatus 
# see https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/Abatus_cordatus/Abatus_cordatus_res.html
compt <- 0

for (j in 1:repli){
  
  # random sampling of n=repli sea urchin lengths (contained between 2.5-4.5cm)
  L_rand <- runif(n=1, min=2.5,max=4.5)
  
  # random sampling of n=repli temperatures (contained within the T° range of the season)
  temp_rand <- runif(n=1, min=min(values(temp_layer), na.rm=T),max=max(values(temp_layer), na.rm=T))
  
  for (i in f_i) { # calculate growth performance for these environmental replicates 
    compt <- compt + 1
    s_M <- get_tj(p=result_abatus,f=i)$s_M
    
    powers_calc <- get_powers_j(p=result_abatus,L=L_rand,e=i, s_M)
    pG_empty[compt] <- powers_calc[5]  / result_abatus$k_M / t_corr(p=result_abatus, Ti=(273.15+temp_rand))
  }
  compt <- 0
  
  pG_stack[,j] <- pG_empty
}

head(pG_stack)

pG_mean <- apply(pG_stack, mean,MARGIN=1, na.rm=T)
pG_sd <- apply(pG_stack, sd,MARGIN=1, na.rm=T)

head(pG_mean)
head(pG_sd)

# 
food_vector <- f_i
physio_par <- pG_mean
se <- pG_sd
max_physio_par <- quantile(physio_par, 0.90)
min_physio_par <- quantile(physio_par, 0.10)
physio_par <-  (physio_par-min_physio_par)/(max_physio_par-min_physio_par) # scale
#physio_par[physio_par<0]=0
#physio_par[physio_par>1]=1
SE <- se/(max_physio_par-min_physio_par) # std scaled

plot(food_vector,physio_par, ylab="growth performance", xlab="food")

# ## poly_prior_new
# #----------------
library(betareg)
library(AICcmodavg)
library(nlme)
library(minpack.lm)

phydata <- data.frame(
  food_vector,
  food_vector_sq = food_vector^2,
  physio_par,
  SE,
  tau= 1/(SE)^2
)
# don't look higher than ^2
head(phydata)

phydata$psi <- 1-pnorm(0,phydata$physio_par,phydata$SE)#pnorm gives the distribution function
phydata$psi[phydata$psi>0.9999995]=0.9999995
phydata$psi[phydata$psi<0.0000005]=0.0000005
head(phydata)

plot(food_vector,phydata$psi, xlab="food")

# AIC- graph fitting:
ctrl <- betareg.control(maxit=1000)

#run several models and choose the one for which the AIC score is the lowest (i.e. best model)
beta_poly <- betareg(psi~food_vector+food_vector_sq, phydata)
aicc_poly <- AICc(beta_poly, second.ord = FALSE) # does not work because data not contrasted enough, to be tried with another DEB output?
#
beta_lin <- betareg(psi~food_vector, phydata)
aicc_lin <- AICc(beta_lin, second.ord = FALSE)

# Choose the best one...
#------------------------
aicc_all <- data.frame(aicc_poly,aicc_lin)
aicc_all <- sort(aicc_all,decreasing = FALSE)
head(aicc_all)

#
if (aicc_all[1]==aicc_lin){ # save the a,b,c coefficients of the selected polynomial model
  a <- coef(beta_lin)[1]
  b <- coef(beta_lin)[2]
} else if(aicc_all[1]==aicc_poly){
  a <- coef(beta_poly)[1]
  b <- coef(beta_poly)[2]
  c <- coef(beta_poly)[3]
}

source("scripts/ex1_globals.r")
print("create bPrior")
bPrior <- matrix(c(
  a, 1.0E-3,	# mean and tau for a0
  b, 1.0E-3,	#  a1
  c, 1.0E-3),	#a2
  byrow=TRUE, ncol=2)

# Calculate MCMC model from the previous polynom
print("Calculate MCMC model from the previous polynom")
phymodel_data <- list(psi=phydata$psi, food_vector=phydata$food_vector, bPrior=bPrior)
phymodel_par <- c('a0', 'a1', 'a2', 'phi')
model2 <- jags.model("scripts/phy_model2.jags", data = phymodel_data, n.chains=settings$chains, n.adapt=settings$tuning)
update(model2, settings$burnin) # burn in the model

phymodel_results <- coda.samples(model2, phymodel_par, settings$samples*settings$thin, thin=settings$thin)
head(phymodel_results)

#extracting coefficients from the model results
#and then calculate the mean and se:
mo_data <- as.matrix(phymodel_results)
head(mo_data)
a0 <- mo_data[,1]
a1 <- mo_data[,2]
a2 <- mo_data[,3]

prob <- data.frame(matrix(NA, nrow = length(a0), ncol = length(food_vector)))

library(foreach)
foreach (i= 1:length(food_vector), j = 1:length(a0)) %do% {
  prob[j,i] <- a0[j] + (a1[j]*food_vector[i]) + (a2[j] * food_vector[i]^2) # $
  prob[j,i] <- boot::inv.logit(prob[j,i])
}

# store the average values of the parameters
# these averages will be used as priors for the Bayesian integrated model
a0_m <- mean(a0)
a1_m <- mean(a1)
a2_m <- mean(a2)
a0_sd <- sd(a0)
a1_sd <- sd(a1)
a2_sd <- sd(a2)

coef_data <- as.data.frame(rbind(a0_m,a1_m,a2_m))
names(coef_data)[1] <- paste("mean")
coef_data$sd <- rbind(a0_sd,a1_sd,a2_sd)
write.csv(coef_data, file = "results_examples/physio_data_aug.csv")
