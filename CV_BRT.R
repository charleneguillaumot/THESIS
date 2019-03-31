######################################################################################################
## Paper "Broad-scale species distribution models applied to data-poor areas", Progress in Oceanography
## 04/2019
## Guillaumot Charlène < charleneguillaumot21@gmail.com>
## SCRIPT for application of spatial cross-validation procedures
## use of Boosted Regression Trees
########################################################################################################
## Library and function
library(raster)
library(dismo)
library(gbm) ##
library(geosphere)
library(ape)

source("scripts/Function_ja.r")

## Parameters
cv.boot <- 100 # Number of cross-validation replicates (sample of background data + define spatial transects that split the dataset in training and test subsets)
prof_max_species <- -1500 # maximum depth until which the species is found. The model predictions will be restrained to shallower depths
# BRT calibration parameters
tc = 4   # tree complexity
lr = 0.007 # learning rate
bf = 0.75  # bag fraction

P_method <- "C2" # choose the cross-validation procedure to apply 
# R -> random
# B -> block.method
# C2 | C3 | C4 | C6 -> clock

## Assess the parameter of the K-fold cross-validation procedure
if(any(P_method == c("R", "B"))){ 
  P_K <- 4
} else {
  P_K <- as.numeric(substr(P_method, 2, 2))
} 

# ODONTASTER
#----------------
data <- read.csv("occurences/O_validus_occ_modif.csv", sep=";",dec=".", header=T)
odontaster.occ <- data[,c(1,2)]

#======================================================================================#
# Load environmental data
#-----------------------------------------
predictors_2005_2012 <- raster::stack("envi/predictors_2005_2012_ANT_final_withcurrent.grd")
# can be found in Fabri-Ruiz et al. 2017
predictors_2005_2012 <- dropLayer(predictors_2005_2012,16) 

# Limit the studied and projection area: crop the environmental maps to maximum depth
#--------------------------------------------------------
source("scripts/delim.area.R") 
extent_map <- extent(predictors_2005_2012)
predictors_2005_2012_crop_species <- delim.area (predictors = predictors_2005_2012,longmin = extent_map[1],longmax = extent_map[2], latmin =extent_map[3],latmax = extent_map[4],interval = c(0,prof_max_species))
predictors <- predictors_2005_2012_crop_species

#======================================================================================#
# KDE lAYER 
#-----------------------------------------
KDE <- raster("envi/bias.grd")
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

# Model preparation
#----------------------------------
stack.pred<-subset(predictors,1);values(stack.pred)<-NA
testvaluesp<-rep(NA,nrow(unique(odontaster.occ))) ; testvaluesa<-testvaluesp

eval_self <- matrix(nrow=18,ncol=cv.boot,data=NA)
colnames(eval_self)<-seq(1,ncol(eval_self),1)
rownames(eval_self)<-c("AUC", "maxSSS","COR","pCOR","TSS","test_gp", "valid_test_pres","nb_NA_test","nb_test","moran_resi","pval_resi","sd_moran_resi","moran_extract","pval_extract","sd_moran_extract","borne_sup","borne_inf","ntrees") 

eval_cv <- matrix(NA, 7, cv.boot*P_K, 
                  dimnames = list(c("AUC", "Cor", "maxSSS", "valid_test_data", "prop_cv", "TSS","valid_test_pres"), NULL))

# Stores of the contributions
contTr <- matrix(NA, dim(predictors)[3], cv.boot, dimnames = list(names(predictors), NULL))

# Convert the raster file of the environmental descriptors into data.frame to earn computing time
predDF <- as.data.frame(predictors) 
predDF <- na.omit(predDF) # delete NA pixel
numNA <- attr(predDF, "na.action") 
predR <- matrix(NA, nrow(predDF), cv.boot) 


## Get environmental descriptors values recorded on pixels containing presence or background data 
#-----------------------------------------------------------
presvals_glob <- extract(predictors, unique(odontaster.occ)) # extract value for presence data

# Random sampling of background data
background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
colnames(background_data) <- colnames(odontaster.occ)
pseudoabs_glob <- extract(predictors,background_data)

# Extract environmental values for presence data 
odontaster.occ <- cbind(odontaster.occ, IsP = 1, extract(predictors, odontaster.occ))

stack.pred<-subset(predictors,1);values(stack.pred)<-NA

# BRT Launching and replicates loop
#-----------------------------------------------------------------------------
for (j in 1:cv.boot){
  
  # Random sampling of background data for each replicate 
  background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 1000, prob=values(KDE)[!is.na(values(KDE))]))
  colnames(background_data) <- colnames(odontaster.occ)[1:2]
  background_data <- cbind(background_data, IsP = 0, extract(predictors, background_data))
  
  # Gather information of presence and background data 
  dat1 <- rbind(background_data, odontaster.occ)
  
  # Split the dataset into "folds", except for the random cross-validation procedure
  idP <- which(dat1$IsP == 1) # id of the presences
  MyFold <- rep(NA, nrow(dat1)) # empty box to store the folds
  
  if(P_method == "R"){ # If "Random"
    MyFold <- NULL}
    
    if (P_method == "B"){ # If - "B", "C2", "C3","C4","C6"
    # Runs the 'get.block' function from Function_ja.R
    blockF <- get.block(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
      
    # Extracts the folds
    MyFold[idP] <- blockF$occ.grp
    MyFold[-idP] <- blockF$bg.grp
      
    # Plots to check
    if(j%%10 == 0){
      MyPngName <- paste0("output/", P_method, "/Folds_", P_method, "_bo", j, ".png")
      png(MyPngName)
        plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue", "green", "gold2")[as.factor(MyFold)]  )
      dev.off()
      #points(ant_map, type="l")
    }}
  
    if(P_method == "C2"){
    
    source("clock2_cg.R")
    clock2F <- clock2(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
    
    # Extracts the folds
    MyFold[idP] <- clock2F$occ.grp
    MyFold[-idP] <- clock2F$bg.coords.grp
    plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue")[as.factor(MyFold)]  )
    
    # Plots to check
    if(j%%10 == 0){
      MyPngName <- paste0("output/", P_method, "/Folds_", P_method, "_bo", j, ".png")
      png(MyPngName)
      plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue")[as.factor(MyFold)]  )
      dev.off()}}
      
    if(P_method == "C3"){
    source("clock3-cg.R")
    clock3F <- clock3(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
    
    # Extracts the folds
    MyFold[idP] <- clock3F$occ.grp
    MyFold[-idP] <- clock3F$bg.coords.grp
    plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black")[as.factor(MyFold)]  )
    
    # Plots to check
    if(j%%10 == 0){
      MyPngName <- paste0("output/", P_method, "/Folds_", P_method, "_bo", j, ".png")
      png(MyPngName)
      plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black")[as.factor(MyFold)]  )
      dev.off()}}
      
    if(P_method == "C4"){
    source("clock4_cg.R")
    clock4F <- clock4(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
    
    # Extracts the folds
    MyFold[idP] <- clock4F$occ.grp
    MyFold[-idP] <- clock4F$bg.coords.grp
    plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black","purple")[as.factor(MyFold)])
  
    # Plots to check
    if(j%%10 == 0){
      MyPngName <- paste0("output/", P_method, "/Folds_", P_method, "_bo", j, ".png")
      png(MyPngName)
      plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black","purple")[as.factor(MyFold)])
      dev.off()}}
 
  if(P_method == "C6"){
    source("clock6_cg.R")
    clock6F <- clock6(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
    
    # Extracts the folds
    MyFold[idP] <- clock6F$occ.grp
    MyFold[-idP] <- clock6F$bg.coords.grp
    plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black","purple","orange","green")[as.factor(MyFold)])
    # on plot certain cv pour verifier
    if(j%%10 == 0){
      MyPngName <- paste0("output/", P_method, "/Folds_", P_method, "_bo", j, ".png")
      png(MyPngName)
      plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black","purple","orange","green")[as.factor(MyFold)])
      dev.off()
    }}
  
  #----------------------
  # LAUNCH BRT ## /!\ gbm.step_v2
  #----------------------
  model.res <- gbm.step_v2(data = dat1, 
                           gbm.x = 4:ncol(dat1), # all environmental descriptors
                           gbm.y = 3, # lat, long et Isp
                           n.folds = P_K, #### /!\ for random CV
                           fold.vector = MyFold, #### for other spatial CV with folds
                           family = "bernoulli",
                           tree.complexity = tc,
                           learning.rate = lr,
                           bag.fraction = bf)
  
  
  #------------------------------------------------------------
  # Model outputs  
  #------------------------------------------------------------
  p<-predict(predictors,model.res,n.trees=model.res$gbm.call$best.trees,type="response", na.rm=F)
  p <- mask(p,subset(predictors,1))
  stack.pred<-stack(stack.pred,p) # stack all the maps replicates 
  
  # Get number of trees
  eval_self["ntrees", j] <- model.res$n.trees
  
  ##########
  ## SELF ## (self = on all data)
  ##########
  
  # AUC score
  eval_self["AUC", j] <- model.res$self.statistics$discrimination
  
  # Point biserial correlation 
  eval_self["COR", j] <- model.res$self.statistics$correlation


  ########
  ## CV ##
  ########
    j_cv <- ((j-1) * P_K+1):(j*P_K) # count to fill the result tables
  
  eval_cv["AUC", j_cv] <- model.res$cv.roc.matrix
  eval_cv["Cor", j_cv] <- model.res$cv.cor.matrix
  eval_cv["maxSSS", j_cv] <- model.res$cv.th.matrix
  eval_cv["valid_test_data", j_cv] <- model.res$cv.corr.class*100
  eval_cv["prop_cv", j_cv] <- model.res$cv.length*100
  eval_cv["TSS", j_cv] <- model.res$tss.cv
  eval_cv["valid_test_pres", j_cv] <- model.res$cv.length*100
  
  RI <- summary(model.res, plotit = F) # extract the contribution
  contTr[match(RI$var, rownames(contTr)), j] <- RI[,"rel.inf"]
  
}


### I. Predictions
# Average prediction and convert to raster format
predRtot <- rep(NA, length(predictors[[1]]))
predRtot[-numNA] <- apply(predR, 1, mean, na.rm=T)
predRtot <- setValues(predictors[[1]], predRtot)
writeRaster(predRtot, paste0("output/", P_method, "/pred_mean_", P_method, ".asc") )

mean_stack <- raster::calc(stack.pred, mean, na.rm=T)
sd_stack <- raster::calc(stack.pred,sd, na.rm=T)

writeRaster(mean_stack, "output/C2/mean_raster_predNA.asc")
### II. SELF
esM <- apply(eval_self, 1, mean, na.rm=T)
esSD <- apply(eval_self, 1, sd, na.rm=T)
esTot <- paste(round(esM, 3), round(esSD, 3), sep = " ± ")
names(esTot) <- names(esM)

### III. CV
#eval_cv["valid_test_pres",] <- eval_cv["valid_test_pres",]*100
ecM <- apply(eval_cv, 1, mean, na.rm=T)
ecSD <- apply(eval_cv, 1, sd, na.rm=T)
ecTot <- paste(round(ecM, 3), round(ecSD, 3), sep = " ± ")
names(ecTot) <- names(ecM)

# % of data included in the "test data"
r1 <- round(range(eval_cv["prop_cv"])*100)
r1 <- paste0(paste(r1, collapse = "-"), "%")

## Export results
Ncol <- c(R = "Standard CV Random splitting", 
          B = "Spatial CV Block method",
          C2 = "2-fold Clock method",
          C3 = "3-fold Clock method",
          C4 = "4-fold Clock method",
          C6 = "6-fold Clock method")

ResF <- data.frame(c(ecTot["AUC"], ecTot["valid_test_data"], r1, ecTot["Cor"], esTot["ntrees"], ecTot["TSS"], ecTot["valid_test_pres"]))
rownames(ResF) <- c("AUC", "Correctly classified test data (% of total data)", "Test data (% of total dataset)", "COR", "ntrees", "TSS", "Correctly classified test data (% of presence data)")
colnames(ResF) <- Ncol[P_method]

write.csv(ResF, paste0("output/", P_method, "/Table2_", P_method, ".csv"))

ResF

write.csv(eval_cv, paste0("output/", P_method, "/eval_cv_", P_method, ".csv"))


### IV. Contribution
CtM <- apply(contTr, 1, mean)
CtSD <- apply(contTr, 1, sd)
CtTot <- paste(round(CtM, 3), round(CtSD, 3), sep = " ± ")
names(CtTot) <- names(CtM)

CtTot <- data.frame(CtTot)
colnames(CtTot) <- Ncol[P_method]
write.csv(CtTot, paste0("output/", P_method, "/Contribution_", P_method, ".csv"))


### contri MESS
stack1_mean <- raster("output/C2/mean_raster_predNA.asc")
MESS <- raster("output/mess_O_validus_1500.asc")
maxSSS_mean <- 0.177673

stack1_mean_favo <- stack1_mean > maxSSS_mean
values(stack1_mean_favo)
stack1_mean_favo <- reclassify(stack1_mean_favo,cbind(FALSE,0))
stack1_mean_favo <- reclassify(stack1_mean_favo,cbind(TRUE,1))
plot(stack1_mean_favo)

# Produce raster layers to study MESS/mean SDM overlap
favo_MESSOK_raster <- stack1_mean_favo==1 & MESS==1
favo_MESSoutside <- stack1_mean_favo==1 & MESS==0
defavo_MESSOK <- stack1_mean_favo==0 & MESS==1
defavo_MESSoutside <- stack1_mean_favo==0 & MESS==0

# Percentages 
favo_MESSOK_prop <- 100*length(which(values(favo_MESSOK_raster)==1))/ length(na.omit(values(favo_MESSOK_raster)))
favo_MESSoutside_prop <- 100*length(which(values(favo_MESSoutside)==1))/ length(na.omit(values(favo_MESSoutside)))
defavo_MESSOK_prop <- 100*length(which(values(defavo_MESSOK)==1))/ length(na.omit(values(defavo_MESSOK)))
defavo_MESSoutside_prop <- 100*length(which(values(defavo_MESSoutside)==1))/ length(na.omit(values(defavo_MESSoutside)))

length(which(values(MESS==0)))/length(subset(predictors,1))
MESS2 <- reclassify(MESS,cbind(1,NA))
length(which(!is.na(values(MESS2))))*100 /length(which(!is.na(values(subset(predictors,1       )))))
