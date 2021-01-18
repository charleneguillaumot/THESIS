clock3 <- function(occ, bg.coords){
  # the spatial sampling aims at defining areas containing either test or training data 
  # for the CLOCK-3 method, 2 diagonal transects are defined, cutting the Southern Ocean circle into 3 areas. Two of the 3 areas contains the data that will be used to train the model, the remaining one, the data that will be used to test the model after predictions 
  
  #initialise vectors that are used to generate the spatial sampling structure
  random_longA <- seq(0,179,1)
  random_longB <- seq(-180,-1,1)
  random_longC <- seq(0,180,1)
  random_longD <- seq(-179,0,1)
  random_long <- c(random_longA,random_longB,random_longC,random_longD)
  
  # sample a number between -180 and 180 to define the random sampling transect 
  tirage <- sample(seq(181,541,1),1)
  random_long_tirage <- random_long[tirage]
  random_long_tirage_right <- random_long[tirage+120]
  random_long_tirage_left <- random_long[tirage-120]
  
## define training and test groups (composed of both presence and background data)
occ.grp <- rep(NA, nrow(occ))
bg.coords.grp <- rep(NA, nrow(bg.coords))
  

## define training and test groups (composed of both presence and background data)
presence_tot <- occ
training_presences.occ <- NA;training_presences.occ_left <-NA;training_presences.occ_right <-NA
training_backgr.occ <- NA; training_backgr.occ_left <- NA;training_backgr.occ_right <-NA
test_presences.occ <- NA; training_presences.occ_inf_left <- NA; training_presences.occ_inf_right <- NA; training_presences.occ_supp_left <- NA; training_presences.occ_supp_right <- NA;training_backgr.occ_supp_left <- NA; training_backgr.occ_supp_right <- NA; training_backgr.occ_inf_left <- NA; training_backgr.occ_inf_right <- NA;
background_data <- bg.coords

# training presences
#---------------------
if (random_long_tirage_left < 0){
  training_presences.occ_left <- which(presence_tot[,1] > random_long_tirage_left & presence_tot[,1] < 0 )
} else {
  training_presences.occ_left <- which(presence_tot[,1] > random_long_tirage_left)
}

if(random_long_tirage_right>0){
  training_presences.occ_right <- which(presence_tot[,1] <random_long_tirage_right & presence_tot[,1] > 0)
} else {
  training_presences.occ_right <- which(presence_tot[,1] <random_long_tirage_right)
}

if(random_long_tirage_right>0 & random_long_tirage_left>0){
  training_presences.occ_supp_left <- which(presence_tot[,1] <random_long_tirage & presence_tot[,1] > -179)
  training_presences.occ_supp_right <- which(presence_tot[,1] >random_long_tirage & presence_tot[,1] < 0)
} 

if(random_long_tirage_right<0 & random_long_tirage_left<0){
  training_presences.occ_inf_left <- which(presence_tot[,1] <random_long_tirage & presence_tot[,1] > 0)
  training_presences.occ_inf_right <- which(presence_tot[,1] >random_long_tirage & presence_tot[,1] < 180)
} 

training_presence_left_all <- c(training_presences.occ_left,training_presences.occ_inf_left,training_presences.occ_supp_left)
training_presence_right_all <- c(training_presences.occ_right,training_presences.occ_inf_right,training_presences.occ_supp_right)

for (i in training_presence_left_all){
  occ.grp[i] <- 1
}
for (i in training_presence_right_all){
  occ.grp[i] <- 2
}

occ.grp[which(is.na(occ.grp))]=3

# Background
#---------------------
if (random_long_tirage_left < 0){
  training_backgr.occ_left <- which(background_data[,1] > random_long_tirage_left & background_data[,1] < 0 )
} else {
  training_backgr.occ_left <- which(background_data[,1] > random_long_tirage_left)
}

if(random_long_tirage_right>0){
  training_backgr.occ_right <- which(background_data[,1] <random_long_tirage_right & background_data[,1] > 0)
} else {
  training_backgr.occ_right <- which(background_data[,1] <random_long_tirage_right)
}

if(random_long_tirage_right>0 & random_long_tirage_left>0){
  training_backgr.occ_supp_left <- which(background_data[,1] <random_long_tirage & background_data[,1] > -179)
  training_backgr.occ_supp_right <- which(background_data[,1] >random_long_tirage & background_data[,1] < 0)
} 

if(random_long_tirage_right<0 & random_long_tirage_left<0){
  training_backgr.occ_inf_left <- which(background_data[,1] <random_long_tirage & background_data[,1] > 0)
  training_backgr.occ_inf_right <- which(background_data[,1] >random_long_tirage & background_data[,1] < 180)
} 

training_backgr_left_all <- c(training_backgr.occ_left,training_backgr.occ_inf_left,training_backgr.occ_supp_left)
training_backgr_right_all <- c(training_backgr.occ_right,training_backgr.occ_inf_right,training_backgr.occ_supp_right)

for (i in training_backgr_left_all){
  bg.coords.grp[i] <- 1
}
for (i in training_backgr_right_all){
  bg.coords.grp[i] <- 2
}

bg.coords.grp[which(is.na(bg.coords.grp))]=3

out <- list(occ.grp=occ.grp,bg.coords.grp= bg.coords.grp)
return(out)
}
