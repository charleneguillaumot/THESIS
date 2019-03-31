clock2 <- function(occ, bg.coords){
# the spatial sampling aims at defining areas containing either test or training data 
# for the CLOCK-2 method, a single diagonal transect is defined, cutting the Southern Ocean circle into 2 areas. One of the 2 areas contains the data that will be used to train the model, the other, the data that will be used to test the model after predictions 

#initialise vectors that are used to generate the spatial sampling structure
random_long <- seq(-180,180,1)
random_long_opp <- rep(NA,361)

for (i in 1:361){
  if (random_long[i] >= 0){
    random_long_opp[random_long > 0] <- random_long[random_long > 0]-180
  } else{
    random_long_opp[random_long <= 0] <- random_long[random_long <= 0]+180
  }
}

# sample a number between -180 and 180 to define the random sampling transect 
random_long_tirage <- sample(random_long,1)
random_long_tirage_t <- random_long_opp[random_long==random_long_tirage]

## define training and test groups (composed of both presence and background data)

occ.grp <- rep(NA, nrow(occ))
bg.coords.grp <- rep(NA, nrow(bg.coords))
  
if (random_long_tirage>0){
  
  for (i in 1:length(occ[,1])){
    if ((occ[i,1])<random_long_tirage &(occ[i,1])>random_long_tirage_t){
      occ.grp[i] <- 1 }
    else {
      occ.grp[i] <- 2
      }
    }
  
  for (i in 1:length(bg.coords[,1])){
    if ((bg.coords[i,1])<random_long_tirage &(bg.coords[i,1])>random_long_tirage_t){
      bg.coords.grp[i] <- 1 }
    else {
        bg.coords.grp[i] <- 2
    }
  }}

if (random_long_tirage<0) {
  for (i in 1:length(occ[,1])){
    if ((occ[i,1])>random_long_tirage &(occ[i,1])<random_long_tirage_t){
      occ.grp[i] <- 1 }
    else {
      occ.grp[i] <- 2
    }
  }
  for (i in 1:length(bg.coords[,1])){
    if ((bg.coords[i,1])>random_long_tirage &(bg.coords[i,1])<random_long_tirage_t){
      bg.coords.grp[i] <- 1 }
    else {
      bg.coords.grp[i] <- 2
    }
  }
}
  
out <- list(occ.grp=occ.grp,bg.coords.grp= bg.coords.grp)
return(out)
}
