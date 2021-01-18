#############################################################################################################
## Paper "Influence of environmental descriptor choice in modelling Southern Ocean benthic species distributions", Diversity and Distribution
## 12/2018
## Guillaumot Charlène < charleneguillaumot21@gmail.com>
## SCRIPT to create extreme event layers (minimal extreme event)
## Downloadable at https://github.com/AustralianAntarcticDivision/blueant
## Tutorial and metadata https://australianantarcticdivision.github.io/blueant/articles/SO_SDM_data.html
#############################################################################################################

# The function uses a RasterStack for which each layer contains the environmental values of each month 
# Per pixel, it estimates the number of months that contain an extreme event 
# For each pixel, the median and its variability (Median Absolute Deviation, MAD) are assessed. The code then evaluates whether the value is lower than the median and its variability. This is evaluated for each month. The number of months being a minimal extreme event per pixel is finally counted and saved as a raster layer. 

extremeMonths_min <- function(rasterstack){
  
  rastervide <- subset(rasterstack,1); values(rastervide) <- NA
  
  for (i in 1:ncell(rasterstack)){
    compt <- 0
    VEC <- as.vector(values(rasterstack)[i,])
    
    if (length(which(is.na(VEC)))>=length(VEC)) {
      compt <- NA} 
    else {
      VECmin <- VEC - mad(VEC, na.rm=T)
      VECmin <- replace(VECmin, is.na(VECmin),0)
      VECmax <- VEC + mad(VEC, na.rm=T)
      VECmax <- replace(VECmax, is.na(VECmax),0)
      MED <- median(VEC, na.rm=T)
      
      for (j in 1:nlayers(rasterstack)){
        if(VECmin[j]<MED & VECmax[j]<MED){
          compt <- compt+1}
      }
    }
    rastervide[i] <- compt
  }
  return(rastervide)
}