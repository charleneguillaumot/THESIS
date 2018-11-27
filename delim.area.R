delim.area <- function(predictors, longmin, longmax, latmin, latmax, interval=NULL, crslayer = raster::crs(predictors)) {
  
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("The function requires 'raster' package to work, please install it.", call. = FALSE)
  }
  
  data <- predictors
  data <- raster::crop(data, raster::extent(longmin, longmax, latmin, latmax))
  raster::crs(data) <- crslayer
  stackpred.final <- data
  
  if(is.null(interval)==FALSE){
    # transform the NA values in -99999 values
    layertochange <- raster::subset(data, subset = 1)
    layertochange <- raster::reclassify(layertochange, cbind(NA, -99999))
    layertochange_v <- raster::rasterToPoints(layertochange)
    
    nullraster <- layertochange
    raster::values(nullraster) <- 0  # create a null raster to implement rasterize
    
    if(interval[1]>interval[2]){
      for (i in 1:nrow(layertochange_v)) {
        if (layertochange_v[i, 3] <= interval[1] & layertochange_v[i, 3] >= interval[2]) {
          layertochange_v[i, 3] <- layertochange_v[i, 3]
        } else {
          
          layertochange_v[i, 3] <- (-99999)
        }
      }
    }
    
    if(interval[1]<interval[2]){
      for (i in 1:nrow(layertochange_v)) {
        if (layertochange_v[i, 3] >= interval[1] & layertochange_v[i, 3] <= interval[2]) {
          layertochange_v[i, 3] <- layertochange_v[i, 3]
        } else {
          
          layertochange_v[i, 3] <- (-99999)
        }
      }
    }
    
    
    layertochange <- raster::rasterize(as.matrix(layertochange_v[, 1:2]), nullraster, as.vector(layertochange_v[,
                                                                                                                3]))
    layertochange <- raster::reclassify(layertochange, base::cbind(-99999, NA))
    layertochange <- raster::crop(layertochange, raster:: extent(longmin, longmax, latmin, latmax))
    
    stackpred.final <- raster::stack(layertochange, raster::subset(data, c(2:raster::nlayers(data))))
    names(stackpred.final) <- c(names(predictors))
  }
  
  return(stackpred.final)
  
}