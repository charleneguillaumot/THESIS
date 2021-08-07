library(raster)
library(SpaDES.tools)
library(dismo)

continent <- read.csv("results_maps/worldmap-0.01.csv")

#shape <- raster::shapefile("Proposed MPA/mpa_AP_proposal.shp")
#crs(shape)<-crs(x)

# partie utilisée avant quand testait les stats sur des gros pixels !#
#.......................................................................
# rasterzonal <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)
# rasterzonal2 <- SpaDES.tools::splitRaster(rasterzonal,nx=4,ny=24)
# for (i in 1:96){
#   values(rasterzonal2[[i]]) <- i
# }
# rasterzonal <- mergeRaster(rasterzonal2)
# plot(rasterzonal)
# points(continent,type="l")

#---------------------
# NUMBER OF PARTICLES 
#---------------------
BW5_10m_R03_OND_2008_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2008_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2009_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2009_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2010_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2010_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2011_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2011_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2012_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2012_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2013_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2013_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2014_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2014_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2015_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2015_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]

BW5_10m_R03_OND_2016_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW5-001_release_all_10-11-12_year_2016_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]



#### Transform the csv into rasters
#------------------------------------
rastervide <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)

BW5_10m_R03_OND_2008 <-rasterize(BW5_10m_R03_OND_2008_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2008_csv[,3])
BW5_10m_R03_OND_2009 <-rasterize(BW5_10m_R03_OND_2009_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2009_csv[,3])
BW5_10m_R03_OND_2010 <-rasterize(BW5_10m_R03_OND_2010_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2010_csv[,3])
BW5_10m_R03_OND_2011 <-rasterize(BW5_10m_R03_OND_2011_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2011_csv[,3])
BW5_10m_R03_OND_2012 <-rasterize(BW5_10m_R03_OND_2012_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2012_csv[,3])
BW5_10m_R03_OND_2013 <-rasterize(BW5_10m_R03_OND_2013_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2013_csv[,3])
BW5_10m_R03_OND_2014 <-rasterize(BW5_10m_R03_OND_2014_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2014_csv[,3])
BW5_10m_R03_OND_2015 <-rasterize(BW5_10m_R03_OND_2015_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2015_csv[,3])
BW5_10m_R03_OND_2016 <-rasterize(BW5_10m_R03_OND_2016_csv[,c(1,2)],rastervide,BW5_10m_R03_OND_2016_csv[,3])
plot(BW5_10m_R03_OND_2016)
points(continent,type="l")

### MAKE MAPS 
#--------------------------------
writeRaster(BW5_10m_R03_OND_2008,"maps_for_figures/BW5_10m_R03_OND_2008.tiff")
writeRaster(BW5_10m_R03_OND_2009,"maps_for_figures/BW5_10m_R03_OND_2009.tiff")
writeRaster(BW5_10m_R03_OND_2010,"maps_for_figures/BW5_10m_R03_OND_2010.tiff")
writeRaster(BW5_10m_R03_OND_2011,"maps_for_figures/BW5_10m_R03_OND_2011.tiff")
writeRaster(BW5_10m_R03_OND_2012,"maps_for_figures/BW5_10m_R03_OND_2012.tiff")
writeRaster(BW5_10m_R03_OND_2013,"maps_for_figures/BW5_10m_R03_OND_2013.tiff")
writeRaster(BW5_10m_R03_OND_2014,"maps_for_figures/BW5_10m_R03_OND_2014.tiff")
writeRaster(BW5_10m_R03_OND_2015,"maps_for_figures/BW5_10m_R03_OND_2015.tiff")
writeRaster(BW5_10m_R03_OND_2016,"maps_for_figures/BW5_10m_R03_OND_2016.tiff")


#CALCULATE AVERAGE YEARS
#------------------------

BW5_10m_R03_OND_2008 <-raster("maps_for_figures/BW5_10m_R03_OND_2008.tif")
BW5_10m_R03_OND_2009 <-raster("maps_for_figures/BW5_10m_R03_OND_2009.tif")
BW5_10m_R03_OND_2010 <-raster("maps_for_figures/BW5_10m_R03_OND_2010.tif")
BW5_10m_R03_OND_2011 <-raster("maps_for_figures/BW5_10m_R03_OND_2011.tif")
BW5_10m_R03_OND_2012 <-raster("maps_for_figures/BW5_10m_R03_OND_2012.tif")
BW5_10m_R03_OND_2013 <-raster("maps_for_figures/BW5_10m_R03_OND_2013.tif")
BW5_10m_R03_OND_2014 <-raster("maps_for_figures/BW5_10m_R03_OND_2014.tif")
BW5_10m_R03_OND_2015 <-raster("maps_for_figures/BW5_10m_R03_OND_2015.tif")
BW5_10m_R03_OND_2016 <-raster("maps_for_figures/BW5_10m_R03_OND_2016.tif")

plot(BW5_10m_R03_OND_2008)
#points(continent, type="l")

BW5_10m_R03_OND_2008_2016_stack <- stack(BW5_10m_R03_OND_2008,BW5_10m_R03_OND_2009,BW5_10m_R03_OND_2010,BW5_10m_R03_OND_2011,BW5_10m_R03_OND_2012,BW5_10m_R03_OND_2013,BW5_10m_R03_OND_2014,BW5_10m_R03_OND_2015,BW5_10m_R03_OND_2016)

BW5_10m_R03_OND_2008_2016_mean<- raster::calc(BW5_10m_R03_OND_2008_2016_stack, mean,na.rm=T)

plot(BW5_10m_R03_OND_2008_2016_mean)
#points(continent, type="l")

writeRaster(BW5_10m_R03_OND_2008_2016_mean,"maps_for_figures/per_zones/R03_RZ5/BW5_10m_R03_OND_2008_2016_meannarmT.tiff")

#BW5_10m_R03_OND_2008_2016_mean_relative <- BW5_10m_R03_OND_2008_2016_mean/445

#writeRaster(BW5_10m_R03_OND_2008_2016_mean_relative,"maps_for_figures/per_zones/R03/BW5_10m_R03_OND_2008_2016_mean_relative_narmF.tiff")


## MOYENNE GENERALE ANNUELLE
BW5_10m_R0X_OND_2008_2016_mean_relative <- raster("maps_for_figures/per_zones/R0X/BW5_10m_R0X_OND_2008_2016_mean_relative_narmF.tif")
BW5_10m_R0X_OND_2008_2016_mean_relative <- raster("maps_for_figures/per_zones/R0X/BW5_10m_R0X_OND_2008_2016_mean_relative_narmF.tif")
BW5_10m_R0X_OND_2008_2016_mean_relative <- raster("maps_for_figures/per_zones/R0X/BW5_10m_R0X_OND_2008_2016_mean_relative_narmF.tif")
BW5_10m_R0X_OND_2008_2016_mean_relative <- raster("maps_for_figures/per_zones/R0X/BW5_10m_R0X_OND_2008_2016_mean_relative_narmF.tif")

stack_BW5 <- stack(BW5_10m_R0X_OND_2008_2016_mean_relative,BW5_10m_R0X_OND_2008_2016_mean_relative,BW5_10m_R0X_OND_2008_2016_mean_relative,BW5_10m_R0X_OND_2008_2016_mean_relative)

#BW5_10m_R0X_all_seasons_2008_2016_mean<- raster::calc(stack_BW5, mean, na.rm=T)
BW5_10m_R0X_all_seasons_2008_2016_mean<- raster::calc(stack_BW5, mean)

writeRaster(BW5_10m_R0X_all_seasons_2008_2016_mean,"maps_for_figures/per_zones/R0X/BW5_10m_R0X_all_seasons_2008_2016_mean_narmF.tiff")



## PLOT FREQUENCE
#-----------------
stack_BW5_10m_R0X_all_season_2008 <- stack(BW5_10m_R0X_OND_2008,BW5_10m_R0X_OND_2008,BW5_10m_R0X_OND_2008,BW5_10m_R0X_OND_2008)
BW5_10m_R0X_all_season_2008 <- calc(stack_BW5_10m_R0X_all_season_2008,mean, na.rm=T)
plot(BW5_10m_R0X_all_season_2008)
#BW5_10m_R0X_all_season_2008_re <- reclassify(BW5_10m_R0X_all_season_2008, cbind(NA, -99))
BW5_10m_R0X_all_season_2008[values(BW5_10m_R0X_all_season_2008)>=0] <- 1
plot(BW5_10m_R0X_all_season_2008)
  








#### Schoener simple comparison
#------------------------------
nicheOverlap(BW3_10m_R0X_winter,BW3_200m_R0X_winter, stat="D")

# Spearman correlation 
cor(na.omit(cbind(values(BW3_10m_R0X_winter),values(BW3_200m_R0X_winter))), method = "spearman")

### BIG GRIDS ET WILCOXON 
#--------------------------------
layerToCompare1 <-BW3_10m_R0X_winter
layerToCompare2 <-BW3_200m_R0X_winter
layer <- "ALLareas"
raster::zonal(layerToCompare1,rasterzonal,sum)
raster::zonal(layerToCompare2,rasterzonal,sum)
x <- (rbind(raster::zonal(layerToCompare1,rasterzonal,sum)[,2],
          raster::zonal(layerToCompare2,rasterzonal,sum)[,2]))
length(which(!is.nan(x[1,]) != !is.nan(x[2,]) ))

x <- x[,!is.nan(x[1,]) & !is.nan(x[2,]) ]
xt <- t(x)

wilcox.test(xt)

# compare little grids 
layerToCompare1 <-BW3_10m_R0X_winter
layerToCompare2 <-BW3_200m_R0X_winter
x <- (cbind(values(layerToCompare1),values(layerToCompare2)))
head(x)
x <- x[!is.na(x[,1]) & !is.na(x[,2]), ]
wilcox.test(x)

par(mar=c(2, 4, 4, 2))
barplot(x, beside=T, col=c("orange","black"), width=c(5,5),  border=NA, xpd=T, ylab="Number of particules per big-grid", main= paste("Compare winter-winter 10m BW3 ",layer,"big-grids with no particles were removed"))
legend("topright",legend=c("winter","winter"), col=c("orange","black"), pch=20)

par(mar=c(2, 4, 4, 2))
par(mfrow=c(1,2))
plot(layerToCompare1,zlim=c(0, 25000))
points(continent,type="l")
plot(layerToCompare2,zlim=c(0, 25000))
points(continent,type="l")



### MAKE PLOTS 
#--------------------------------
stackR0X_BW3_ALLareas_allyear <- stack(BW3_10m_R03_allyear,BW3_10m_R03_allyear,BW3_10m_R03_allyear,BW3_10m_R03_allyear,BW3_10m_R03_allyear,BW3_10m_R03_allyear)

plot(BW3_10m_ALLareas_allyear)
plot(BW3_10m_R03_allyear,zlim=c(0, 25000));points(continent,type="l")
plot(BW3_10m_R03_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R03_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R03_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R03_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R03_allyear, add=T,zlim=c(0, 25000))

plot(BW3_10m_R03_allyear)
plot(BW3_10m_R03_allyear)
plot(BW3_10m_R03_allyear)
plot(BW3_10m_R03_allyear)
plot(BW3_10m_R03_allyear)

writeRaster(BW3_10m_ALLareas_allyear,"results_maps/tiff/BW3_10m_ALLareas_allyear.tiff")




