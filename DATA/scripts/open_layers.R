library(raster)
library(SpaDES.tools)
continent <- read.csv("results_maps/worldmap-0.01.csv")

#shape <- raster::shapefile("Proposed MPA/mpa_AP_proposal.shp")
#crs(shape)<-crs(x)

rasterzonal <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)
rasterzonal2 <- SpaDES.tools::splitRaster(rasterzonal,nx=4,ny=24)
for (i in 1:96){
  values(rasterzonal2[[i]]) <- i
}
rasterzonal <- mergeRaster(rasterzonal2)
plot(rasterzonal)
points(continent,type="l")

#---------------------
# NUMBER OF PARTICLES 
#---------------------
BW3_10m_R03_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]
BW3_10m_R04_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-04.csv")[,-1]
BW3_10m_R05_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-05.csv")[,-1]
BW3_10m_R06_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-06.csv")[,-1]
BW3_10m_R07_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-07.csv")[,-1]
BW3_10m_R08_spring_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_07-08-09_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-08.csv")[,-1]

BW3_10m_R03_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]
BW3_10m_R04_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-04.csv")[,-1]
BW3_10m_R05_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-05.csv")[,-1]
BW3_10m_R06_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-06.csv")[,-1]
BW3_10m_R07_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-07.csv")[,-1]
BW3_10m_R08_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_10-11-12_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-08.csv")[,-1]

#### Transform the csv into rasters
#------------------------------------
rastervide <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)

BW3_10m_R03_spring <-rasterize(BW3_10m_R03_spring_csv[,c(1,2)],rastervide,BW3_10m_R03_spring_csv[,3])
BW3_10m_R04_spring <-rasterize(BW3_10m_R04_spring_csv[,c(1,2)],rastervide,BW3_10m_R04_spring_csv[,3])
BW3_10m_R05_spring <-rasterize(BW3_10m_R05_spring_csv[,c(1,2)],rastervide,BW3_10m_R05_spring_csv[,3])
BW3_10m_R06_spring <-rasterize(BW3_10m_R06_spring_csv[,c(1,2)],rastervide,BW3_10m_R06_spring_csv[,3])
BW3_10m_R07_spring <-rasterize(BW3_10m_R07_spring_csv[,c(1,2)],rastervide,BW3_10m_R07_spring_csv[,3])
BW3_10m_R08_spring <-rasterize(BW3_10m_R08_spring_csv[,c(1,2)],rastervide,BW3_10m_R08_spring_csv[,3])

BW3_10m_R03_summer <-rasterize(BW3_10m_R03_summer_csv[,c(1,2)],rastervide,BW3_10m_R03_summer_csv[,3])
BW3_10m_R04_summer <-rasterize(BW3_10m_R04_summer_csv[,c(1,2)],rastervide,BW3_10m_R04_summer_csv[,3])
BW3_10m_R05_summer <-rasterize(BW3_10m_R05_summer_csv[,c(1,2)],rastervide,BW3_10m_R05_summer_csv[,3])
BW3_10m_R06_summer <-rasterize(BW3_10m_R06_summer_csv[,c(1,2)],rastervide,BW3_10m_R06_summer_csv[,3])
BW3_10m_R07_summer <-rasterize(BW3_10m_R07_summer_csv[,c(1,2)],rastervide,BW3_10m_R07_summer_csv[,3])
BW3_10m_R08_summer <-rasterize(BW3_10m_R08_summer_csv[,c(1,2)],rastervide,BW3_10m_R08_summer_csv[,3])

# merge
BW3_10m_ALLareas_spring <- merge(BW3_10m_R03_spring,BW3_10m_R04_spring,BW3_10m_R05_spring,BW3_10m_R06_spring,BW3_10m_R07_spring,BW3_10m_R08_spring)

BW3_10m_ALLareas_summer <- merge(BW3_10m_R03_summer,BW3_10m_R04_summer,BW3_10m_R05_summer,BW3_10m_R06_summer,BW3_10m_R07_summer,BW3_10m_R08_summer)


#### Schoener simple comparison
#------------------------------
library(dismo)
nicheOverlap(BW3_10m_R03_spring,BW3_10m_R03_summer, stat="D")
nicheOverlap(BW3_10m_R04_spring,BW3_10m_R04_summer, stat="D")
nicheOverlap(BW3_10m_R05_spring,BW3_10m_R05_summer, stat="D")
nicheOverlap(BW3_10m_R06_spring,BW3_10m_R06_summer, stat="D")
nicheOverlap(BW3_10m_R07_spring,BW3_10m_R07_summer, stat="D")
nicheOverlap(BW3_10m_R08_spring,BW3_10m_R08_summer, stat="D")
nicheOverlap(BW3_10m_ALLareas_spring,BW3_10m_ALLareas_summer, stat="D")

# Spearman correlation 
cor(na.omit(cbind(values(BW3_10m_R03_spring),values(BW3_10m_R03_summer))), method = "spearman")
cor(na.omit(cbind(values(BW3_10m_R04_spring),values(BW3_10m_R04_summer))), method = "spearman")


### BIG GRIDS ET WILCOXON 
#--------------------------------
layerToCompare1 <-BW3_10m_R08_spring
layerToCompare2 <-BW3_10m_R08_summer
layer <- "ALLareas"
raster::zonal(layerToCompare1,rasterzonal,mean)
raster::zonal(layerToCompare2,rasterzonal,mean)
x <- (rbind(raster::zonal(layerToCompare1,rasterzonal,mean)[,2],
          raster::zonal(layerToCompare2,rasterzonal,mean)[,2]))
length(which(!is.nan(x[1,]) != !is.nan(x[2,]) ))

x <- x[,!is.nan(x[1,]) & !is.nan(x[2,]) ]
xt <- t(x)

wilcox.test(xt)

par(mar=c(2, 4, 4, 2))
barplot(x, beside=T, col=c("orange","black"), width=c(5,5),  border=NA, xpd=T, ylab="Number of particules per big-grid", main= paste("Compare spring-summer 10m BW3 ",layer,"big-grids with no particles were removed"))
legend("topright",legend=c("spring","summer"), col=c("orange","black"), pch=20)

par(mar=c(2, 4, 4, 2))
par(mfrow=c(1,2))
plot(layerToCompare1,zlim=c(0, 25000))
points(continent,type="l")
plot(layerToCompare2,zlim=c(0, 25000))
points(continent,type="l")

### MAKE PLOTS 
#--------------------------------
stackR0X_BW3_ALLareas_allyear <- stack(BW3_10m_R03_allyear,BW3_10m_R04_allyear,BW3_10m_R05_allyear,BW3_10m_R06_allyear,BW3_10m_R07_allyear,BW3_10m_R08_allyear)

plot(BW3_10m_ALLareas_allyear)
plot(BW3_10m_R03_allyear,zlim=c(0, 25000));points(continent,type="l")
plot(BW3_10m_R04_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R05_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R06_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R07_allyear, add=T,zlim=c(0, 25000))
plot(BW3_10m_R08_allyear, add=T,zlim=c(0, 25000))

plot(BW3_10m_R04_allyear)
plot(BW3_10m_R05_allyear)
plot(BW3_10m_R06_allyear)
plot(BW3_10m_R07_allyear)
plot(BW3_10m_R08_allyear)

writeRaster(BW3_10m_ALLareas_allyear,"results_maps/tiff/BW3_10m_ALLareas_allyear.tiff")




