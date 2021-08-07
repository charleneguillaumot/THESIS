library(raster)
library(SpaDES.tools)
continent <- read.csv("results_maps/worldmap-0.01.csv")


BW3_10m_R03_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03.csv")[,-1]
BW3_10m_R04_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-04.csv")[,-1]
BW3_10m_R05_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-05.csv")[,-1]
BW3_10m_R06_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-06.csv")[,-1]
BW3_10m_R07_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-07.csv")[,-1]
BW3_10m_R08_summer_csv <- read.csv("results_maps/PER_RELEASE/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-08.csv")[,-1]

#BW1_10m_R0X_summer_csv <- read.csv("results_maps/Number_of_part_BW1-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

#BW2_10m_R0X_summer_csv <- read.csv("results_maps/Number_of_part_BW2-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

#BW3_10m_R0X_summer_csv <- read.csv("results_maps/Number_of_part_BW3-001_release_all_01-02-03_year_all_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]


#### Transform the csv into rasters
#------------------------------------
rastervide <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)

BW3_10m_R03_summer <-rasterize(BW3_10m_R03_summer_csv[,c(1,2)],rastervide,BW3_10m_R03_summer_csv[,3])
BW3_10m_R04_summer <-rasterize(BW3_10m_R04_summer_csv[,c(1,2)],rastervide,BW3_10m_R04_summer_csv[,3])
BW3_10m_R05_summer <-rasterize(BW3_10m_R05_summer_csv[,c(1,2)],rastervide,BW3_10m_R05_summer_csv[,3])
BW3_10m_R06_summer <-rasterize(BW3_10m_R06_summer_csv[,c(1,2)],rastervide,BW3_10m_R06_summer_csv[,3])
BW3_10m_R07_summer <-rasterize(BW3_10m_R07_summer_csv[,c(1,2)],rastervide,BW3_10m_R07_summer_csv[,3])
BW3_10m_R08_summer <-rasterize(BW3_10m_R08_summer_csv[,c(1,2)],rastervide,BW3_10m_R08_summer_csv[,3])

#BW1_10m_R0X_summer <-rasterize(BW1_10m_R0X_summer_csv[,c(1,2)],rastervide,BW1_10m_R0X_summer_csv[,3])
#BW2_10m_R0X_summer <-rasterize(BW2_10m_R0X_summer_csv[,c(1,2)],rastervide,BW2_10m_R0X_summer_csv[,3])
#BW3_10m_R0X_summer <-rasterize(BW3_10m_R0X_summer_csv[,c(1,2)],rastervide,BW3_10m_R0X_summer_csv[,3])

par(mfrow=c(1,3))
plot(BW3_10m_R03_summer)
points(continent,type="l")

plot(BW3_10m_R04_summer)
points(continent,type="l")

plot(BW3_10m_R05_summer)
points(continent,type="l")

# Transform into ascii and plot on Qgis 
writeRaster(BW3_10m_R03_summer,"prod_maps_results/BW3/summer/BW3_10m_R03_summer_Depthparticles.tif")
writeRaster(BW3_10m_R04_summer,"prod_maps_results/BW3/summer/BW3_10m_R04_summer_Depthparticles.tif")
writeRaster(BW3_10m_R05_summer,"prod_maps_results/BW3/summer/BW3_10m_R05_summer_Depthparticles.tif")
writeRaster(BW3_10m_R06_summer,"prod_maps_results/BW3/summer/BW3_10m_R06_summer_Depthparticles.tif")
writeRaster(BW3_10m_R07_summer,"prod_maps_results/BW3/summer/BW3_10m_R07_summer_Depthparticles.tif")
writeRaster(BW3_10m_R08_summer,"prod_maps_results/BW3/summer/BW3_10m_R08_summer_Depthparticles.tif")


#### Transform into big grids
#---------------------------
# BW1_10m_R03_spring_big <- aggregate(BW1_10m_R03_spring, fact=25, fun='sum')
# BW1_10m_R04_spring_big <- aggregate(BW1_10m_R04_spring, fact=25, fun='sum')
# BW1_10m_R05_spring_big <- aggregate(BW1_10m_R05_spring, fact=25, fun='sum')
# BW1_10m_R06_spring_big <- aggregate(BW1_10m_R06_spring, fact=25, fun='sum')
# BW1_10m_R07_spring_big <- aggregate(BW1_10m_R07_spring, fact=25, fun='sum')
# BW1_10m_R08_spring_big <- aggregate(BW1_10m_R08_spring, fact=25, fun='sum')

BW1_10m_R0X_summer_big <- aggregate(BW1_10m_R0X_summer, fact=25, fun='sum')
BW2_10m_R0X_summer_big <- aggregate(BW2_10m_R0X_summer, fact=25, fun='sum')
BW3_10m_R0X_summer_big <- aggregate(BW3_10m_R0X_summer, fact=25, fun='sum')

BW1_10m_R0X_summer_big <- aggregate(BW1_10m_R0X_summer, fact=25, fun='mean')
BW2_10m_R0X_summer_big <- aggregate(BW2_10m_R0X_summer, fact=25, fun='mean')
BW3_10m_R0X_summer_big <- aggregate(BW3_10m_R0X_summer, fact=25, fun='mean')

par(mfrow=c(1,2))
plot(BW1_10m_R0X_summer)
points(continent,type="l")

plot(BW1_10m_R0X_summer_big)
plot(rasterToPolygons(BW1_10m_R0X_summer_big), add=TRUE, border='black', lwd=0.2) 
points(continent,type="l")

# stacker tout 
#--------------
stack_10m_spring_big <- stack(BW1_10m_R03_spring_big,BW1_10m_R04_spring_big,BW1_10m_R05_spring_big,BW1_10m_R06_spring_big,BW1_10m_R07_spring_big,BW1_10m_R08_spring_big)
#sum_10m_spring_big <- calc(stack_10m_spring_big, sum, na.rm=T)
plot(sum_10m_spring_big)
plot(rasterToPolygons(sum_10m_spring_big), add=TRUE, border='black', lwd=0.2) 
points(continent,type="l")
#tableau <- rasterToPoints(stack(sum_10m_spring_big,stack_10m_spring_big))
#tableau


stack_all <- stack(BW1_10m_R0X_summer_big,BW2_10m_R0X_summer_big,BW3_10m_R0X_summer_big)
tableau <- rasterToPoints(stack_all)
head(tableau)

## Boxplot 
tableau_modif <- tableau[,-c(1,2)]
colnames(tableau_modif) <- c("BW1","BW2","BW3")
head(t(tableau_modif))
barplot(t(tableau_modif), beside=T)
tableau_modif

barplot(t(tableau_modif), beside=T, col=c("blue","red","yellow"))



tableau_gg2 <- c(tableau_modif[,1],tableau_modif[,2],tableau_modif[,3])
id <- c(rep("BW1", nrow(tableau_modif)),rep("BW2", nrow(tableau_modif)),rep("BW3", nrow(tableau_modif)))
tableau_gg <- as.data.frame(cbind(tableau_gg2, id, seq (1,length(tableau_gg2),1)))
tableau_gg
colnames(tableau_gg) <- c("nb_part","BW","x")

library(ggplot2)

ggplot(data=tableau_gg, aes(x=x,y=nb_part, fill=BW))+
  geom_bar(stat="identity")


# se faire une ptite grille vierge avec les numéros de pixels 
raster_vierge <- BW1_10m_R03_spring_big; values(raster_vierge) <- seq (1, 200, 1)
plot(raster_vierge)
plot(rasterToPolygons(raster_vierge), add=TRUE, border='black', lwd=0.2) 
points(continent,type="l")
