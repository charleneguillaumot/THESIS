library(raster)





BW6_10m_R0X_JFM_2008_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2008_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2009_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2009_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2010_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2010_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2011_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2011_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2012_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2012_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2013_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2013_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2014_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2014_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2015_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2015_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]

BW6_10m_R0X_JFM_2016_csv <- read.csv("results_maps/nouveaux/BW4-5-6/Number_of_part_BW6-001_release_all_01-02-03_year_2016_month_01-02-03-04-05-06-07-08-09-10-11-12_depth_-3000-5000_depth_at_start_-10_zone_R-03-R-04-R-05-R-06-R-07-R-08.csv")[,-1]


#BW6_10m_R0X_JFM_2008 <-raster("maps_for_figures/BW6_10m_R0X_JFM_2008.tif")


rastervide <- raster(nrows=240, ncol=480,resolution=c(0.25,0.125),xmn=-110,xmx=10,ymn=-80,ymx=-50)

BW6_10m_R0X_JFM_2008 <-rasterize(BW6_10m_R0X_JFM_2008_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2008_csv[,3])
BW6_10m_R0X_JFM_2009 <-rasterize(BW6_10m_R0X_JFM_2009_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2009_csv[,3])
BW6_10m_R0X_JFM_2010 <-rasterize(BW6_10m_R0X_JFM_2010_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2010_csv[,3])
BW6_10m_R0X_JFM_2011 <-rasterize(BW6_10m_R0X_JFM_2011_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2011_csv[,3])
BW6_10m_R0X_JFM_2012 <-rasterize(BW6_10m_R0X_JFM_2012_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2012_csv[,3])
BW6_10m_R0X_JFM_2013 <-rasterize(BW6_10m_R0X_JFM_2013_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2013_csv[,3])
BW6_10m_R0X_JFM_2014 <-rasterize(BW6_10m_R0X_JFM_2014_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2014_csv[,3])
BW6_10m_R0X_JFM_2015 <-rasterize(BW6_10m_R0X_JFM_2015_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2015_csv[,3])
BW6_10m_R0X_JFM_2016 <-rasterize(BW6_10m_R0X_JFM_2016_csv[,c(1,2)],rastervide,BW6_10m_R0X_JFM_2016_csv[,3])


## PLOT FREQUENCE
#-----------------
stack_BW6_10m_R0X_all_season_2008 <- stack(BW6_10m_R0X_AMJ_2008,BW6_10m_R0X_JAS_2008,BW6_10m_R0X_OND_2008,BW6_10m_R0X_JFM_2008)
BW6_10m_R0X_all_season_2008 <- calc(stack_BW6_10m_R0X_all_season_2008,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2008)
#BW6_10m_R0X_all_season_2008_re <- reclassify(BW6_10m_R0X_all_season_2008, cbind(NA, -99))
BW6_10m_R0X_all_season_2008[values(BW6_10m_R0X_all_season_2008)>=0] <- 1
plot(BW6_10m_R0X_all_season_2008)

stack_BW6_10m_R0X_all_season_2009 <- stack(BW6_10m_R0X_AMJ_2009,BW6_10m_R0X_JAS_2009,BW6_10m_R0X_OND_2009,BW6_10m_R0X_JFM_2009)
BW6_10m_R0X_all_season_2009 <- calc(stack_BW6_10m_R0X_all_season_2009,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2009)
#BW6_10m_R0X_all_season_2009_re <- reclassify(BW6_10m_R0X_all_season_2009, cbind(NA, -99))
BW6_10m_R0X_all_season_2009[values(BW6_10m_R0X_all_season_2009)>=0] <- 1
plot(BW6_10m_R0X_all_season_2009)

stack_BW6_10m_R0X_all_season_2010 <- stack(BW6_10m_R0X_AMJ_2010,BW6_10m_R0X_JAS_2010,BW6_10m_R0X_OND_2010,BW6_10m_R0X_JFM_2010)
BW6_10m_R0X_all_season_2010 <- calc(stack_BW6_10m_R0X_all_season_2010,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2010)
#BW6_10m_R0X_all_season_2010_re <- reclassify(BW6_10m_R0X_all_season_2010, cbind(NA, -99))
BW6_10m_R0X_all_season_2010[values(BW6_10m_R0X_all_season_2010)>=0] <- 1
plot(BW6_10m_R0X_all_season_2010)

stack_BW6_10m_R0X_all_season_2011 <- stack(BW6_10m_R0X_AMJ_2011,BW6_10m_R0X_JAS_2011,BW6_10m_R0X_OND_2011,BW6_10m_R0X_JFM_2011)
BW6_10m_R0X_all_season_2011 <- calc(stack_BW6_10m_R0X_all_season_2011,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2011)
#BW6_10m_R0X_all_season_2011_re <- reclassify(BW6_10m_R0X_all_season_2011, cbind(NA, -99))
BW6_10m_R0X_all_season_2011[values(BW6_10m_R0X_all_season_2011)>=0] <- 1
plot(BW6_10m_R0X_all_season_2011)

stack_BW6_10m_R0X_all_season_2012 <- stack(BW6_10m_R0X_AMJ_2012,BW6_10m_R0X_JAS_2012,BW6_10m_R0X_OND_2012,BW6_10m_R0X_JFM_2012)
BW6_10m_R0X_all_season_2012 <- calc(stack_BW6_10m_R0X_all_season_2012,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2012)
#BW6_10m_R0X_all_season_2012_re <- reclassify(BW6_10m_R0X_all_season_2012, cbind(NA, -99))
BW6_10m_R0X_all_season_2012[values(BW6_10m_R0X_all_season_2012)>=0] <- 1
plot(BW6_10m_R0X_all_season_2012)

stack_BW6_10m_R0X_all_season_2013 <- stack(BW6_10m_R0X_AMJ_2013,BW6_10m_R0X_JAS_2013,BW6_10m_R0X_OND_2013,BW6_10m_R0X_JFM_2013)
BW6_10m_R0X_all_season_2013 <- calc(stack_BW6_10m_R0X_all_season_2013,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2013)
#BW6_10m_R0X_all_season_2013_re <- reclassify(BW6_10m_R0X_all_season_2013, cbind(NA, -99))
BW6_10m_R0X_all_season_2013[values(BW6_10m_R0X_all_season_2013)>=0] <- 1
plot(BW6_10m_R0X_all_season_2013)


stack_BW6_10m_R0X_all_season_2014 <- stack(BW6_10m_R0X_AMJ_2014,BW6_10m_R0X_JAS_2014,BW6_10m_R0X_OND_2014,BW6_10m_R0X_JFM_2014)
BW6_10m_R0X_all_season_2014 <- calc(stack_BW6_10m_R0X_all_season_2014,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2014)
#BW6_10m_R0X_all_season_2014_re <- reclassify(BW6_10m_R0X_all_season_2014, cbind(NA, -99))
BW6_10m_R0X_all_season_2014[values(BW6_10m_R0X_all_season_2014)>=0] <- 1
plot(BW6_10m_R0X_all_season_2014)

stack_BW6_10m_R0X_all_season_2015 <- stack(BW6_10m_R0X_AMJ_2015,BW6_10m_R0X_JAS_2015,BW6_10m_R0X_OND_2015,BW6_10m_R0X_JFM_2015)
BW6_10m_R0X_all_season_2015 <- calc(stack_BW6_10m_R0X_all_season_2015,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2015)
#BW6_10m_R0X_all_season_2015_re <- reclassify(BW6_10m_R0X_all_season_2015, cbind(NA, -99))
BW6_10m_R0X_all_season_2015[values(BW6_10m_R0X_all_season_2015)>=0] <- 1
plot(BW6_10m_R0X_all_season_2015)

stack_BW6_10m_R0X_all_season_2016 <- stack(BW6_10m_R0X_AMJ_2016,BW6_10m_R0X_JAS_2016,BW6_10m_R0X_OND_2016,BW6_10m_R0X_JFM_2016)
BW6_10m_R0X_all_season_2016 <- calc(stack_BW6_10m_R0X_all_season_2016,mean, na.rm=T)
plot(BW6_10m_R0X_all_season_2016)
#BW6_10m_R0X_all_season_2016_re <- reclassify(BW6_10m_R0X_all_season_2016, cbind(NA, -99))
BW6_10m_R0X_all_season_2016[values(BW6_10m_R0X_all_season_2016)>=0] <- 1
plot(BW6_10m_R0X_all_season_2016)

stack_layers <- stack(BW6_10m_R0X_all_season_2008,BW6_10m_R0X_all_season_2009,BW6_10m_R0X_all_season_2010,BW6_10m_R0X_all_season_2011,BW6_10m_R0X_all_season_2012,BW6_10m_R0X_all_season_2013,BW6_10m_R0X_all_season_2014,BW6_10m_R0X_all_season_2015,BW6_10m_R0X_all_season_2016)
stack_layers_freq <- calc(stack_layers, sum, na.rm=T)
plot(stack_layers_freq)
writeRaster(stack_layers_freq, "stack_layers_freq_BW6.tiff")

#### look for the upper variability year
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "Oranges")
#my.palette.blue <- rev(brewer.pal(n = 9, name = "Blues"))
palettecol <-colorRampPalette(c("deepskyblue", "darkseagreen","lightgreen","green","yellow","gold","orange", "red","firebrick"))(100)

plot(stack(BW6_10m_R0X_all_season_2008,BW6_10m_R0X_all_season_2009,BW6_10m_R0X_all_season_2010,BW6_10m_R0X_all_season_2011,BW6_10m_R0X_all_season_2012), col=palettecol)
plot(stack(BW6_10m_R0X_all_season_2013,BW6_10m_R0X_all_season_2014,BW6_10m_R0X_all_season_2015,BW6_10m_R0X_all_season_2016), col=palettecol)
plot(BW6_10m_R0X_all_season_2008, col=palettecol, main="2008")
plot(BW6_10m_R0X_all_season_2009, col=palettecol, main="2009")
plot(BW6_10m_R0X_all_season_2010, col=palettecol, main="2010")
plot(BW6_10m_R0X_all_season_2011, col=palettecol, main="2011")
plot(BW6_10m_R0X_all_season_2012, col=palettecol, main="2012")
plot(BW6_10m_R0X_all_season_2013, col=palettecol, main="2013")
plot(BW6_10m_R0X_all_season_2014, col=palettecol, main="2014")
plot(BW6_10m_R0X_all_season_2015, col=palettecol, main="2015")
plot(BW6_10m_R0X_all_season_2016, col=palettecol, main="2016")

