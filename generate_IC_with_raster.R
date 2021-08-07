  setwd("/home/RBINS.BE/vduliere/R/PART/BALLAST/BALLAST_IC/")
  # 
  
  rm(list=ls()) # clear all previously defined variables
  
  # les librairie 
  library(ncdf4)
  library(maps)
  library(ggplot2)
  library(ggmap)
  library(RColorBrewer)
  library(bazar)
  library(raster)
  library(reshape2)
  library(plyr)
  library(data.table)
  library(raster)
  
  
  library (rgdal)
  library (rgeos)
  library(maptools)
  library(tmap)
  
  #---------------------------------------------------
  # Function
  #---------------------------------------------------
  "%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0 
  #--------------------------------------------------
  
  file='Minimum_distance2.nc'
  
  nm_to_m<-1.852 # conversion from nautical miles to meters
  dist1_nm=40#2#40#190  * nm_to_m # nm
  dist2_nm=60#20#60#210 
  
  lon1=-82
  #lon1=-90
  lon2=-10
  #lon2=20
  min_depth=200#200#10 #200
  div_lat=6
  div_lon=4
  
  Rt<-6356.800 # earth radius
  fact1=3.14159265*Rt/180.
  fact2=3.14159265/180.
  dist1=dist1_nm*nm_to_m
  dist2=dist2_nm*nm_to_m
  nc <- nc_open(file)
  lat <- ncvar_get(nc,'latitude')
  i_lat<-array(1:length(lat))
  lon <- ncvar_get(nc,'longitude')
  i_lon<-array(1:length(lon))
  depth <- ncvar_get(nc,'depth')
  dist_min <- ncvar_get(nc,'min_dist')
  nc_close(nc)
  
  lat_temp=depth
  lon_temp=depth
  i_lat_temp=depth
  i_lon_temp=depth
  
  for (it in (1:length(lon))) {lat_temp[it,]<-lat[]}
  for (it in (1:length(lon))) {i_lat_temp[it,]<-i_lat[]}
  for (it in (1:length(lat))) {lon_temp[,it]<-lon[]}
  for (it in (1:length(lat))) {i_lon_temp[,it]<-i_lon[]}
  
  
  data2<-data.frame(lat=c(lat_temp),lon=c(lon_temp),dist_min=c(dist_min),depth=c(depth),i_index=c(i_lat_temp),j_index=c(i_lon_temp))
  
  data2<-data2[data2$lon>lon1,]
  data2<-data2[data2$lon<lon2,]
  data2<-data2[data2$depth>min_depth,]
  data2<-data2[data2$dist_min>dist1,]
  data2<-data2[data2$dist_min<dist2,]
  #data2<-data2[((data2$i_index+data2$j_index) %% 4 ==0),]
  data2<-data2[(data2$i_index %% div_lon ==0),]
  data2<-data2[(data2$j_index %% div_lat ==0),]
  print(length(data2[,1]))
  

  # South America
  lat_point_zone<-c(-50,-46,-54.5,-56,-56,-56,-52)
  lon_point_zone<-c(-68,-67,-67,-64,-66,-71,-59)
  point_zone<-c('R-01','R-01','R-01','R-01','R-01','R-01','R-02')
  # WAP
  lat_point_zone<-append(lat_point_zone,c(-71,-70,-62,-63,-65,-63,-68.5,-61,-61))
  lon_point_zone<-append(lon_point_zone,c(-78,-75,-59,-60,-67,-63,-68,-60,-55))
  point_zone<-append(point_zone,c('R-03','R-03','R-03','R-03','R-03','R-03','R-03','R-03','R-03'))
  # EAP
  lat_point_zone<-append(lat_point_zone,c(-63,-65,-67,-68,-70,-72,-70,-71.6,-72.3,-73.9,-71.4))
  lon_point_zone<-append(lon_point_zone,c(-53,-56,-59,-59,-58,-58,-57,-58,-58.4,-60,-55))
  point_zone<-append(point_zone,c('R-04','R-04','R-04','R-04','R-04','R-04','R-04','R-04','R-04','R-04','R-04'))
  # Weddell Sea
  lat_point_zone<-append(lat_point_zone,c(-74.6,-75,-73,-70)) #,-73,-73,-70))
  lon_point_zone<-append(lon_point_zone,c(-60,-53,-40,-11)) #,-39,-22,-11))
  point_zone<-append(point_zone,c('R-05','R-05','R-05','R-05'))#,'R-05','R-05','R-05'))
  # the 3 islands areas (South Georgia and South Sandwich Island (R-06); Coronation Island (R-07); Scotia Arc (R-08))
  lat_point_zone<-append(lat_point_zone,c(-53.6,-54.2,-60,-61,-58))
  lon_point_zone<-append(lon_point_zone,c(-42,-37,-45,-44,-26))
  point_zone<-append(point_zone,c('R-06','R-06','R-07','R-07','R-08'))
  
  dist_zone=array(100000.,dim=c(length(data2$lon),length(lat_point_zone)))
  
  for (it in (1:length(lat_point_zone))){
    for (ip in (1:length(data2$lat))){
      dlat=(abs(data2$lat[ip]-lat_point_zone[it]))*fact1
      dlat=dlat*dlat
            
      dlon=(abs(data2$lon[ip]-lon_point_zone[it]))*fact1*cos(data2$lat[ip]*fact2)
      dlon=dlon*dlon
      dist_zone[ip,it]=sqrt(dlat+dlon)
      
    }
  }
  for (ip in (1:length(data2$lat))){
    data2$zone[ip]=point_zone[which(dist_zone[ip,]==min(dist_zone[ip,]))]
    data2$state[ip]=2
  }
  
  # remove zone R-01 and R-02
  data2<-data2[data2$zone!="R-01",]
  data2<-data2[data2$zone!="R-02",]
  
  
  label=array(1:length(data2[,1]))
  out=data.frame(zone=data2$zone,label=label,data2$lat,data2$lon)
  csvpath <- ""
  csvname <- paste("Release_location_",dist1_nm,"-",dist2_nm,"_every_",div_lat,'-',div_lon,"_mindepth_",min_depth,"_AP_noSA.csv",sep="")
  csvfile <- paste(csvpath, csvname, sep="")
  write.table(na.omit(out),csvfile, row.names=FALSE, sep=",")
  
  lonmin <- -90
  lonmax <- 10
  latmin <- -80
  latmax <- -45
  lat_orient = -90.
  lon_orient = -45.
  
  lonlim <- c(lonmin,lonmax)
  latlim <- c(latmin,latmax)
  
  
  # world <- map_data("world")
  # worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
  #   geom_path() +
  #   scale_y_continuous(breaks=(-2:2) * 30) +
  #   scale_x_continuous(breaks=(-4:4) * 45) 
  # worldmap_tot_first <- worldmap + coord_map("ortho",
  #   orientation=c(lat_orient, lon_orient, 0),xlim=lonlim,ylim=latlim) +
  #   geom_count(data=data2,aes(x=data2$lon,y=data2$lat,colour=data2$zone, #"red",
  #                          group=data2$state))+scale_size_area(max_size = 1) +
  #   scale_colour_hue(c(200,270))# + scale_colour_gradientn(colours=rainbow(4)) #
  # worldmap_tot_first
  # ggsave(paste("Release_location_",dist1_nm,"-",dist2_nm,"_every_",div_lat,'-',div_lon,"_mindepth_",min_depth,"_AP_noSA.jpeg",sep=""),width = 30,
  #        height = 24, dpi=100, units=c("cm"))
  # dev.off()

  dir="~/R/PART/BALLAST/SHP/MPA_split_raatd/"
  PG <- readShapePoly(paste(dir,"proposed_mpa.shp",sep=""))
PG2<- readShapePoly(paste(dir,"mpa_AP_proposal.shp",sep=""))
file=paste(dir,"WDPA_RAATD_marine-shapefile-polygons_NAME_Yorktown.shp",sep="")
PG <-readShapePoly(paste(dir,"WDPA_RAATD_marine-shapefile-polygons_NAME_Yorktown.shp",sep=""))
PG <- readShapePoly(paste(dir,"WDPA_RAATD_marine-shapefile-polygons_NAME_South\ Orkney\ Islands\ Southern\ Shelf\ Marine\ Protected\ Area.dbf",sep=""))
  library(ggplot2)
  library (rgdal)
  library (rgeos)
  library(maptools)
  #PG<-readOGR(paste(dir,"proposed_mpa.shp",sep=""),layer=mpa)
  AG <- fortify(PG)


  world <- map_data("world")
  worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_path() +
    scale_y_continuous(breaks=(-2:2) * 30) +
    scale_x_continuous(breaks=(-4:4) * 45) 
  worldmap_tot_first <- worldmap + coord_map("ortho",
    orientation=c(lat_orient, lon_orient, 0),xlim=lonlim,ylim=latlim) +
    geom_count(data=data2,aes(x=data2$lon,y=data2$lat,colour=data2$zone,group=data2$state))+scale_size_area(max_size = 1) + 
    
    scale_colour_hue(c(200,270))+
    geom_polygon(data=AG, aes(long, lat, group = group), 
      colour = alpha("darkgray", 1/2), size = 0.7, fill = 'skyblue', alpha = .3)# + scale_colour_gradientn(colours=rainbow(4)) # 
  worldmap_tot_first
  ggsave(paste("Release_location_",dist1_nm,"-",dist2_nm,"_every_",div_lat,'-',div_lon,"_mindepth_",min_depth,"_AP_noSA.jpeg",sep=""),width = 30,
         height = 24, dpi=100, units=c("cm"))
  dev.off()


  
