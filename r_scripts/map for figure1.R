rm(list=ls())
# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
# load packages
library(rgdal)  
library (raster)   #
library(sf)
library(maps)
library(tidyverse)
library(mapdata)
library(stars)
select=dplyr::select
map=maps::map

northeast_us=tolower(c("Maine", "Vermont", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New York", "Pennsylvania", 
"New Jersey", "Delaware",  "Maryland"))
states <- map_data("state") 
roi = states[states$region %in% northeast_us,]
bounds=c(min(roi$long)-.2,min(roi$lat)-.2,max(roi$long)+.2,max(roi$lat)+.2)
xmin = min(bounds[c(3,1)]); xmax = max(bounds[c(3,1)])
ymin = min(bounds[c(4,2)]); ymax= max(bounds[c(4,2)])

#load nlcd data
img_path="~/Dropbox/Fall_2014/postdoc/inat project/data_current/NLCD_2016_LandCover/NLCD_2016_Land_Cover_L48_20190424.img"
nlcd=read_stars(img_path)
# test=nlcd_st[,120000:161190,1:52212]
# plot(test)

#load the data for the map of the northeastern US
test <- map("world2Hires",region=c("USA","Canada"),fill=T,plot=F)
test$x <- test$x-360
test2 = st_as_sf(test, crs = 4326)
state_sf=states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
grid_pts=read_sf("/Users/colleen/Dropbox/gis\ files/northeast_coords_within.shp")
grid_pts_crs=st_transform(grid_pts,crs = 4326)
lat_lon=data.frame(st_coordinates(grid_pts_crs)) %>% rename(longitude=X,latitude=Y)%>% 
    mutate(state=maps::map.where("state",longitude,latitude),
           country= maps::map.where("world",longitude,latitude))
#
df=lat_lon %>% filter(!is.na(country) & !is.na(state))

#make map
grid_sites=ggplot(data=test2) +
    geom_sf(fill='honeydew1') +
    geom_sf(data = states, fill = NA)+
    coord_sf(crs = st_crs(4326),
             xlim = c(-82,-66.7), 
             ylim = c(37.5,47.5), expand = FALSE)+
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Longitude") + 
    ylab("Latitude")+
    geom_point(data = df, aes(y = latitude, x= longitude),size=.3)

grid_sites

#next pick one of the landscapes to 
check_me=grid_pts[1:10,]
ggplot(check_me)+geom_sf()
nrow(grid_pts)



#first with random pt
square = st_sfc(st_buffer(st_point(c(0, 1500000)), 400,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
plot(nlcd[square], reset = FALSE,main="")

## downsample set to c(0,0,1)
plot(square, col = NA, border = 'red', add = TRUE, lwd = 2)

#now with buffers
buffer600=read_sf("/Users/colleen/Dropbox/gis\ files/600m_buffer.shp")
buffer1500=read_sf("/Users/colleen/Dropbox/gis\ files/1500m_buffer.shp")
buffer2400=read_sf("/Users/colleen/Dropbox/gis\ files/2400m_buffer.shp")

i=500
focal600=buffer600[i,]
focal1500=buffer1500[i,]
focal2400=st_sfc(buffer2400, crs = st_crs(nlcd))
plot(nlcd[focal2400],reset=F,main="")


one_point=grid_pts[i,]
test_me=st_transform(st_buffer(one_point$geometry,600,endCapStyle="SQUARE",joinStyle="MITRE"),crs=st_crs(nlcd))
plot(nlcd[test_me])

plot(st_buffer(one_point$geometry,400))
plot(st_point(one_point))

square = st_sfc(st_buffer(st_point(c(0, 1500000)), 200,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
square = st_sfc(st_buffer(st_point(c(1801116,1500000)), 400,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
plot(nlcd[square], reset = FALSE,main="",axes=T)
plot(st_sfc(st_point(c(0, 1500000)),crs=st_crs(nlcd)),add=T,cex=3,pch=16)

plot(square,axes=T,add=T)
one_point$geometry
plot(nlcd,axes=T)
plot(st_sfc(st_point(c(0, 1500000)),crs=st_crs(nlcd)),add=T,cex=3,pch=16)
plot(st_sfc(st_point(c(100, 1500000)),crs=st_crs(nlcd)),add=T,cex=3,pch=16)

plot(st_transform(one_point$geometry,crs=st_crs(nlcd)),add=T,cex=13)

####old
lat_lon=data.frame(st_coordinates(grid_pts_crs)) %>% rename(longitude=X,latitude=Y)%>% 
    mutate(state=maps::map.where("state",longitude,latitude),
           country= maps::map.where("world",longitude,latitude))
#shape <- readOGR(dsn = "/Users/colleen/Dropbox/gis\ files/northeast_coords_within.shp")

buffer600=read_sf("/Users/colleen/Dropbox/gis\ files/600m_buffer.shp")
buffer1500=read_sf("/Users/colleen/Dropbox/gis\ files/1500m_buffer.shp")
buffer2400=read_sf("/Users/colleen/Dropbox/gis\ files/2400m_buffer.shp")

crs_600m=st_transform(buffer600,crs = 4326)


###get rid of obs in the wrong states and convert to a shape file
eastern_states = tolower(c("district of columbia",'alabama','connecticut','delaware',
                           "Florida", "Georgia", "Louisiana", "Maine", "Maryland", "Massachusetts", "Mississippi", "New Hampshire", "New Jersey", "New York", "North Carolina", "Pennsylvania", "Rhode Island", "South Carolina", "Vermont", "Virginia", "West Virginia"))
df=lat_lon %>% filter(!is.na(country) & !is.na(state))

shape <- readOGR(dsn = "/Users/colleen/Dropbox/gis\ files/600m_buffer.shp")
proj4string(shape)
a=CRS("+init=epsg:5070")
a@projargs==proj4string(shape)
proj4string(shape) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
range(shape@data$Latitude)
plot(shape)
ggplot()+geom_sf(data=crs_600m)#+
    coord_sf(crs = st_crs(4326),
             xlim =c(-75.5,-75), 
             ylim = c(38.5,42), expand = FALSE)
ggplot()+geom_sf(data=buffer600)+
    coord_sf(crs = st_crs(3347,
                          xlim = c(-75.5,-75), 
                          ylim = c(38.5,42), expand = FALSE))

ggplot()+
    geom_sf(data=buffer2400)+
    geom_sf(data=buffer1500)+
    geom_sf(data=buffer600)+
    coord_sf(xlim=c(1650000,1700000),ylim=c(2250000,2300000),expand=F)

ggplot()+geom_sf(data=buffer600)+
    coord_sf(xlim=c(1600000,1800000),ylim=c(2200000,2400000),datum=sf::st_crs(5070),expand=F)
ggplot()+geom_sf(data=buffer600)+
    coord_sf(datum=sf::st_crs(5070))

ggplot()+geom_sf(data=crs_600m)+
    coord_sf(xlim=c(-80,-79.5),ylim=c(40.5,41),expand=F)

#xmin: 1273437 ymin: 1835196 xmax: 2250019 ymax: 3004870
st_geometry(buffer600)

xmin_w=1273437
xmin=70
xmax_w=2250019
xmax=82

weird_xs_norm=(xmin_w:xmax_w-xmin_w)/(xmin_w)
val=76
x_norm=(val-xmin)/(xmin)
the_index=which(round(weird_xs_norm,6)==round(x_norm,6))
my_seq=xmin_w:xmax_w
my_seq[the_index]
   
xmin_w=1835196
xmin=70
xmax_w=3004870
xmax=82

weird_xs_norm=(xmin_w:xmax_w-xmin_w)/(xmin_w)
val=76
x_norm=(val-xmin)/(xmin)
the_index=which(round(weird_xs_norm,6)==round(x_norm,6))
my_seq=xmin_w:xmax_w
my_seq[the_index]

 


#make map
grid_sites_over_region=ggplot(data=test2) +
    geom_sf(fill='honeydew1') +
    geom_sf(data = states, fill = NA)+
    coord_sf(crs = st_crs(4326),
             xlim = c(-82,-66.7), 
             ylim = c(37.5,47.5), expand = FALSE)+
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Longitude") + 
    ylab("Latitude")+
    geom_point(data = df, aes(y = latitude, x= longitude),size=.3)

#next map showing how i measured forest edge at three different spatial scales
#nlcd data wiht forest in green, outlined in black to show forest edge
# and red squares for the different scales


ggplot(data=buffer1500)+
    geom_sf(fill='aliceblue')+
    geom_sf(data=buffer600)+
    coord_sf(crs=st_crs(4326),
             xlim = c(-75.5,-75), 
             ylim = c(38.5,39), expand = FALSE)

#make map of just nj
usa <- map("world2Hires",region=c("USA"),fill=T,plot=F)
ggplot(data=test2) +
    geom_sf(fill='antiquewhite') +
    geom_sf(data = states, fill = NA)+
    coord_sf(crs = st_crs(4326),
             xlim = c(-76,-72), 
             ylim = c(38.5,42), expand = FALSE)+
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Longitude") + 
    ylab("Latitude")+
    geom_point(data = df, aes(y = latitude, x= longitude),size=.3)

##add buffer
# ggplot(data=test2) +
#     geom_sf(fill='antiquewhite') +
#     geom_sf(data = states, fill = NA)+
#     coord_sf(crs = st_crs(4326),
#              xlim = c(-82,-66.7), 
#              ylim = c(37.5,47.5), expand = FALSE)+
#     theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     xlab("Longitude") + 
#     ylab("Latitude")+
#     geom_point(data = df, aes(y = latitude, x= longitude),size=.3)+
#     geom_sf(data=crs_600m)

