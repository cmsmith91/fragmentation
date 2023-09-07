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

