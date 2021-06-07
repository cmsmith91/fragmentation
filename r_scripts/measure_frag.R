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
library(furrr)
select=dplyr::select
map=maps::map

#plan: crop nlcd data to just be in the northeast US
#load nlcd data
img_path="~/Dropbox/Fall_2014/postdoc/inat project/data_current/NLCD_2016_LandCover/NLCD_2016_Land_Cover_L48_20190424.img"
nlcd=read_stars(img_path)

#get shapefile of northeastern US
northeast_us=tolower(c("Maine", "Vermont", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New York", "Pennsylvania", 
                       "New Jersey", "Delaware",  "Maryland"))

northeast_sf= st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
    filter(ID %in% northeast_us)
northeast_sf_crs=st_transform(northeast_sf,crs=st_crs(nlcd))

#crop the nlcd data to be just states in the northeastern US
nlcd_cropped=nlcd %>% st_crop(northeast_sf_crs)
plot(nlcd_cropped)

#now make a grid over northeast us shapefile
#transform to crs in meters: utm zone are in meters - most of northeast is in 18n
#espg=26918
northeast_sf_crs_meters=st_transform(northeast_sf,crs=26918)
# ggplot()+geom_sf(data=northeast_sf_crs_meters)

#make a grid - takes a while to run so try 
big_grid=st_make_grid(northeast_sf_crs_meters,5000,what = "centers") 
grid=big_grid[northeast_sf_crs_meters] %>% st_transform(crs=st_crs(nlcd))

#plot to make sure it looks ok
#just plot nj for now...
nj=northeast_sf_crs %>% filter(ID=='new jersey')
grid_nj=grid[nj]
ggplot(data=northeast_sf_crs %>% filter(ID=='new jersey')) + geom_sf(fill='lightblue')+
    geom_sf(data=grid_nj,size=.06)

#loop through each grid point and make a square buffer around it
grid[1]
square = st_sfc(st_buffer(grid[1], 1200,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
plot(nlcd[square], reset = FALSE,main="",axes=T)

for(i in 1:6){
    square = st_sfc(st_buffer(grid[i], 1200,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
    
    plot(nlcd[square], reset = FALSE,main="",axes=T)
    
}

#loop through and calculate % forest, water, matrix and NA values
# for each of the three scales

#convert raster to forest/matrix/water
landscape=nlcd[square]
data.frame(landcsape) %>% rename()

landscape %>% map(function(x) min(x))
yy = adrop(landscape)
landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
landscape_table=table(landscape_data[[1]])
landscape_table[landscape_table !=0]

#group the nlcd data into categories
forest_cats=c("Woody Wetlands","Deciduous Forest",'Mixed Forest',"Evergreen Forest")
water="Open Water"
wetland="Emergent Herbaceuous Wetland"
null_vals=c("Unclassified" ,"")
matrix_cats=unique(names(landscape_table))[!unique(names(landscape_table)) %in% c(forest_cats,water,wetland,null_vals)]

i=2

#first, 2400 m
lu_ls=list()
plan(multisession, workers = 6)
landscape_sizes=c(2400,1500,600)
lu_summs=landscape_sizes %>% purrr::map(function(landscape_size){
    1:length(grid)%>% future_map(possibly(function(i){
        square = st_sfc(st_buffer(grid[i], landscape_size/2,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
        landscape=nlcd[square]
        
        landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
        landscape_table=table(landscape_data[[1]])
        nonzero=landscape_table[landscape_table !=0]
        total_pixels =sum(landscape_table)#for checking that landscapes are the right size
        
        data.frame(
            
            prop_forest=sum(nonzero[names(nonzero) %in% forest_cats])/total_pixels,
            prop_matrix=sum(nonzero[names(nonzero) %in% matrix_cats])/total_pixels,
            prop_wetland=sum(nonzero[names(nonzero) ==wetland])/total_pixels,
            prop_water=sum(nonzero[names(nonzero) ==water])/total_pixels,
            prop_null=sum(nonzero[names(nonzero)  %in% null_vals])/total_pixels
        ) %>% 
            mutate(total_pixels=total_pixels,grid_index=i,size_m=landscape_size)
        
    },otherwise=NULL) )
})

    
summary(landscape_data)

