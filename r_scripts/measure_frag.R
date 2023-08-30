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
library(landscapemetrics)
library(factoextra)
library(cluster)
library(NbClust)
select=dplyr::select
map=maps::map

# crop nlcd data to just be in the northeast US
# load nlcd data
# img_path = "~/Dropbox/Fall_2014/postdoc/inat project/data_current/NLCD_2016_LandCover/NLCD_2016_Land_Cover_L48_20190424.img"
img_path = '/Volumes/Seagate/nlcd/NLCD_2016_Land_Cover_L48_20190424.img'

nlcd=read_stars(img_path)


#get shapefile of northeastern US
northeast_us=tolower(c("Maine", "Vermont", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New York", "Pennsylvania", 
                       "New Jersey", "Delaware",  "Maryland"))

northeast_sf= st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
    filter(ID %in% northeast_us)
northeast_sf_crs=st_transform(northeast_sf,crs=st_crs(nlcd))

# To get the sum of each land-use in the northeasern USA, it is much faster to use QGIS than R
# to use QGIS:
# load into QGIS the shapefile of the northeastern USA and the raster NLCD data 
# create a buffer of distance 0 around the northeast shapefile
# and cut the raster data using the northeast shapefile:
# # # use Processing toolbox - GDAL - Raster Extraction - clip raster by mask layer
# then do raster layer unique values report

# code for downloading northeast shapefile
# st_write(northeast_sf_crs, "processed_data/northeast_shapefile.shp")

nlcd_table = read_csv('processed_data/nlcd unique values.csv')


#crop the nlcd data to be just states in the northeastern US
nlcd_cropped=nlcd %>% st_crop(northeast_sf_crs)

#make nlcd key for matrix habitat (in 'everything else' category)
nlcd_key = data.frame(category = c('developed','developed','developed','developed','barren','shrubland','shrubland','herbaceous','herbaceous','herbaceous','herbaceous','agriculture','agriculture'),
                      value = c(21,22,23,24,31,51,52,71,72,73,74,81,82))

#filter nlcd_table to just be matrix (ie everything in nlcd_key)
(nlcd_matrix = nlcd_key %>%
    left_join(nlcd_table) %>%
    filter(!is.na(count)) %>%
    group_by(category) %>%
    summarize(total_count=sum(count)) %>%
    mutate(percentage = 100*total_count/sum(total_count)))

#group shrubland, barren and herbaceous together
nlcd_matrix %>%
    filter(category  %in% c('barren',"shrubland", "herbaceous")) %>%
    summarize(sum(percentage))

# plot(nlcd_cropped)
# check = st_as_stars(nlcd_cropped==1)
# table(check)

#now make a grid over northeast us shapefile
#transform to crs in meters: utm zone are in meters - most of northeast is in 18n
#espg=26918
northeast_sf_crs_meters=st_transform(northeast_sf,crs=26918)
# ggplot()+geom_sf(data=northeast_sf_crs_meters)

#make a grid 
big_grid=st_make_grid(northeast_sf_crs_meters,5000,what = "centers") 
grid=big_grid[northeast_sf_crs_meters] %>% st_transform(crs=st_crs(nlcd))

#plot to make sure it looks ok
#just plot nj for now...
nj=northeast_sf_crs %>% 
    filter(ID=='new jersey')
grid_nj=grid[nj]
ggplot(data=northeast_sf_crs %>% filter(ID=='new jersey')) + geom_sf(fill='lightblue')+
    geom_sf(data=grid_nj,size=.06)

##
#northeast_sf_crs
## try looping through states, converting to stars object and getting amount of total habitat
## and the amount of matrix habitat
# northeast_data = st_as_stars(nlcd[nj])
# table(northeast_data)

#loop through each grid point and make a square buffer around it
grid[1]
square = st_sfc(st_buffer(grid[1], 1200,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
plot(nlcd[square], reset = FALSE,main="",axes=T)

# for(i in 1:6){
#     square = st_sfc(st_buffer(grid[i], 1200,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
#     
#     plot(nlcd[square], reset = FALSE,main="",axes=T)
#     
# }

#loop through and calculate % forest, water, matrix and NA values
# for each of the three scales

#convert raster to forest/matrix/water
landscape=nlcd[square]

landscape %>% purrr::map(function(x) min(x))
yy = adrop(landscape)
landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
landscape_table=table(landscape_data[[1]])
landscape_table[landscape_table !=0]

#group the nlcd data into categories
forest_cats=c("Woody Wetlands","Deciduous Forest",'Mixed Forest',"Evergreen Forest")
water="Open Water"
wetland="Emergent Herbaceuous Wetlands"
null_vals=c("Unclassified" ,"")
matrix_cats=unique(names(landscape_table))[!unique(names(landscape_table)) %in% c(forest_cats,water,wetland,null_vals)]


##
i=2

#first, 2400 m
lu_ls=list()
plan(multisession, workers = 6)
landscape_sizes=c(2400,1500,600)


all_the_squares=(20:80)^2
seq(400,6400,by=(6400-400)/4)

#get new lanscape sizes
vec=c()
for(numb in c(1900 ,3400, 4900)){
    i=which(c(1900 ,3400, 4900)==numb)
    new=all_the_squares[abs(all_the_squares-numb)==min(abs(all_the_squares-numb))]
    vec[i]=new
}
landscape_sizes=c(20, 44, 58, 70, 80)*30
plot(landscape)

#for testing the map functions:
landscape_size = 1320


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
        
    },otherwise=NULL) ,.options=furrr_options(seed=T))
})


    
i=2

lu_summs[[i]]

#next measure edge
landscape_size=landscape_sizes[1]
the_edges=landscape_sizes %>% purrr::map(function(landscape_size){
    1:length(grid)%>% future_map(possibly(function(i){
    square = st_sfc(st_buffer(grid[i], landscape_size/2,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
    landscape=nlcd[square]
    landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
    
    #group all the forest categories
    landscape_data[[1]][landscape_data[[1]] %in% forest_cats & !is.na(landscape_data[[1]])] <- "Deciduous Forest" # group forest
    
    #group water and wetland
    landscape_data[[1]][landscape_data[[1]] == wetland & !is.na(landscape_data[[1]])] <- "Open Water"
    
    #group everything else
    landscape_data[[1]][landscape_data[[1]] %in% matrix_cats & !is.na(landscape_data[[1]])] <- "Cultivated Crops" 
    
    
    lsm_c_ed(landscape_data)

    },otherwise=NULL),.options=furrr_options(seed=T))
})
the_edges[[2]][10]
#the tricky thing is, is that we don't care about forest edge with water
#to get just forest edge subtract the water edge from the sum  of all edges
# all edges = sum(value)/2

#categories: 12=ice snow;42=evergreen forest;83 does not exist; #seems like it added 1 to the true categories?
i=1
square = st_sfc(st_buffer(grid[i], landscape_size/2,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
landscape=nlcd[square]
plot(landscape)
df=the_edges[[1]][[1]]
ls
fm_edge=the_edges %>% purrr::map(function(ls){
    ls %>% purrr::map(possibly(function(df){
        if(12 %in% df$class){
            fm_edge=sum(df$value)/2-df[df$class==12,]$value
            }
        else{
            sum(df$value)/2
        }
        
        
    },otherwise=NULL)) 
})



par(mfrow=c(1,3))
for(i in 1:3){
    test1=fm_edge[[i]]
    a=lu_summs[[i]] %>% bind_rows %>% mutate(edge_density=test1)
    keep_me=which(a$prop_water+a$prop_wetland <.05)
    with(a[keep_me,], plot(prop_forest,edge_density,ylim=c(0,230)))
    
}
#get rid of landscapes overlapping the edge of the raster data [or where there is null raster data]
#or where there is greater than 5% water/wetland
remove_me=lu_summs %>% purrr::map(function(ls){
    df=ls %>% bind_rows %>% filter(prop_null !=0 | (prop_wetland+prop_water)>0.05) 
    df$grid_index
}) %>% unlist %>% unique()

lu_df=1:length(lu_summs) %>% purrr::map(function(i){
    test1=fm_edge[[i]] %>% unlist
    a=lu_summs[[i]] %>% bind_rows %>% mutate(edge_density=test1)
    keep_me=which(!a$grid_index %in% remove_me)
    
    a[keep_me,]
})
lu_df2=lu_df %>% purrr::map(function(df){
    df$f_cat=NA
    df$edge_cat=NA
    df$focal_landscape=F
    for(f_cat in c(.1,.3,.5) ){
        f_range=c(f_cat-.025,f_cat+0.025)
        
        #filter to only be landsacpes with f_cat amount of forest plus or minus 2.5 percent
        df_ordered=df %>% 
            filter(prop_forest > f_range[1] & prop_forest <f_range[2]) %>%
            arrange(edge_density)
        
        #get the indices of the 15 landcapes with the least edge
        # and the 15 with the most
        low_edge_is=df_ordered[1:15,]$grid_index
        high_edge_is=df_ordered[(nrow(df_ordered)-14):nrow(df_ordered),]$grid_index
        
        #modify df: label as a focal landscape, and give which forest size category its in
        df[df$grid_index %in% low_edge_is | df$grid_index %in% high_edge_is,]$f_cat<-f_cat
        df[df$grid_index %in% low_edge_is | df$grid_index %in% high_edge_is,]$focal_landscape<-T
        
        #modify df: label as high or 
        df[df$grid_index %in% low_edge_is,]$edge_cat<-'low'
        df[df$grid_index %in% high_edge_is,]$edge_cat<-'high'
        
    }
    df
    
})

# plot to make sure everything looks ok
for(df in lu_df2){
    with(df,plot(prop_forest,edge_density,col=ifelse(focal_landscape,'red','black')))
}

#get the maps of each of the focal_landscapes 
plan(multisession,workers=5)
hab_maps=1:length(lu_df2) %>% purrr::map(function(i){
    df=lu_df2[[i]]
    landscape_size=landscape_sizes[[i]]
    get_these=df[df$focal_landscape,]$grid_index
    
    get_these %>% future_map(possibly(function(j){
        square = st_sfc(st_buffer(grid[j], landscape_size/2,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
        landscape=nlcd[square]
        
        landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
        
        #group all the forest categories
        landscape_data[[1]][landscape_data[[1]] %in% forest_cats & !is.na(landscape_data[[1]])] <- "Deciduous Forest" # group forest
        
        #group water and wetland
        landscape_data[[1]][landscape_data[[1]] == wetland & !is.na(landscape_data[[1]])] <- "Open Water"
        
        #group everything else
        landscape_data[[1]][landscape_data[[1]] %in% matrix_cats & !is.na(landscape_data[[1]])] <- "Cultivated Crops" 
        
        #now convert from a matrix to a df with x-y coordinates
        long_df=1:ncol(landscape_data[[1]]) %>% purrr::map_dfr(function(x){
            1:nrow(landscape_data[[1]]) %>% purrr::map_dfr(function(y){
                data.frame(x=x,y=y,category=as.character(landscape_data[[1]][y,x]))
            })
        }) %>% filter(!is.na(category))
        
        #change labeels to be shorter and more general
        if('Cultivated Crops' %in% long_df$category) long_df[long_df$category=='Cultivated Crops',]$category='matrix'
        if('Deciduous Forest' %in% long_df$category)long_df[long_df$category=='Deciduous Forest',]$category='forest'
        if('Open Water' %in% long_df$category) long_df[long_df$category=='Open Water',]$category<-'water'
        
        long_df %>%
            mutate(grid_index=j,landscape_size=landscape_size)
        
        },otherwise=NULL),.options=furrr_options(seed=T))
    
    })
future:::ClusterRegistry("stop")

        
hab_maps %>% purrr::map(function(ls){
    with(ls[[1]],plot(x,y,col=as.factor(category)))

    })  

# saveRDS(hab_maps,'processed_data/hab_maps03dec2021.rds')
# saveRDS(lu_df2,'processed_data/lu_info03dec2021.rds')

#old code

head(lu_summs[[1]] %>% bind_rows)
length(test1 %>% unlist)

nrow(lu_summs[[1]] %>% bind_rows)
test1 %>% purrr::map(function(df) is_empty(df))

test1[[989]]

the_edges[[1]][[989]]

square = st_sfc(st_buffer(grid[i], landscape_size/2,endCapStyle="SQUARE",joinStyle="MITRE"), crs = st_crs(nlcd))
landscape=nlcd[square]
landscape_data=st_as_stars(landscape) # convert from stars_proxy object to stars object with all the data
ls_tb=table(landscape_data[[1]])
sum(ls_tb[ls_tb !=0])
landscape_data %>% pull(1)
head(matrix(landscape_data[[1]]))#[,1:3] #%in% forest_cats
attr(landscape_data[[1]],'colors') <- NULL
landscape_data[[1]]<-landscape_data[[1]] %in% forest_cats & !is.na(landscape_data[[1]])
plot(landscape_data)
is.na(landscape_data[[1]])

#convert to forest/non-forest:
landscape
f_landscape
forest_nonforest=landscape_data[[1]] %in% forest_cats
forest_nonforest=landscape_data %>% st_apply(c('x','y'),function(x) x %in% forest_cats)
plot(forest_nonforest,breaks='equal')
lsm_c_ed

