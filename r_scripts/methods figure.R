rm(list=ls())
pardefault <- par()

library(tidyverse)
library(vroom)
library(furrr)
library(RColorBrewer)
library(mapdata)
library(sf)
library(ggspatial)
map=maps::map

#cols for plotting
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_list2 = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_df2 = data.frame('species'=1:length(col_list2),col = col_list2,stringsAsFactors = F)
col_df1=data.frame('species'=1:length(col_list2),col = rev(col_list2),stringsAsFactors = F)

# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
if(getwd() !="~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation") setwd("~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation")
alt_path="~/Dropbox/Fall_2014/Research/Fragmentation/data/slurm-out"

hab_maps=readRDS('processed_data/hab_maps03dec2021.rds')
a=hab_maps %>% unlist(recursive=F)
df=a[[1]]
hab_maps_unlisted=hab_maps %>% unlist(recursive=F) %>%
    purrr::map(function(df){
        max_xy=sqrt(nrow(df))
        df$x=rep(1:max_xy,each=max_xy)
        df$y=rep(1:max_xy,max_xy)
        return(df)
    })
edge=readRDS('processed_data/lu_info03dec2021.rds')
focal_edge=edge %>% 
    map_dfr(function(df) df %>% filter(focal_landscape)) %>%
    mutate(landscape_id=1:n())


#get path names of all the communities
path_names=paste0(alt_path,"/",list.files(alt_path)) %>% purrr::map(function(file_path){
    path_name_vec=paste0(file_path,'/',list.files(file_path))
    path_name_vec[!grepl('-env-',path_name_vec) & !grepl('prop',path_name_vec) & grepl("steps100000",path_name_vec)]
    
}) %>% unlist

#get size_xy and density
path=path_names[1]
(size_xy=as.numeric(sub("_density.*","",sub(".*_sizexy","",path))))
(density=as.numeric(sub("_rep.*","",sub(".*_density","",path))))

#code for converting from pixel_indices to x and y
ind_coord=1:(density*size_xy^2)
pixel_coord=floor((ind_coord-1)/density)+1
x=((pixel_coord-1)%%size_xy)+1
y=floor((pixel_coord-1)/size_xy)+1
coord_df=data.frame(ind_coord=ind_coord,pixel_coord=pixel_coord,x,y)

#get one high edge landscape and one low from 2400
new_df=focal_edge %>% filter(size_m==2400,f_cat==.3)
high_index=new_df[new_df$edge_cat=='high',]$grid_index[7]
low_index=new_df[new_df$edge_cat=='low',]$grid_index[5]

#for ease of plotting get landscapes with no water
no_water=which(1:length(hab_maps_unlisted) %>% map_lgl(function(i) !'water' %in% hab_maps_unlisted[[i]]$category))
no_water_info=no_water %>% map_dfr(function(i){
    df=hab_maps_unlisted[[i]] 
    data.frame(grid_index=df$grid_index[1],landscape_size=df$landscape_size[1])
}) %>% left_join(focal_edge %>% select(grid_index,f_cat,edge_cat))
high_nowater=no_water_info %>% filter(landscape_size==2400,f_cat==0.5,edge_cat=='high')
low_nowater=no_water_info %>% filter(landscape_size==2400,f_cat==0.5,edge_cat=='low')
highedge_gridindex=high_nowater$grid_index[1]
lowedge_gridindex=low_nowater$grid_index[2]

df=hab_maps_unlisted[[1]]
highedge_index=which(hab_maps_unlisted %>% map_lgl(function(df){
    df$landscape_size[1]==2400 & df$grid_index[1] ==highedge_gridindex
}))
lowedge_index=which(hab_maps_unlisted %>% map_lgl(function(df){
    df$landscape_size[1]==2400 & df$grid_index[1] ==lowedge_gridindex
}))
highedge_ls=hab_maps_unlisted[[highedge_index]] %>% mutate(color=ifelse(category=="forest",'darkgreen','white'))
if(nrow(highedge_ls[highedge_ls$category=='water',]) !=0) highedge_ls[highedge_ls$category=='water',]$color <-'blue'


lowedge_ls=hab_maps_unlisted[[lowedge_index]] %>% mutate(color=ifelse(category=="forest",'darkgreen','white'))
if(nrow(lowedge_ls[lowedge_ls$category=='water',]) !=0) lowedge_ls[lowedge_ls$category=='water',]$color <-'blue'

par(mfrow=c(1,2),pty="s",cex=.55,yaxt="n",xaxt='n')
with(lowedge_ls,plot(x,y,col=color,pch=15,xlab='',ylab=''))
with(highedge_ls,plot(x,y,col=color,pch=15,xlab='',ylab=''))


#now upload 1 community 
#get path names of all the communities
com1=vroom(path_names[3],delim=" ",col_names=c('species'))

#
n_distinct(com1$species)

#
sp_map=com1
unique_sp=sp_map %>% group_by(species) %>% summarize(n=n()) %>% arrange(desc(n))
unique_cols=unique_sp %>%
     mutate(species2 = 1:nrow(unique_sp)) %>% left_join(col_df1 %>% rename(species2=species))
plot_df1=sp_map %>% left_join(unique_cols %>% select(species,col)) %>%
    mutate(ind_coord=1:n()) %>% left_join(coord_df)

sp_map2=com1
unique_sp=sp_map2 %>% group_by(species) %>% summarize(n=n()) %>% arrange(desc(n))
unique_cols=unique_sp %>%
    mutate(species2 = 1:nrow(unique_sp)) %>% left_join(col_df1 %>% rename(species2=species))
plot_df2=sp_map2 %>% left_join(unique_cols %>% select(species,col)) %>%
    mutate(ind_coord=1:n()) %>% left_join(coord_df)

par(mfrow=c(1,2),cex=.4,pch=15,pty="s",xaxt='n',yaxt='n')
with(plot_df1,plot(jitter(x),jitter(y),col=col,xlab='',ylab=''))
with(plot_df2,plot(jitter(x),jitter(y),col=col,xlab='',ylab=''))

#now simulate forest lost and plot the new communities

loss_plot1=plot_df1 %>% left_join(lowedge_ls)
loss_plot1[loss_plot1$category!='forest',]$col = 'white'
loss_plot1[loss_plot1$category!='forest',]$species = 0


loss_plot2=plot_df2 %>% left_join(highedge_ls)
loss_plot2[loss_plot2$category!='forest',]$col = 'white'
loss_plot2[loss_plot2$category!='forest',]$species = 0


#number of species before and after
lowedge_before=n_distinct(plot_df1$species)
highedge_before=n_distinct(plot_df2$species)
lowedge_after=n_distinct(loss_plot1[loss_plot1$species !=0,]$species)
highedge_after=n_distinct(loss_plot2[loss_plot2$species !=0,]$species)


par(mfrow(1,2))
with(loss_plot1,plot(jitter(x),jitter(y),col=col,xlab='',ylab=''))
with(loss_plot2,plot(jitter(x),jitter(y),col=col,xlab='',ylab=''))

#pdf('figures/methods_fig-24july2022.pdf',width=12,height=8)
par(mfcol=c(2,3),mar=c(.2,.2,3.5,.2),pty="s",pch=15,cex=.59,yaxt="n",xaxt='n',cex.main=4)
#forest landscapes
with(lowedge_ls,plot(x,y,col=color,pch=15,xlab='',ylab='',main='low edge'))
with(highedge_ls,plot(x,y,col=color,pch=15,xlab='',ylab='',main='high edge'))

#communities before forest loss
with(plot_df1,plot(jitter(x),jitter(y),col=col,xlab='',ylab='',main=paste0(lowedge_before,' species')))
with(plot_df2,plot(jitter(x),jitter(y),col=col,xlab='',ylab='',main=paste0(highedge_before,' species')))

#coms after forest loss
with(loss_plot1,plot(jitter(x),jitter(y),col=col,xlab='',ylab='',main=paste0(lowedge_after,' species')))
with(loss_plot2,plot(jitter(x),jitter(y),col=col,xlab='',ylab='',main=paste0(highedge_after,' species')))
#dev.off()


par(pardefault)

#####
##now add map
#get shapefile of northeastern US
northeast_us=tolower(c("Maine", "Vermont", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New York", "Pennsylvania", 
                       "New Jersey", "Delaware",  "Maryland"))

northeast_sf= st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
    filter(ID %in% northeast_us)
northeast_sf_crs_meters=st_transform(northeast_sf,crs=26918)
#make a grid 

big_grid=st_make_grid(northeast_sf_crs_meters,5000,what = "centers") 
big_grid2=st_make_grid(northeast_sf_crs_meters,10000,what = "centers") 

grid=big_grid[northeast_sf_crs_meters] 
grid2=st_transform(grid,crs = 4326)




# load the data for the map of the northeastern US
test <- map("world2Hires",region=c("USA","Canada"),fill=T,plot=F)
test$x <- test$x-360
test2 = st_as_sf(test, crs = 4326)
#get state lines
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

par(mfrow=c(1,1))
# RColorBrewer::display.brewer.all(colorblindFriendly = T)
cols_for_map=adjustcolor(RColorBrewer::brewer.pal(5,'Dark2'),.5)
adjustcolor(cols_for_map,.5)
########
region_map=ggplot(data=test2) +
    geom_sf(fill='antiquewhite') +
    geom_sf(data = states, fill = NA)+
    coord_sf(crs = st_crs(4326),
             xlim = c(-82,-66.7), 
             ylim = c(36.5,48.2), expand = FALSE)+
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Longitude") + 
    ylab("Latitude")

#format grid data to add to map
lat_lon2=data.frame(st_coordinates(grid2)) %>% rename(longitude=X,latitude=Y)%>% 
    mutate(state=maps::map.where("state",longitude,latitude),
           country= maps::map.where("world",longitude,latitude),
           grid_index=1:nrow(.))
lat_lon_focal2400=focal_edge %>% filter(size_m==2400) %>%
    select(f_cat,edge_cat,grid_index) %>%
    left_join(lat_lon2) %>%
    mutate(color=ifelse(edge_cat=='low',cols_for_map[5],cols_for_map[2]), #low edge is green and high edge is orange/red
           shape=ifelse(grid_index %in% c(highedge_gridindex,lowedge_gridindex),17,16),
           size=ifelse(grid_index %in% c(highedge_gridindex,lowedge_gridindex),4,2.5))
# 
grid_sites_over_region2400=region_map+    
        geom_point(data = lat_lon_focal2400 , aes(x=longitude,y=latitude,
                   fill=color),shape=21,
                   size=2)+
    theme(text = element_text(size=35),
          panel.background = element_rect(fill = "aliceblue"),
          axis.text.x = element_text(size = 17),
          axis.text.y = element_text(size = 17),
          
          legend.position=c(.84,.1),
          
          # Change background color underneath legend keys
          legend.key = element_rect("transparent", color=NA),
          legend.background = element_rect('transparent', color=NA, size=1),
          
          # Change font size of legend text, title 
          legend.text = element_text(size=16),
          legend.title = element_text(size=12),
          #remove legend margin
          legend.margin = margin(0, 0, 0, 0, "cm"))+
    #modify manually fill of points in map, plus the color, fill, shape of points in legend/guide
    scale_fill_manual(name = "", labels = c('low edge',"high edge"),
                      values=c(cols_for_map[5],cols_for_map[2]),
                      guide = guide_legend(override.aes = list(size = c(3.2,3.2),
                                                               fill=c(cols_for_map[5],cols_for_map[2]),
                                                               shape=21
                                                               
                      )))
          
# pdf('figures/map_fig1.pdf',height=10)
grid_sites_over_region2400
# dev.off()




##old 
lat_lon_focal1500=focal_edge %>% filter(size_m==1500) %>%
    select(f_cat,edge_cat,grid_index) %>%
    left_join(lat_lon2) %>%
    mutate(color=ifelse(edge_cat=='low',cols_for_map[3],cols_for_map[5]))
lat_lon_focal600=focal_edge %>% filter(size_m==600) %>%
    select(f_cat,edge_cat,grid_index) %>%
    left_join(lat_lon2) %>%
    mutate(color=ifelse(edge_cat=='low',cols_for_map[1],cols_for_map[2]))

(grid_sites_over_region2400=region_map+    
        geom_point(data = lat_lon_focal2400, aes(x=longitude,y=latitude),size=.9,col=lat_lon_focal2400$color))
(grid_sites_over_region1500=region_map+    
        geom_point(data = lat_lon_focal1500, aes(x=longitude,y=latitude),size=.9,col=lat_lon_focal1500$color))
(grid_sites_over_region600=region_map+    
        geom_point(data = lat_lon_focal600, aes(x=longitude,y=latitude),size=.9,col=lat_lon_focal600$color))

        
        
      
