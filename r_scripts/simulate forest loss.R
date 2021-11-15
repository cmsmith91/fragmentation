rm(list=ls())
library(tidyverse)
library(vroom)
library(furrr)

# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
if(getwd() !="~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation") setwd("~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation")
alt_path="~/Dropbox/Fall_2014/Research/Fragmentation/data/slurm-out"

hab_maps=readRDS('processed_data/hab_maps10nov2021.rds')
hab_maps_unlisted=hab_maps %>% unlist(recursive=F)
edge=readRDS('processed_data/lu_info10nov2021.rds')
focal_edge=edge %>% 
    map_dfr(function(df) df %>% filter(focal_landscape)) %>%
    mutate(landscape_id=1:n())

#get path names of all the communities
path_names=paste0(alt_path,"/",list.files(alt_path)) %>% map(function(file_path){
    path_name_vec=paste0(file_path,'/',list.files(file_path))
    path_name_vec[!grepl('-env-',path_name_vec)]
    
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


#set the number of coms per landscape
n_coms=50#start with 5
n_landscapes=sum(hab_maps %>% map_int(length))
replicates_needed=n_landscapes*n_coms #number of landscapes times n_coms


# # get row indices for the three different landscape sizes
# n_sizes=length(hab_maps)
# size_index=2
# 1:n_sizes %>% map(function(size_index){
#     total_pixels=hab_maps[[size_index]][[1]]$landscape_size[1]
#     
# })

#assign n_coms path names to each 
paths_split=path_names[1:replicates_needed] %>% split(rep(1:n_landscapes,each=n_coms))
i=1;
com_file_name=paths_split[[1]][1]

plan(multisession,workers=6)
output=1:n_landscapes %>% future_map(function(i){
    com_paths=paths_split[[i]]
    f_map=hab_maps_unlisted[[i]]
    ls_size=nrow(f_map)
    ls_size_xy=sqrt(ls_size)
    
    com_paths %>% map_dfr(function(com_file_name){
        com=vroom(com_file_name,delim=" ",col_names=c('species'))
        
        #cut the community out to be the right size
        #first figure out waht coordinates (starting fromthe cetner) are in a community of that size
        midway_pt=(size_xy+1)/2
        upper=floor(midway_pt + ls_size_xy/2)
        lower = floor(midway_pt -ls_size_xy/2)+1
        
        #cut out the community
        row_indices=with(coord_df,which(x<=upper &  x >= lower & y<=upper & y>=lower ))
        sp_map=com[row_indices,] %>% bind_cols(coord_df[row_indices,])
        
        min_xy=min(sp_map$x) # lower the x and y-values the x and
        sp_map = sp_map%>% mutate(xy=paste0(x,y),
                                  x_lowered=x-min_xy,
                                  y_lowered=y-min_xy) 
        
        #now join the community with its associated forest map
        combo = left_join(sp_map,f_map %>% rename(
            x_lowered = x,y_lowered =y
        ),by = c("x_lowered", "y_lowered")) 
        
        #exclude species in water
        if('water' %in% f_map$category){
            combo[which(combo$category=='water'),]$species <- 0 #get rid of any individuals in the water
        }
        
        loss_happens=combo
        
        #cookie cutter
        loss_happens[which(loss_happens$category=='matrix'),]$species <- 0 #individuals in the matrix go extinct
        
        #record initial richness and final richness
        initial_rich=n_distinct(combo$species[combo$species !=0])
        final_rich=n_distinct(loss_happens$species[loss_happens$species !=0])
        
        data.frame(
            grid_index=f_map$grid_index[1],
            landscape_size=f_map$landscape_size[1],
            initial_rich=initial_rich,
            final_rich=final_rich,
            community=com_file_name
        ) %>% mutate(percent_loss=(initial_rich-final_rich)/initial_rich*100)
        
    
        })

    
    },.options=furrr_options(seed=T)) %>% bind_rows
future:::ClusterRegistry("stop")

output_summ=output %>% group_by(grid_index,landscape_size) %>%
    summarize(mean_percent_loss=mean(percent_loss),
              median_percent_loss=median(percent_loss)) %>% 
    left_join(focal_edge) %>% 
    arrange(size_m)
par(mfrow=c(1,3))
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    with(plot_df,
         plot(prop_forest,median_percent_loss,col=as.factor(edge_cat),
              main=paste0('size = ', size_m[1]))
         )
}
