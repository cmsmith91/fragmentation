rm(list=ls())
library(tidyverse)
library(vroom)
library(furrr)
map=purrr::map

# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
if(getwd() !="~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation") setwd("~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation")
alt_path="~/Dropbox/Fall_2014/Research/Fragmentation/data/slurm-out"

hab_maps=readRDS('processed_data/hab_maps18nov2021.rds')
hab_maps_unlisted=hab_maps %>% unlist(recursive=F)
edge=readRDS('processed_data/lu_info18nov2021.rds')
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
n_coms=100#start with 5
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
            size_m=f_map$landscape_size[1],
            initial_rich=initial_rich,
            final_rich=final_rich,
            community=com_file_name
        ) %>% mutate(percent_loss=(initial_rich-final_rich)/initial_rich*100)
        
    
        })

    
    },.options=furrr_options(seed=T)) %>% bind_rows
future:::ClusterRegistry("stop")
# saveRDS(output,"processed_data/loss_simulatin_output.rds")
output=readRDS("processed_data/loss_simulatin_output.rds")

output_summ=output %>% group_by(grid_index,size_m) %>%
    summarize(mean_percent_loss=mean(percent_loss),
              median_percent_loss=median(percent_loss)) %>% 
    left_join(focal_edge) %>% 
    arrange(size_m) %>%
    mutate(size_ha=size_m^2*0.0001)



# pdf("figures/loss_sim_output1.pdf",width=13)
par(mfrow=c(1,3),pch=16,cex=1.6)
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    with(plot_df,
         plot(prop_forest,median_percent_loss,col=as.factor(edge_cat),
              main=paste0('size = ', size_ha[1],' ha'),ylim=c(0,60))
         )
}
# dev.off()
plot_df=output_summ %>% filter(size_m==600)
legend_labs=unique(plot_df$f_cat)[order(unique(plot_df$f_cat))]
legend_labs2=paste0(legend_labs*100,'%')

pdf("figures/loss_sim_output2.pdf",width=13)
par(mfrow=c(1,3),pch=16,cex=1.6,cex.main=1.3)
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    with(plot_df,
         plot(edge_density,100-median_percent_loss,col=as.factor(f_cat),
              main=paste0(size_ha[1],' ha'),ylim=c(0,100),xlim=c(10,200),
              ylab='% of species remaining',xlab='edge density (m/ha)')
    )
    if(my_size==600){
        legend('bottomright',cex=.8,legend=legend_labs2,title = 'forest remaining',pch=16,col=as.factor(legend_labs),bty='n')
    }
}
dev.off()


# dev.off()

#now try plotting edge density v species loss with forest amount being different colors
# pdf("figures/loss_sim_output1.pdf",width=13)
par(mfrow=c(1,3),pch=16,cex=1.3)
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    with(plot_df,
         plot(prop_forest,median_percent_loss,col=as.factor(edge_cat),
              main=paste0('size = ', size_m[1]),ylim=c(0,60))
    )
}
# dev.off()

#now try plotting edge density v species loss with forest amount being different colors
pdf("figures/loss_sim_output3.pdf",width=13)
par(mfrow=c(1,3),pch=16,cex=1.3,cex.lab=1.3)
unique(output_summ$f_cat)
for(my_propforest in legend_labs){
    plot_df=output_summ %>% filter(f_cat==my_propforest)
    plot_ls=plot_df %>% split(.$edge_cat) 
    plot_df_highedge=plot_df %>% filter(edge_cat=='high')
    plot_df_lowedge=plot_df %>% filter(edge_cat=='low')
    
    percentremaining_frag=100-plot_df_highedge$median_percent_loss
    percentremaining_cont=100-plot_df_lowedge$median_percent_loss
    percent_diff=100*(percentremaining_frag-percentremaining_cont)/percentremaining_frag
    

    with(plot_df_lowedge,
         plot(size_ha,percent_diff,ylab='Species remaining: fragmented-continuous',
              xlab="Size (ha)",main=paste0(my_propforest*100,'% forest remaining'))
    )
}
dev.off()

