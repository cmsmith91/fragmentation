rm(list=ls())
library(tidyverse)
library(vroom)
library(furrr)
map=purrr::map

#set the working directory 
# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
if(getwd() !="~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation") setwd("~/Dropbox/Fall_2014/Research/Fragmentation/fragmentation")

#alt_path is the directory where the neutral-model community files are located
alt_path="~/Dropbox/Fall_2014/Research/Fragmentation/data/slurm-out"

length(list.files(paste0(alt_path,'/',"18065402")))

#load the forest-loss landscapes
hab_maps=readRDS('processed_data/hab_maps03dec2021.rds')
hab_maps_unlisted=hab_maps %>% unlist(recursive=F)

#load the data file with the land-use info (gives info about how much forest edge, area, etc)
edge=readRDS('processed_data/lu_info09march2022.rds')
focal_edge=edge %>% 
    map_dfr(function(df) df %>% filter(focal_landscape)) %>%
    mutate(landscape_id=1:n())

#get and save full path names of all the communities
path_names=paste0(alt_path,"/",list.files(alt_path)) %>% map(function(file_path){
    path_name_vec=paste0(file_path,'/',list.files(file_path))
    path_name_vec[!grepl('-env-',path_name_vec)]
    
}) %>% unlist

#extract the size of the communities and the density from the file names
path=path_names[1]
(size_xy=as.numeric(sub("_density.*","",sub(".*_sizexy","",path))))
(density=as.numeric(sub("_rep.*","",sub(".*_density","",path))))

#this code 
#code for converting from pixel_indices to x and y
ind_coord=1:(density*size_xy^2)
pixel_coord=floor((ind_coord-1)/density)+1
x=((pixel_coord-1)%%size_xy)+1
y=floor((pixel_coord-1)/size_xy)+1
coord_df=data.frame(ind_coord=ind_coord,pixel_coord=pixel_coord,x,y)


#set the number of communities per landscape
n_coms=100 #start with 5
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
tail(paths_split)
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
        
        #also record mean population sizes
        initial_n=median(table(combo$species[combo$species !=0]))
        final_n=median(table(loss_happens$species[loss_happens$species !=0]))
        
        data.frame(
            grid_index=f_map$grid_index[1],
            size_m=f_map$landscape_size[1],
            initial_rich=initial_rich,
            final_rich=final_rich,
            initial_n=initial_n,
            final_n=final_n,
            community=com_file_name
        ) %>% mutate(percent_loss=(initial_rich-final_rich)/initial_rich*100)
        
    
        })

    
    },.options=furrr_options(seed=T)) %>% bind_rows
future:::ClusterRegistry("stop")
# saveRDS(output,"processed_data/loss_simulation_n100_output.rds")
# output=readRDS("processed_data/loss_simulatin_n70_output.rds")

output_summ=output %>% group_by(grid_index,size_m) %>%
    summarize(mean_percent_loss=mean(percent_loss),
              median_percent_loss=median(percent_loss)) %>% 
    left_join(focal_edge) %>% 
    arrange(size_m) %>%
    mutate(size_ha=size_m^2*0.0001)



# pdf("figures/loss_sim_output_16march2022.pdf",width=13)
par(mfrow=c(1,5),pch=16,cex=1.6)
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

#get slope of edge/percent remaining
my_size=600
size_df=unique(output_summ$size_m) %>% map_dfr(function(my_size){
    df=output_summ %>% filter(size_m==my_size)
    unique(df$f_cat) %>% map_dfr(function(my_f_cat){
        new_df = df %>% filter(f_cat==my_f_cat) %>%mutate(percent_remaining=100-median_percent_loss)
        mod=with(new_df,lm(percent_remaining~edge_density))
        edge_relationship=coef(mod)
        conf_ints=confint(mod,"edge_density")
        lower=conf_ints[1]
        upper=conf_ints[2]
        
        tibble(size_m=my_size,size_ha=my_size^2/10000,f_cat=my_f_cat,
               edge_slope=edge_relationship[2],mod_intercept=edge_relationship[1],
               lower_ci=lower,upper_ci=upper)
    })
    
}) %>% mutate(size_ha_jit=jitter(size_ha))

my_size=600
my_cols=RColorBrewer::brewer.pal(8,'Paired')[c(7,2,4)]
my_cols_lines=my_cols
my_cols=adjustcolor(my_cols,alpha=.67)

#pdf("figures/loss_sim_output2-08march2022.pdf",width=18)
cex_lab=1.3
par(mfrow=c(1,5),pch=16,cex=1.4,cex.main=1.3,cex.lab=cex_lab,mar=c(4.1,2.5,2,.5))
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    size_plot=size_df %>% filter(size_m==my_size)
    
    
    
    if(my_size==600){
       par(oma=c(.5,1.5,1,1))
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],
                  main=paste0(size_ha[1],' ha'),ylim=c(0,100),xlim=c(10,200),
                  ylab='',xlab='')
        )
        mtext('% species remaining', side = 2, line = 0, outer = T,cex=1.9,adj=.6)
        legend('bottomright',cex=1,legend=legend_labs2,title = 'forest remaining',pch=16,col=my_cols[as.factor(legend_labs)],bty='n')
    }
    if(my_size==1740){
        
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],
                  main=paste0(size_ha[1],' ha'),ylim=c(0,100),xlim=c(10,200),
                  ylab='',xlab='edge density (m/ha)')
        )
        
    }
    if(!my_size %in% c(600,1740)){
        
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],
                  main=paste0(size_ha[1],' ha'),ylim=c(0,100),xlim=c(10,200),
                  ylab='',xlab='')
        )
        
        
    }
    #get lines for each value of forest area
    #start with 10% forest
    for(my_fcat in unique(plot_df$f_cat)){
        a=size_plot[size_plot$f_cat==my_fcat,]$mod_intercept
        b=size_plot[size_plot$f_cat==my_fcat,]$edge_slope
        
        
        #get min and max values for making sequence btw min and max
        fcat_points=plot_df %>% filter(f_cat==my_fcat) 
        min_x=min(fcat_points$edge_density); max_x=max(fcat_points$edge_density)
        seq_x=seq(min_x,max_x,by=(max_x-min_x)/1000)
        
        seq_y=a+b*seq_x
        line_col=my_cols_lines[which(legend_labs==my_fcat)]
        
        lines(seq_x,seq_y,col=line_col,lwd=2)
        
        }
    
}
# dev.off()


# pdf('figures/size_v_effectedge.pdf')
par(mfrow=c(1,1),pch=16,cex=1.6,mar=c(4.5,4,2,2))
with(size_df,plot(size_ha_jit,edge_slope,col=my_cols_lines[as.factor(f_cat)],ylim=c(0,.27),ylab="edge effect",xlab="landscape size (ha)"))
with(size_df,arrows(x0=size_ha_jit,y0=lower_ci,y1=upper_ci,angle=90,length=0,col=my_cols_lines[as.factor(f_cat)]))
legend('topright',cex=1,legend=legend_labs2,title = 'forest remaining',pch=16,col=my_cols_lines[as.factor(legend_labs)],bty='n')
#dev.off()

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
# pdf("figures/loss_sim_output3.pdf",width=13)
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
# dev.off()


# pdf("figures/loss_sim_output4.pdf",width=13)
par(mfrow=c(1,3),pch=16,cex=1.3,cex.lab=1.3)
unique(output_summ$f_cat)
for(my_propforest in legend_labs){
    plot_df=output_summ %>% filter(f_cat==my_propforest)
    plot_ls=plot_df %>% split(.$edge_cat) 
    plot_df_highedge=plot_df %>% filter(edge_cat=='high')
    plot_df_lowedge=plot_df %>% filter(edge_cat=='low')
    
    percentremaining_frag=100-plot_df_highedge$median_percent_loss
    percentremaining_cont=100-plot_df_lowedge$median_percent_loss
    species_diff=percentremaining_frag-percentremaining_cont
    edge_diff=plot_df_highedge$edge_density-plot_df_lowedge$edge_density
    
    effect_size=species_diff/edge_diff
    edge_eff_median=effect_size %>% split(plot_df_lowedge$size_ha) %>% map_dbl(median)
    with(plot_df_lowedge,
         plot(size_ha,effect_size,ylab='Effect size of fragmentation',
              xlab="Size (ha)",ylim=c(0,.45),main=paste0(my_propforest*100,'% forest remaining'))
    )
    points(x=unique(plot_df_lowedge$size_ha),y=edge_eff_median,col='red')
}
# dev.off()

# pdf("figures/forestvedge.pdf",width=12)
par(mfrow=c(1,3),pch=16,cex=1.3,cex.lab=1.3,mar=c(4,4,2,1))
for(i in 3:1){
    plot_df_a=edge[[i]] %>% arrange(focal_landscape) %>%
        mutate(prop_forest2=prop_forest^2)
    size_ha=(plot_df_a$size_m[1]^2)*0.0001
    
    
    fit_line=with(plot_df_a,lm(edge_density~prop_forest+prop_forest2))
    prediction_x=seq(min(plot_df_a$prop_forest),max(plot_df_a$prop_forest),by=(max(plot_df_a$prop_forest)-min(plot_df_a$prop_forest))/1000)
    prediction_x2=prediction_x^2
    a=coef(fit_line)[[1]]
    b=coef(fit_line)[[2]]*prediction_x
    c=coef(fit_line)[[3]]*prediction_x2
    
    x_lab=ifelse(i %in% c(3,2),'% forest',"") 
    prediction_y=a+b+c
    with(plot_df_a,
         plot(prop_forest*100,edge_density,
              col=as.factor(focal_landscape),
              main=paste0(size_ha," ha"),ylim=c(0,200),
              xlab=x_lab,ylab='edge density (m/ha)'
              ))
    lines(x=prediction_x*100,y=prediction_y,col=2,type='l',lwd=3)
}

# dev.off()

