rm(list=ls())
setwd("~/Dropbox/Fall_2014/Research/Fragmentation/data/neutralmod_output_22may2021")
setwd("~/Dropbox/Fall_2014/Research/Fragmentation/data/neutralmod_output_25may2021")

library(truncnorm)
library(tidyverse)
library(vroom)
library(vegan)
library(furrr)
map=purrr::map; select=dplyr::select

## extract parameter info from the file name
files=list.files()
global=as.numeric(gsub("_size.*","",gsub(".*_global",'',files)))
rep=as.numeric(gsub('.*rep','',gsub('timesteps.*','',files)))
spec=as.numeric(gsub("_global.*","",gsub('.*spec','',files)))
timestep=as.numeric(gsub('\\.txt.*','',gsub('.*timesteps','',files)))
density=as.numeric(gsub('_rep.*','',gsub('.*density','',files)))
size_pix=as.numeric(gsub('_den.*','',gsub('.*sizexy','',files)))
size_ha=(size_pix*30)^2/10000
size_m=size_pix*30
df=data.frame(spec=spec,global=global,size_pix=size_pix,size_ha=size_ha,size_m=size_m,
              density=density,rep=rep,timestep=timestep,filename=files)
param_df=df
check=param_df %>% group_by(spec,global) %>% summarize(n=n_distinct(rep))
check %>% filter(global !=.0001 | spec==0.000009) %>% 
    mutate(to_do=5-n) %>% filter(to_do != 0)
check %>% filter(global !=.0001 | spec==0.000009) %>%
    arrange(global)

## load the bee data
# load list of forest-associated bees
fbee=read.csv("~/Dropbox/Fall_2014/Research/Fragmentation/data/forest_bees.csv") %>% 
    select(genus_species)

# load site data
# only include sites sampled in both years 
# (year==2018 bc all sites sampled in 2018 were sampled in 2017 but not vice versa) to be all sites; filter to just sites sampled in both 2017 and 2018:
sites=read.csv("~/Dropbox/Fall_2014/Research/Fragmentation/data/ForestBeeSiteInfo_final.csv") %>%
    filter(year==2018) %>% select(-X)
#load specimen data 
fbees_allsites=read.csv("~/Dropbox/Fall_2014/Research/Fragmentation/data/forestbee1718_spec.csv") %>% 
    mutate(genus_species=paste(genus,species,sep="_"),
           latitude=round(latitude,4),longitude=round(longitude,4)) %>%
    filter(site %in% sites$site & genus_species %in% fbee$genus_species)
# combine site and specimen data, calculate, abundance and diversity (pooled for both 2017 and 2018)
local_all = fbees_allsites %>% 
    group_by(site, genus_species) %>% summarize(n=n()) %>%
    summarize(rich = length(unique(genus_species)),abund=sum(n),
              shan=exp(diversity(n,'shannon')),simp=diversity(n,'invsimpson')) %>%
    left_join(sites)

#Fit a GLM with richness as the response and forest cover as the predictor
with(local_all,plot(forest500,rich))
rich_mod=glm(rich~forest500,family='poisson',data=local_all)
rich_mod_lm=lm(rich~forest500,data=local_all)

# Use the model to predict bee richness in a landscape with 100% forest cover 
(predicted_richness=exp(as.numeric(coef(rich_mod)[1]+coef(rich_mod)[2]*1)))
(predicted_richness2=as.numeric(coef(rich_mod_lm)[1]+coef(rich_mod_lm)[2]*1))

#get the ci around the model prediction
# https://fromthebottomoftheheap.net/2017/05/01/glm-prediction-intervals-i/
# a Wald confidence interval on the fitted function based on the standard 
# errors of the estimates of the model coefficients
p <- predict(rich_mod, newdata = data.frame(forest500 = 1), se.fit=TRUE)
predicted_richness=p$fit
(upper_rich=exp(p$fit+p$se.fit*2))
(lower_rich=exp(p$fit-p$se.fit*2))

p <- predict(rich_mod, newdata = data.frame(forest500 = 1), se.fit=TRUE)
exp(p$fit+p$se.fit*2)

# fit_resp  = ilink(fit_link)
# right_upr = ilink(fit_link + (2 * se_link)),
# right_lwr = ilink(fit_link - (2 * se_link)))
sigma(rich_mod_lm)
exp(sigma(rich_mod))

#Fit a linear model with abundance as the response variable and 
#forest cover as the predictor variable. 
with(local_all,plot(forest500,abund)) #plot
with(local_all,plot(forest500,log(abund))) #plot

mod_abund=lm(abund~forest500,data=local_all)
# Use the model to predict bee abundance in a landscape with 100% forest cover 
predict_abund=as.numeric(coef(mod_abund)[1]+coef(mod_abund)[2]*1)
sigma_abund=sigma(mod_abund)

coef(summary(mod_abund))


# load simulated data and look at how div changes with each timestep
#calculate the diversity of each community
coms=list.files() %>% map(function(file_name) vroom(file_name,delim="\"",col_names = F))
#Error: Unknown TZ UTC

plan(multisession,workers=6)
div=coms %>% future_map(function(df){
    tab=table(df$X1)
    data.frame(
        richness=length(tab),
               shannon=exp(diversity(tab)),
               simpson=diversity(tab,'invsimpson'))
    
})
focal_spec=as.numeric(names(table(spec)[order(table(spec),decreasing=T)])[1])
focal_global=as.numeric(names(table(global)[order(table(global),decreasing=T)])[1])

div_df=div %>% bind_rows %>% bind_cols(df) %>%
    arrange(desc(spec),desc(global)) %>%
    mutate(vary_spec=global==focal_global,vary_global=spec==focal_spec)

par(mfrow=c(3,3),mar=c(4,5,5,1))
for(i in unique(div_df$spec)){
    summ=div_df %>% filter(vary_spec & spec==i) %>% group_by(timestep) %>% 
        summarize(rich=mean(richness),n=n())
    with(summ %>% filter(timestep !=1),
         plot(rich~timestep,main=paste0('spec = ',i,
                                        '\n n = ',n[1])))
}
par(mfrow=c(3,3),mar=c(4,5,5,1))
for(i in unique(div_df$global)){
    summ=div_df %>% filter(vary_global & global==i) %>% group_by(timestep) %>% 
        summarize(rich=mean(richness),n=n())
    with(summ %>% filter(timestep !=1),
         plot(rich~timestep,main=paste0('global = ',i,
                                        '\n n = ',n[1])))
}

plot(div_df$richness~div_df$timestep)
plot(div_df$shannon~div_df$timestep)
plot(div_df$simpson~div_df$timestep)

#subsample and calculate div of each community

# first make a dataframe 'coords' with the  coordinates for each individual in a species lists
# start by getting the 1d indices for each individual
dens=density[1];size_xy=size_pix[1] 
individual_indices_r=1:(size_xy^2*dens) #r: indices range from 1:N
individual_indices_p=(1:(size_xy^2*dens))-1 #python: indices range from 0:(n-1)

# from these, get the indices of pixels in 1d list
pixel_indices=as.integer(individual_indices_p/dens)

#then get indices of rows and columns for each pixel
coords1=data.frame(x=pixel_indices%%size_xy,y=as.integer(pixel_indices/size_xy)) 

# next, set the assumed foraging radius (in m)
foraging_radius=500
#calculate the area of a circle with that radius
#and take the square root to get the length of size of a square in m
foraging_dist_square=sqrt(pi*foraging_radius^2)
foraging_dist_pix=foraging_dist_square/30
half_dist_pixels=floor(foraging_dist_pix/2)

# to calculate richness, I'll sample communities in a square in the center of the grid of size pi*foraging_radius^2
# to calculate beta diversity I'll divide the grid into four equidistant quadrants of size pi*foraging_radius^2

# first the middle square for calculating richness: get coords for the square in the center of the grid
center_pix=(size_m[1]/4*2/30)-1
coords1$center_square=F
coords1[coords1$x< (center_pix+half_dist_pixels) & coords1$x > (center_pix-half_dist_pixels) 
& coords1$y < (center_pix+half_dist_pixels) & coords1$y > (center_pix-half_dist_pixels),]$center_square=T
coords1$col_center=ifelse(coords1$center_square,'red','black')

# second, the four squares for calculating beta diversity: divide grid into four equidistant quadrants
# start by getting the center coordinates 
# of each of the 4 quadrants
first_fourth=size_m[1]/4
third_fourth=first_fourth*3
center_coords=data.frame(x_m=c(first_fourth,first_fourth,third_fourth,third_fourth),
                         y_m=c(first_fourth, third_fourth, first_fourth,third_fourth),type='center') %>%
    mutate(x=(x_m/30)-1,y=(y_m/30)-1)

coords=coords1 %>%
    left_join(center_coords) %>% mutate(col=ifelse(is.na(type),'black','red')) %>%
    mutate(index_r=individual_indices_r)


#loop through each row of center_coords and use this get 
# 4 different quadrants in coords df 
coords$quadrant=NA
for(i in 1:nrow(center_coords)){
    x=center_coords[i,]$x;y=center_coords[i,]$y
    coords[coords$x< (x+half_dist_pixels) & coords$x > (x-half_dist_pixels) 
           & coords$y < (y+half_dist_pixels) & coords$y > (y-half_dist_pixels),]$quadrant=i
    
}
coords$col_quadrants=ifelse(is.na(coords$quadrant),'gray',as.factor(coords$quadrant))
coords$col_center=ifelse(coords$center_square,4,'gray')
# plot to make sure the code worked
cex_pt=.7
 pdf('~/Dropbox/Fall_2014/Research/Fragmentation/data/sampling_areas.pdf',width=13)
par(mfrow=c(1,2),pch=15,cex=.8,mar=c(5.1,5.1,4.1,2.1),cex.lab=2.4,cex.axis=2)
with(coords,plot(x,y,col=col_center,asp=1,main="",xlab='x coordinate',ylab='y coordinate'))
with(coords ,plot(x,y,col=col_quadrants,asp=1,xlab='x coordinate',ylab='y coordinate',
                  main = "")) #double check the code worked
 dev.off()

#now sub-sample from the center square to calculate richness
max_time=200000 #only sample from the communities at the last timestep
set.seed(10)
focal_reps=df %>% filter(timestep==max_time)
center_indices=which(coords$center_square)

file_name=focal_reps$filename[1] #for debugging map function
focal_coms=focal_reps$filename %>% map_dfr(function(file_name){
    full_com=read.table(file_name)
    center_com=full_com[center_indices,]
    bee_abund=rtruncnorm(1, a=0, mean = predict_abund , sd = sigma_abund)
    com=center_com[sample(1:length(center_com),size=bee_abund)]
    sad=table(com)
    data.frame(rich=n_distinct(com),shannon=exp(diversity(sad)),
               simpson=diversity(sad,'invsimpson'),filename=file_name)
    
})


spec_df=focal_coms %>% left_join(df) %>% group_by(spec,global) %>%
    summarize(rich=mean(rich),shan=mean(shannon),simp=mean(simpson))%>%
    mutate(vary_spec=global==focal_global,vary_global=spec==focal_spec)
#pdf('methods_param.pdf',width=11)
par(mfrow=c(1,2),pch=16,cex=1.5)
with(spec_df %>% filter(vary_spec),plot(rich~spec,ylim=c(13,30)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)

with(spec_df %>% filter(vary_global),plot(rich~log10(global),ylim=c(13,30)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
#dev.off()

# for each community, take a sample from each quadrant and 
# get mean pairwise jaccard between each
coms_beta=focal_reps$filename %>% map_dfr(function(file_name){
    full_com=read.table(file_name) 
    
    a=1:4 %>% map_dfr(function(i){
        bee_abund=rtruncnorm(1, a=0, mean = predict_abund , sd = sigma_abund)
        
        quad_indices=which(coords$quadrant==i)
        quad_com=full_com[quad_indices,]
        com=quad_com[sample(1:length(quad_com),size=bee_abund)]
        sad=table(com) 
        data.frame(sad) %>% rename(species=com,n=Freq) %>% 
            mutate(quad=i) 
        
    })
    sp_matrix=a %>% pivot_wider(values_from=n,names_from=species,values_fill=0)%>% select(-quad)
    data.frame(bray=mean(vegdist(sp_matrix)),
               jaccard=mean(vegdist(sp_matrix,method='jaccard')),
               horn=mean(vegdist(sp_matrix,method='horn')),
               morisita=mean(vegdist(sp_matrix,method='morisita')),
               filename=file_name)
    
    
})


beta_df=coms_beta %>% left_join(df) %>% group_by(spec,global) %>%
    summarize(bray=mean(bray),jaccard=mean(jaccard),horn=mean(horn),morisita=mean(morisita))%>%
    mutate(vary_spec=global==focal_global,vary_global=spec==focal_spec)

par(mfrow=c(1,2))
with(beta_df %>% filter(vary_spec),plot(log(jaccard)~spec))
with(beta_df %>% filter(vary_global),plot(log(jaccard)~log10(global)))

par(mfrow=c(1,2))
with(beta_df %>% filter(vary_spec),plot(bray~spec))
with(beta_df %>% filter(vary_global),plot(bray~log10(global)))

par(mfrow=c(1,2))
with(beta_df %>% filter(vary_spec),plot(log(horn)~spec))
with(beta_df %>% filter(vary_global),plot(log(horn)~log10(global)))


foraging_radii=c(500,400,300,200,100)
foraging_radius=400
diff_sads=foraging_radii %>% map(function(foraging_radius){
    #calculate the area of a circle with that radius
    #and take the square root to get the length of size of a square in m
    foraging_dist_square=sqrt(pi*foraging_radius^2)
    foraging_dist_pix=foraging_dist_square/30
    half_dist_pixels=floor(foraging_dist_pix/2)
    
    # to calculate richness, I'll sample communities in a square in the center of the grid of size pi*foraging_radius^2
    # to calculate beta diversity I'll divide the grid into four equidistant quadrants of size pi*foraging_radius^2
    
    # first the middle square for calculating richness: get coords for the square in the center of the grid
    center_pix=(size_m[1]/4*2/30)-1
    coords1$center_square=F
    coords1[coords1$x< (center_pix+half_dist_pixels) & coords1$x > (center_pix-half_dist_pixels) 
            & coords1$y < (center_pix+half_dist_pixels) & coords1$y > (center_pix-half_dist_pixels),]$center_square=T
    coords1$col_center=ifelse(coords1$center_square,'red','black')
    
    
    # second, the four squares for calculating beta diversity: divide grid into four equidistant quadrants
    # start by getting the center coordinates 
    # of each of the 4 quadrants
    first_fourth=size_m[1]/4
    third_fourth=first_fourth*3
    center_coords=data.frame(x_m=c(first_fourth,first_fourth,third_fourth,third_fourth),
                             y_m=c(first_fourth, third_fourth, first_fourth,third_fourth),type='center') %>%
        mutate(x=(x_m/30)-1,y=(y_m/30)-1)
    
    coords=coords1 %>%
        left_join(center_coords) %>% mutate(col=ifelse(is.na(type),'black','red')) %>%
        mutate(index_r=individual_indices_r)
    
    
    #loop through each row of center_coords and use this get 
    # 4 different quadrants in coords df 
    coords$quadrant=NA
    for(i in 1:nrow(center_coords)){
        x=center_coords[i,]$x;y=center_coords[i,]$y
        coords[coords$x< (x+half_dist_pixels) & coords$x > (x-half_dist_pixels) 
               & coords$y < (y+half_dist_pixels) & coords$y > (y-half_dist_pixels),]$quadrant=i
        
    }
    coords$col_quadrants=ifelse(is.na(coords$quadrant),'gray',as.factor(coords$quadrant))
    
    
    #now sub-sample from the center square to calculate richness
    max_time=200000 #only sample from the communities at the last timestep
    focal_reps=param_df %>% filter(timestep==max_time)
    center_indices=which(coords$center_square)
    
    file_name=focal_reps$filename[1] #for debugging map function
    focal_coms=focal_reps$filename %>% map(function(file_name){
        full_com=read.table(file_name)
        center_com=full_com[center_indices,]
        
        sad=table(center_com)
        
        
    })
    focal_coms
})

#reorg so that each com gets a df
comp_sads=1:length(diff_sads[[1]]) %>% map(function(com_index){
    
    ls_sad=1:length(foraging_radii) %>% map(function(j){
        sad=diff_sads[[j]][[com_index]]
        ordered_sad=sad[order(sad,decreasing=T)]/sum(sad)*100
        
        sp_names=names(ordered_sad)
        df=data.frame(matrix(ordered_sad,nrow=1))
        colnames(df)=sp_names
        return(df)
    }) #%>% bind_rows
    df_sad=bind_rows(ls_sad)
    df_sad[is.na(df_sad)] <-0
    return(df_sad)
})
community=comp_sads[[1]]
cols_rgb=col2rgb(c("lightblue", "lightgreen", "pink",'orange','darkmagenta'))
cols_to_use=1:ncol(cols_rgb) %>% map(function(i){
    col_vals=cols_rgb[,i]
    rgb(col_vals[1],col_vals[2],col_vals[3],max = 255, alpha = 80)

        }) %>% unlist

par(mfrow=c(1,1))
for(community in comp_sads){
    #min_abund=min(community);max_abund=max(community)
    break_vec=seq(0,100,5)
    
    for(i in 1:nrow(community)){
        my_col=cols_to_use[i]
        my_com=as.numeric(community[i,])
        if(i==1) hist(my_com,breaks=break_vec,col=my_col,main="",xlab='relative abund',ylim=c(0,30))
        if(i !=1) hist(my_com,breaks=break_vec,add=T,col=my_col,main="",xlab='relative abund')
    }
}
test %>% bind_rows

for(com_index in 1:length(diff_sads[[1]])){
    for(j in 1:length(foraging_radius)){
        #get 
        diff_sads[[j]][[com_index]]
    
        }
    sad_one_radius=diff_sads[[i]]
    for(j in length)
}

## old
par(mfrow=c(2,5))
for(i in unique(div_df$rep)){
    for(j in unique(div_df$spec)){
        summ=div_df %>% filter(rep==i & spec==j) 
        with(summ %>% filter(timestep !=1),
             plot(richness~timestep,main=paste0("spec=",j,", rep=",i)))   
    }
}
plot(div_df$richness~div_df$timestep)
plot(div_df$shannon~div_df$timestep))
plot(div_df$simpson~div_df$timestep))

size_xy=110;density=13
0:(length(sp$V1)-1) %>% map_dfr(function(i){
   
    data.frame(x=(i%%size_xy),y=floor(i/size_xy))
})
0:10 %>% map(function(i){
    
    data.frame(x=(i%%size_xy),y=floor(i/size_xy))
})
#python code for getting coords
# def get_coords(species,size_xy,density):
#     '''this function takes a list of species from the neutral model with
#     size_xy*size_xy pixels and 'density' individuals per pixel and returns
#     the x and y coordinates of each individual. the total
#     size of the list should be size**2*density '''
# 
#     assert len(species) == size_xy**2*density, "length of species list should be size_xy^2*density"
#     com=np.arange(size_xy**2*density)
#     com_pixel_coords=[int(i/density) for i in com]
# 
#     xs=[i%size_xy for i in range(size_xy**2)]
#     ys=[int(i/size_xy) for i in range(size_xy**2)]
# 
#     com_xs=[xs[i] for i in com_pixel_coords]
#     com_ys=[ys[i] for i in com_pixel_coords]
#     locations={'species': species,"x": com_xs,"y": com_ys, 
#         'coordpixel_1d': com_pixel_coords,"cordindividual_1d": com}
#     return(pd.DataFrame(locations))