rm(list=ls())
# setwd("~/Dropbox/Fall_2014/Research/Fragmentation/data/neutralmod_output_final")
setwd("~/Dropbox/Fall_2014/Research/Fragmentation/data/neutral_output_final_again")
setwd("/Users/colleen/Dropbox/Fall_2014/Research/Fragmentation/data/neutralmod_output_13dec2021-4")
setwd("/Users/colleen/Dropbox/Fall_2014/Research/Fragmentation/data/neutralmod_output_13dec2021-4")

library(truncnorm)
library(tidyverse)
library(vroom)
library(vegan)
library(furrr)
map=purrr::map; select=dplyr::select

## get a vector of all the files
files=list.files()
if(grepl('txt',files)[1]==F){
    files=files %>% map(function(f_name){
        paste0(f_name,"/",list.files(f_name))
    }) %>% unlist
}
#get rid of slurm files
files=files[!grepl('slurm',files)]

## make column for the files with proportions (have "_props" in file name)
## vs files with locations of each individual
global=as.numeric(gsub("_size.*","",gsub(".*_global",'',files)))
rep=as.numeric(gsub('.*rep','',gsub('timesteps.*','',files)))
spec=as.numeric(gsub("_global.*","",gsub('.*spec','',files)))
timestep=as.numeric(gsub('\\.txt.*','',gsub('.*timesteps','',files)))
timestep[is.na(timestep)]=as.numeric(gsub('_props.*','',gsub('.*timesteps','',files[is.na(timestep)])))
density=as.numeric(gsub('_rep.*','',gsub('.*density','',files)))
size_pix=as.numeric(gsub('_den.*','',gsub('.*sizexy','',files)))
size_ha=(size_pix*30)^2/10000
size_m=size_pix*30
props=grepl("_props",files)
df=data.frame(spec=spec,global=global,size_pix=size_pix,size_ha=size_ha,size_m=size_m,
              density=density,rep=rep,timestep=timestep,props=props,filename=files)
param_df=df
param_sads=df %>% filter(props)
param_coms=df %>% filter(!props)
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

#get number of individuals, species and % id'd to the species level
fbee_sad=fbees_allsites %>% group_by(genus_species) %>% summarize(n=n())
print(paste("the number of forest-associated bee species in the data is",nrow(fbee_sad)))
print(paste("the number of forest-associated bee individuals in the data is",sum(fbee_sad$n)))
fbee_sad[fbee_sad$genus_species=="Nomada_bidentate_group",]$n/sum(fbee_sad$n)*100
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
predict_abund/(pi*400^2/10000)
predict_abund/(pi*185^2/10000)
21*.09
# load simulated data and look at how div changes with each timestep
#load the props_files
plan(multisession,workers=6)
sads=param_sads$filename %>% 
    future_map(function(file_name) {
        df=vroom(file_name,delim="\"",col_names = F)
        tab=df$X1
        data.frame(
            richness=length(tab),
            shannon=exp(diversity(tab)),
            simpson=diversity(tab,'invsimpson'))
        })

focal_spec=as.numeric(names(table(spec)[order(table(spec),decreasing=T)])[1])
focal_global=as.numeric(names(table(global)[order(table(global),decreasing=T)])[1])

div_df=sads %>% bind_rows %>% bind_cols(param_sads) %>%
    arrange(desc(spec),desc(global)) %>%
    mutate(vary_spec=round(global,5)==round(focal_global,5),vary_global=spec==focal_spec)

div_df %>% filter(vary_spec &vary_global)

par(mfrow=c(1,2),mar=c(4,5,5,1))
for(i in unique(div_df$spec)){
    summ=div_df %>% filter(vary_spec & spec==i) %>% group_by(timestep) %>% 
        summarize(rich=mean(richness),simp=mean(simpson),n=n())
    if(nrow(summ)>0){
        with(summ %>% filter(timestep !=1),
             plot(rich~timestep,ylim=c(0,max(rich)+max(rich)*.10),main=paste0('spec = ',i,
                                            '\n n = ',n[1])))
        with(summ %>% filter(timestep !=1),
             plot(simp~timestep,ylim=c(0,max(simp)+max(simp)*.10),main=paste0('spec = ',i,
                                                                              '\n n = ',n[1])))
    }
    
}

par(mfrow=c(1,2),mar=c(4,5,5,1),pch=16)
globals=unique(div_df$global) 
globals_ordered=globals[order(globals)]
for(i in globals_ordered){
    summ=div_df %>% filter(vary_global & global==i) %>% group_by(timestep) %>% 
        summarize(rich=mean(richness),simp=mean(simpson),n=n())
    with(summ %>% filter(timestep !=1),
         plot(rich~timestep,ylim=c(0,max(rich)+max(rich)*.10),main=paste0('global = ',round(i,6),
                                        '\n n = ',n[1])))
}



##############################
#next, subsample and calculate div of each community

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
foraging_radius=185
#calculate the area of a circle with that radius
#and take the square root to get the length of size of a square in m
foraging_dist_square=sqrt(pi*foraging_radius^2)
foraging_dist_pix=foraging_dist_square/30
half_dist_pixels=ceiling(foraging_dist_pix/2)

# to calculate richness, I'll sample communities in a square in the center of the grid of size pi*foraging_radius^2
# to calculate beta diversity I'll divide the grid into four equidistant quadrants of size pi*foraging_radius^2

# first the middle square for calculating richness: get coords for the square in the center of the grid
# center_pix=(size_m[1]/4*2/30)-1 #shouldn't this be size_pix[1]/2-1
center_pix=size_pix[1]/2-1
coords1$center_square=F
coords1[coords1$x< (center_pix+half_dist_pixels) & coords1$x >= (center_pix-half_dist_pixels) 
        & coords1$y < (center_pix+half_dist_pixels) & coords1$y >= (center_pix-half_dist_pixels),]$center_square=T
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
    coords[coords$x< (x+half_dist_pixels) & coords$x >= (x-half_dist_pixels) 
           & coords$y < (y+half_dist_pixels) & coords$y >= (y-half_dist_pixels),]$quadrant=i
    
}
coords$col_quadrants=ifelse(is.na(coords$quadrant),'gray',as.factor(coords$quadrant))
coords$col_center=ifelse(coords$center_square,4,'gray')
# plot to make sure the code worked
cex_pt=.7
# pdf('~/Dropbox/Fall_2014/Research/Fragmentation/data/sampling_areas.pdf',width=13)
par(mfrow=c(1,2),pch=15,cex=.8,mar=c(5.1,5.1,4.1,2.1),cex.lab=2.4,cex.axis=2)
with(coords,plot(x,y,col=col_center,asp=1,main="",xlab='x coordinate',ylab='y coordinate'))
with(coords ,plot(x,y,col=col_quadrants,asp=1,xlab='x coordinate',ylab='y coordinate',
                  main = "")) #double check the code worked
# dev.off()

sum(coords$center_square)


#now sub-sample from the center square to calculate richness
max_time=150000 #only sample from the communities at the last timestep
set.seed(10)
focal_reps=param_coms %>% filter(timestep==max_time)
center_indices=which(coords$center_square)

file_name=focal_reps$filename[1] #for debugging map function
focal_coms=focal_reps$filename %>% map_dfr(function(file_name){
    full_com=vroom(file_name,delim="\"",col_names = F)
    
    bee_abund=rtruncnorm(1, a=0, mean = predict_abund , sd = sigma_abund)
    half_dist_pixels=ceiling(sqrt(bee_abund/density[1])/2)
    center_indices=which(coords$x< (center_pix+half_dist_pixels) & coords$x >= (center_pix-half_dist_pixels) 
                         & coords$y < (center_pix+half_dist_pixels) & coords$y >= (center_pix-half_dist_pixels))
    
    center_com=full_com[center_indices,]
    com=center_com[sample(1:nrow(center_com),size=bee_abund),]
    sad=table(com)
    data.frame(rich=n_distinct(com),shannon=exp(diversity(sad)),
               simpson=diversity(sad,'invsimpson'),filename=file_name)
    
})


spec_df=focal_coms %>% left_join(param_coms) %>% 
  group_by(spec,global) %>%
    summarize(rich=median(rich),shan=median(shannon),simp=median(simpson),n=n())%>%
    mutate(vary_spec=round(global,5)==round(focal_global,5),vary_global=spec==focal_spec)

#pdf('methods_param.pdf',width=11)
par(mfrow=c(1,2),pch=16,cex=1.5)
with(spec_df %>% filter(vary_spec ),plot(rich~spec,ylim=c(0,40)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
with(spec_df %>% filter(spec==0.000274250 &vary_spec),points(spec,rich,col='blue'))
with(spec_df %>% filter(round(spec,5)==round(0.000128750,5)  &vary_spec),points(spec,rich,col='red'))

with(spec_df %>% filter(vary_global),plot(rich~log10(global),ylim=c(0,40)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
#dev.off()
analyze_me=spec_df %>% filter(vary_spec)
spec_mod=lm(analyze_me$rich~analyze_me$spec)
make_focal_spec=(exp(predicted_richness)-coef(spec_mod)[1])/coef(spec_mod)[2]

#next measure beta div
file_name=focal_reps$filename[[1]]
plan(multisession,workers=6)
beta_div=focal_reps$filename %>% future_map(function(file_name){
    full_com=vroom(file_name,delim="\"",col_names = F)
    sub_coms=1:nrow(center_coords) %>% map_dfr(function(i){
        bee_abund=rtruncnorm(1, a=0, mean = predict_abund , sd = sigma_abund)
        half_dist_pixels=ceiling(sqrt(bee_abund/density[1])/2)
        
        x=center_coords[i,]$x;y=center_coords[i,]$y
        
        center_indices=which(coords$x< (x+half_dist_pixels) & coords$x >= (x-half_dist_pixels) 
                             & coords$y < (y+half_dist_pixels) & coords$y >= (y-half_dist_pixels))
        center_com=full_com[center_indices,] 
        quadrant_com=center_com[sample(1:nrow(center_com),size=bee_abund),] %>% mutate(site=i)
        
    })
    com_matrix=sub_coms %>% group_by(site,X1) %>% summarize(n=n()) %>%
        pivot_wider(values_from=n,names_from=X1,values_fill = 0) 
    com_matrix=com_matrix[,!names(com_matrix)=='site']
    data.frame(filename=file_name,jaccard=mean(vegdist(com_matrix,method='jaccard')))
},.options=furrr_options(seed=T)) %>% bind_rows%>% left_join(param_coms)


beta_df=beta_div %>% group_by(spec,global) %>% 
  summarize(jaccard=median(jaccard),n=n())%>% left_join(spec_df)


#add observed jaccard to plot
head(fbees_allsites)
head(sites)
big_forests=sites[sites$forest500>.8,]$site #just pick sites with >.8 forest cover
fbee_matrix=fbees_allsites %>% filter(site %in% big_forests) %>%
    group_by(site,genus_species) %>% summarize(n=n())%>%
    pivot_wider(values_from=n,names_from=genus_species,values_fill = 0) 
fbee_matrix=fbee_matrix[,!names(fbee_matrix)=='site']
obs_jaccard_vec=as.vector(vegdist(fbee_matrix,method='jaccard'))
obs_jaccard=mean(obs_jaccard_vec)
lower_jacc=quantile(obs_jaccard_vec,.1)
upper_jacc=quantile(obs_jaccard_vec,.9)

#pdf('methods_param.pdf',width=11)
par(mfrow=c(1,2),pch=16,cex=1.5)
with(beta_df %>% filter(vary_spec),plot(jaccard~spec))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
with(beta_df %>% filter(spec==0.000274250),points(jaccard~spec,col='red'))

with(beta_df %>% filter(vary_global),plot(jaccard~log10(global),ylim=c(0.5,.9)))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
#dev.off()


#next, figure out which values of spec and dispersal to draw from
#first richness
pick_spec_rich=spec_df %>% filter(vary_spec) %>% arrange(rich)
in_btw=pick_spec_rich$rich >=lower_rich & pick_spec_rich$rich <=upper_rich#in btw lower and upper richness + 1 on each side
fit_me=in_btw
for(i in 1:length(in_btw)) {
    if(in_btw[i] & length(in_btw[i-1]) !=0) fit_me[i-1] <- T
    if(in_btw[i] & length(in_btw[i+1]) !=0) fit_me[i+1] <- T
}
in_btw_sr=fit_me
coefs_spec_rich=with(pick_spec_rich[fit_me,],coef(lm(rich~spec)))


pick_global_rich=spec_df %>% filter(vary_global)
in_btw=pick_global_rich$rich >=lower_rich & pick_global_rich$rich <=upper_rich#in btw lower and upper richness + 1 on each side
fit_me=in_btw
for(i in 1:length(in_btw)) {
    if(in_btw[i] & length(in_btw[i-1]) !=0) fit_me[i-1] <- T
    if(in_btw[i] & length(in_btw[i+1]) !=0) fit_me[i+1] <- T
}
in_btw_gr=fit_me
coefs_global_rich=with(pick_global_rich[fit_me,],coef(lm(rich~log10(global))))


#get x vals for plotting the bets fit liinw

spec_xs=seq(min(pick_spec_beta$spec),max(pick_spec_beta$spec),by=(max(pick_spec_beta$spec)-min(pick_spec_beta$spec))/1000)

xs_sr=pick_spec_rich[in_btw_sr,]$spec#in_btw_gr
spec_xs_rich=seq(min(xs_sr),max(xs_sr),by=(max(xs_sr)-min(xs_sr))/1000)
spec_ys_rich=coefs_spec_rich[1]+coefs_spec_rich[2]*spec_xs_rich

xs_gr=pick_global_rich[in_btw_gr,]$global
spec_xs_global=log10(seq(min(xs_gr),max(xs_gr),by=(max(xs_gr)-min(xs_gr))/1000))
spec_ys_global=coefs_global_rich[1]+coefs_global_rich[2]*spec_xs_global


#now plot the fitted lines to the data
par(mfrow=c(1,2),pch=16,cex=1.5)
with(spec_df %>% filter(vary_spec) ,plot(rich~spec,ylim=c(0,30)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
lines(spec_xs_rich,spec_ys_rich,type='l',col='blue')
# abline(a=coefs_spec_rich[1],b=coefs_spec_rich[2],col='blue')
# abline(v=final_min_spec)
# abline(v=final_max_spec)

with(spec_df %>% filter(vary_global),plot(rich~log10(global),ylim=c(0,30)))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
# abline(a=coefs_global_rich[1],b=coefs_global_rich[2],col='blue')
lines(spec_xs_global,spec_ys_global,type='l',col='blue')
# abline(v=vals_global_rich[1],col='green')
# abline(v=final_min_global)
# abline(v=final_max_global)

#dev.off()

#solve equations for where speciation and richness cross upper and lower richness
solve_me=function(intercept,slope,y_val){
    (y_val-intercept)/slope
    
}
vals_spec_rich=solve_me(coefs_spec_rich[1],coefs_spec_rich[2],c(lower_rich,upper_rich))
vals_global_rich=solve_me(coefs_global_rich[1],coefs_global_rich[2],c(lower_rich,upper_rich))


###
#next do the same thing but for beta diversity

pick_spec_beta=beta_df %>% filter(vary_spec) %>% arrange(jaccard)

in_btw=pick_spec_beta$jaccard >=lower_jacc & pick_spec_beta$jaccard <=upper_jacc#in btw lower and upper richness + 1 on each side
fit_me=in_btw
for(i in 1:length(in_btw)) {
    if(in_btw[i] & length(in_btw[i-1]) !=0) fit_me[i-1] <- T
    if(in_btw[i] & length(in_btw[i+1]) !=0) fit_me[i+1] <- T
}
in_btw_sj=fit_me
coefs_spec_beta_loglog=with(pick_spec_beta[fit_me,],coef(lm(log(jaccard)~log(spec))))

xs_sj=pick_spec_beta$spec[in_btw_sj]
spec_xs=seq(min(xs_sj),max(xs_sj),by=(max(xs_sj)-min(xs_sj))/1000)
spec_ys=exp(coefs_spec_beta_loglog[1]+log(spec_xs)*coefs_spec_beta_loglog[2])

#do global div next
pick_global_beta=beta_df %>% filter(vary_global) %>% arrange(jaccard)

in_btw=pick_global_beta$jaccard >=lower_jacc & pick_global_beta$jaccard <=upper_jacc#in btw lower and upper richness + 1 on each side
fit_me=in_btw
for(i in 1:length(in_btw)) {
    if(in_btw[i] & length(in_btw[i-1]) !=0) fit_me[i-1] <- T
    if(in_btw[i] & length(in_btw[i+1]) !=0) fit_me[i+1] <- T
}
in_btw_gj=fit_me
coefs_global_beta=with(pick_global_beta[fit_me,],coef(lm(jaccard~log10(global))))

xs_gj=pick_global_beta[in_btw_gj,]$global
glob_xs=log10(seq(min(xs_gj),10^(-1.9),by=(10^(-1.9)-min(xs_gj))/1000))
glob_ys=coefs_global_beta[1]+coefs_global_beta[2]*glob_xs

vals_spec_beta=exp(solve_me(coefs_spec_beta_loglog[1],coefs_spec_beta_loglog[2],c(log(lower_jacc),log(upper_jacc))))
vals_global_beta=solve_me(coefs_global_beta[1],coefs_global_beta[2],c(lower_jacc,upper_jacc))

#plot
graph_lims=c(.38,.93)
par(mfrow=c(1,2),pch=16,cex=1.5)
with(beta_df %>% filter(vary_spec& spec !=0.00022575),
     plot(jaccard~spec,ylim=graph_lims))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
# abline(a=coefs_spec_beta[1],b=coefs_spec_beta[2],col='blue')
lines(spec_xs,spec_ys,type='l',col='blue')
# abline(v=vals_spec_beta[1],col='blue',lt=2)
# abline(v=vals_spec_beta[2],col='blue',lt=2)
# abline(v=final_min_spec)
# abline(v=final_max_spec)

with(beta_df %>% filter(vary_global),
     plot(jaccard~log10(global),ylim=graph_lims,xlim=c(-3.1,-1.7)))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
lines(glob_xs,glob_ys,type='l',col='blue')
# abline(a=coefs_global_beta[1],b=coefs_global_beta[2],col='blue')
# abline(v=vals_global_beta[1],col='blue',lt=2)
# abline(v=vals_global_beta[2],col='blue',lt=2)
# abline(v=final_min_global)
# abline(v=final_max_global)



#make pdf for supplement
pdf('param_exploration.pdf',width=11)
par(mfrow=c(2,2),pch=16,cex=1.5,mar=c(2,4,1,.5))
with(spec_df %>% filter(vary_spec & spec !=0.00022575) ,
     plot(rich~spec,ylim=c(0,30),xlab="",ylab='species richness'))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
lines(spec_xs_rich,spec_ys_rich,type='l',col='blue')

par(mar=c(2,2,1,.5))
with(spec_df %>% filter(vary_global),
     plot(rich~log10(global),ylim=c(0,30),xlab="",ylab=""))
abline(h=exp(predicted_richness),col='red')
abline(h=(lower_rich),col='red',lt=2)
abline(h=(upper_rich),col='red',lt=2)
lines(spec_xs_global,spec_ys_global,type='l',col='blue')

#plot long-distanece dispersal
par(mar=c(4.3,4,1,.5))
with(beta_df %>% filter(vary_spec& spec !=0.00022575),
     plot(jaccard~spec,ylim=graph_lims,xlab='speciation'))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
lines(spec_xs,spec_ys,type='l',col='blue')

par(mar=c(4.3,2,1,.5))
with(beta_df %>% filter(vary_global),
     plot(jaccard~log10(global),ylim=graph_lims,xlim=c(-3.1,-1.7),
          xlab="log10(long-range dispersal)",ylab=""))
abline(h=obs_jaccard,col='red')
abline(h=(lower_jacc),col='red',lt=2)
abline(h=(upper_jacc),col='red',lt=2)
lines(glob_xs,glob_ys,type='l',col='blue')
dev.off()


#summarize everything
print(paste0('for the communities to match the richness of bee communities, the speciation value
      should be between ',vals_spec_rich[1],' and ',vals_spec_rich[2]))
print(paste0('for the communities to match the beta diversity of bee communities, the speciation value
      should be between ',vals_spec_beta[1],' and ',vals_spec_beta[2]))

print(paste0('for the communities to match the richness of bee communities, the global value
      should be between 10**',vals_global_rich[1],' and 10**',vals_global_rich[2]))
print(paste0('for the communities to match the beta diversity of bee communities, the global value
      should be between 10**',vals_global_beta[1],' and 10**',vals_global_beta[2]))


#we'll use the [minimum-bounding box] of these values
#speciation to use
final_min_spec=max(vals_spec_rich[1],vals_spec_beta[1])
final_max_spec=min(vals_spec_rich[2],vals_spec_beta[2])
final_max_global=min(vals_global_rich[1],vals_global_beta[1])
final_min_global=max(vals_global_rich[2],vals_global_beta[2])

print(paste0("the minimum speciation value we'll use is ",final_min_spec," and
             the maximum speciation value we'll use is ",final_max_spec))
print(paste0("the minimum global value we'll use is 10^",final_min_global," and
             the maximum global value we'll use is 10^",final_max_global))

# (lower_rich-coefs_spec_rich[1])/coefs_spec_rich[2]
#plot everything

min(spec)
unique(beta_df$spec)

check=with(spec_df %>% filter(vary_spec), lm(rich~spec))
(exp(predicted_richness)-coef(check)[1])/coef(check[2])


exp(predicted_richness)-coef(check)[1]
6.24059 /48842.59259 


# ##check change in beta_div over time
# beta_div_time=param_coms$filename %>% future_map(function(file_name){
#     full_com=vroom(file_name,delim="\"",col_names = F)
#     sub_coms=1:nrow(center_coords) %>% map_dfr(function(i){
#         bee_abund=rtruncnorm(1, a=0, mean = predict_abund , sd = sigma_abund)
#         half_dist_pixels=ceiling(sqrt(bee_abund/density[1])/2)
#         
#         x=center_coords[i,]$x;y=center_coords[i,]$y
#         
#         center_indices=which(coords$x< (x+half_dist_pixels) & coords$x >= (x-half_dist_pixels) 
#                              & coords$y < (y+half_dist_pixels) & coords$y >= (y-half_dist_pixels))
#         center_com=full_com[center_indices,] 
#         quadrant_com=center_com[sample(1:nrow(center_com),size=bee_abund),] %>% mutate(site=i)
#         
#     })
#     com_matrix=sub_coms %>% group_by(site,X1) %>% summarize(n=n()) %>%
#         pivot_wider(values_from=n,names_from=X1,values_fill = 0) 
#     com_matrix=com_matrix[,!names(com_matrix)=='site']
#     data.frame(filename=file_name,jaccard=mean(vegdist(com_matrix,method='jaccard')))
# },.options=furrr_options(seed=T)) %>% bind_rows%>% left_join(param_coms)
# 
# 
# beta_df_time=beta_div_time %>% group_by(spec,global,timestep) %>% summarize(jaccard=mean(jaccard),n=n())%>% left_join(spec_df)
# 
# ##plot
# 
# par(mfrow=c(1,2),mar=c(4,5,5,1))
# for(i in unique(beta_df_time$spec)){
#     summ=beta_df_time %>% filter(vary_spec & spec==i) %>% group_by(timestep) %>% 
#         summarize(jacc=mean(jaccard))
#     if(nrow(summ)>0){
#         with(summ %>% filter(timestep !=1),
#              plot(jacc~timestep,main=paste0('spec = ',i)))
#        
#     }
#     
# }
# 
# par(mfrow=c(1,2),mar=c(4,5,5,1))
# globals=unique(div_df$global) 
# globals_ordered=globals[order(globals)]
# for(i in globals_ordered){
#     summ=beta_df_time %>% filter(vary_global & global==i) %>% group_by(timestep) %>% 
#         summarize(jacc=mean(jaccard))
#     with(summ %>% filter(timestep !=1),
#          plot(jacc~timestep,main=paste0('global = ',round(i,6))))
# }
# 
# 
# 
