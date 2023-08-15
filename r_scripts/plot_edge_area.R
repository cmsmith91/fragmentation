rm(list=ls())
# setwd("~/Dropbox/Fall_2014/research/fragmentation/data")
library(tidyverse)
map=purrr::map; select=dplyr::select
lu_df=readRDS('processed_data/lu_info03dec2021.rds') %>% 
    map(function(df)df %>% mutate(size_ha=size_m^2*.0001,edge_cat2=ifelse(is.na(edge_cat),'none',edge_cat)))

#how many landscapes at each spatial scale
lu_df %>% map(nrow)
plotting_cols=cols_for_map=c(adjustcolor(RColorBrewer::brewer.pal(5,'Dark2')[c(2,5)],.7),
                             adjustcolor('black',.4))

lu_df_all =lu_df %>% bind_rows

df=lu_df[[1]]
legend_labs=c('high edge','low edge','other')

# pdf('figures/edge_v_area.pdf',width=16)
tiff('figures/edge_v_area.tiff', units="in", width=18, height=8, res=100, compression = 'lzw')

cex_mtext=1.9
par(pch=16,mfrow=c(1,5),cex.main=2.5,cex.lab=2,cex.axis=2,mar=c(3.5,3,2,2),oma=c(2,2,1,1))
for(df in lu_df){
    plot_df=df %>% arrange(focal_landscape)
    with(plot_df,plot(prop_forest*100,edge_density,main=paste0(round(df$size_ha[1]),' ha'),
                      xlab='',ylab='',ylim=c(0,190),
                      pch=c(16,17)[as.factor(focal_landscape)],cex=c(1,2)[as.factor(focal_landscape)],
                      col=plotting_cols[as.factor(edge_cat2)]))
    if(df$size_m[1]==600)  mtext('edge density (m/ha)', side = 2, line = 0, outer = T,cex=cex_mtext,adj=.5)
    if(df$size_m[1]==1740) mtext('% forest', side = 1, line = 0, outer = T,cex=cex_mtext)
    if(df$size_m[1]==2400) legend('topright', inset = -0.03,cex=2,legend=legend_labs,title = '',pch=c(17,17,16)[as.factor(legend_labs)],col=plotting_cols[as.factor(legend_labs)],bty='n')


}
 dev.off()

#results summary
(max_edge_land=lu_df_all %>% filter(edge_density==max(edge_density)))
paste("Across the", nrow(lu_df_all),"landscapes we sampled in the northeastern U.S., forest edge 
density varied between a minimum of",min(lu_df_all$edge_density),"m/ha, in completely forested and deforested
landscapes, to a maximum of",round(max_edge_land$edge_density[1],1),"m/ha in",nrow(max_edge_land), "xxxxx ha landscape with xxx
      % forest") 

#with forest edge at its maximum when forests occupied a little over xx%  of the landscape, 
df=lu_df[[1]] %>% mutate(forest_squared=prop_forest^2)

lu_df %>% map_dbl(function(df){
    coef_quad_mod=with(df %>% mutate(forest_squared=prop_forest^2),coef(lm(edge_density~prop_forest+forest_squared)))
    xs=seq(0,1,.0001)
    ys=coef_quad_mod[1]+xs*coef_quad_mod[2]+xs^2*coef_quad_mod[3]
    
    #plot(xs,ys)
    xs[ys==max(ys)]
    
})

##
(smallest50=lu_df[[1]] %>% filter(f_cat==.1) %>% summarize(min_edge=min(edge_density),max_edge=max(edge_density)))
var_smallest=smallest50$max_edge/smallest50$min_edge

(largest50=lu_df[[5]] %>% filter(f_cat==.1) %>% summarize(min_edge=min(edge_density),max_edge=max(edge_density)))
var_largest=largest50$max_edge/largest50$min_edge


paste("At the smallest spatial scale (36 ha) forest edge density varied by a factor of", round(var_smallest,1) ,
"in landscapes with 10% (± 2.5%) forest habitat, whereas at the largest spatial 
scale (576 ha), it varied by a factor of", round(var_largest,1) ,"in landscapes with the same percentage of forest habitat.") 

