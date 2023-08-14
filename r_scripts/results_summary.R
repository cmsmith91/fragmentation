rm(list=ls())
library(tidyverse)
map=purrr::map

#load data
output=readRDS("processed_data/loss_simulation_n100_output.rds")
edge=readRDS('processed_data/lu_info03dec2021.rds')
focal_edge=edge %>% 
    map_dfr(function(df) df %>% filter(focal_landscape)) %>%
    mutate(landscape_id=1:n())

#format data
output_summ=output %>% group_by(grid_index,size_m) %>%
    summarize(mean_percent_loss=mean(percent_loss),
              median_percent_loss=median(percent_loss),mean_finalN=mean(final_n),median_finalN=median(final_n)) %>% 
    left_join(focal_edge) %>% 
    arrange(size_m) %>%
    mutate(size_ha=size_m^2*0.0001)%>%mutate(median_percent_remaining=100-median_percent_loss)

set.seed(99)
# Get the slope and intercepts of the relationships between 
# edge density and % species remaining.
my_size = unique(output_summ$size_m)[1]
my_f_cat = unique(output_summ$f_cat)[1]
size_df=unique(output_summ$size_m) %>% map_dfr(function(my_size){
    df=output_summ %>% filter(size_m==my_size)
    unique(df$f_cat) %>% map_dfr(function(my_f_cat){
        new_df = df %>% filter(f_cat==my_f_cat) %>%
            mutate(percent_remaining=100-median_percent_loss)
        mod=with(new_df,lm(percent_remaining~edge_density))
        edge_relationship=coef(mod)
        conf_ints=confint(mod,"edge_density")
        lower=conf_ints[1]
        upper=conf_ints[2]
        
        tibble(size_m=my_size,size_ha=my_size^2/10000,
               f_cat=my_f_cat,
               edge_slope=edge_relationship[2],
               mod_intercept=edge_relationship[1],
               lower_ci=lower,upper_ci=upper)
    })
    
}) %>% mutate(size_ha_jit=jitter(size_ha))


#measure change as average diff
diff_df=output_summ %>% split(.$size_ha) %>% map_dfr(function(df){
    df %>% split(.$f_cat) %>%  map_dfr(function(df2){
        diff_mod=with(df2 %>% mutate(edge_cat=factor(edge_cat,levels=c('low','high'))),
             lm(median_percent_remaining~edge_cat)
             )
        avg_diff=coef(diff_mod)[2]
        cis=confint(diff_mod,'edge_cathigh')
        
        
        data.frame(size_ha=df2$size_ha[1],f_cat=df2$f_cat[1],
                   avg_diff=avg_diff,lower_ci_diff=cis[1],upper_ci_diff=cis[2])
        })
}) %>% left_join(size_df)

with(diff_df,plot(size_ha_jit,avg_diff,col=as.factor(f_cat),pch=16,ylim=c(0,8)))
with(diff_df,arrows(size_ha_jit,lower_ci_diff,y1=upper_ci_diff,angle=90,length=0,col=as.factor(f_cat)))

summary(lm(avg_diff~size_ha + f_cat,data=diff_df))


#min and max differences btw high and low-edge landscapes
max(diff_df$avg_diff);min(diff_df$avg_diff)
diff_df %>% filter(f_cat==0.1) %>% summarize(mean(avg_diff))
diff_df %>% filter(f_cat==0.3) %>% summarize(mean(avg_diff))
diff_df %>% filter(f_cat==0.5) %>% summarize(mean(avg_diff))

#min and max differences at the different scales
diff_df %>% filter(size_ha==min(size_ha)) %>% summarize(mean(avg_diff))
diff_df %>% filter(size_ha==max(size_ha)) %>% summarize(mean(avg_diff))

#differences in the slope between landscape sizes
size_df %>% filter(size_ha==min(size_ha))%>% summarize(mean(edge_slope))
size_df %>% filter(size_ha==max(size_ha))%>% summarize(mean(edge_slope))

#absolute range?
range(size_df$edge_slope)
#what about differences in population size?
output_summ %>% group_by(edge_cat) %>% summarize(mean(mean_finalN))


my_cols=RColorBrewer::brewer.pal(8,'Paired')[c(7,2,4)]
my_cols_lines=my_cols
my_cols=adjustcolor(my_cols,alpha=.67)

plot_df=output_summ %>% filter(size_m==600)
legend_labs=unique(plot_df$f_cat)[order(unique(plot_df$f_cat))]
legend_labs2=paste0(legend_labs*100,'%')
my_pch = c(16,17,15)

# pdf('figures/size_v_effectedge3.pdf',width=12)
#tiff('figures/size_v_effectedge.tiff', units="in", width=12, height=8, res=500, compression = 'lzw')
par(mfrow=c(1,2),cex=1.6,mar=c(4.5,4,2,.9),oma=c(.2,.2,.2,.2))
with(size_df,plot(size_ha_jit,edge_slope,col=my_cols_lines[as.factor(f_cat)],ylim=c(0,.37),ylab="edge density effect: slope",xlab="landscape size (ha)",pch=my_pch[as.factor(f_cat)]))
with(size_df,arrows(x0=size_ha_jit,y0=lower_ci,y1=upper_ci,angle=90,length=0,col=my_cols_lines[as.factor(f_cat)]))
legend('topright',cex=1,legend=legend_labs2,title = 'forest remaining',col=my_cols_lines[as.factor(legend_labs)],bty='n',pch=my_pch[as.factor(legend_labs)])

with(diff_df,plot(size_ha_jit,avg_diff,col=my_cols_lines[as.factor(f_cat)],ylim=c(0,10),ylab="edge density effect: mean difference",xlab="landscape size (ha)",pch=my_pch[as.factor(f_cat)]))
with(diff_df,arrows(size_ha_jit,lower_ci_diff,y1=upper_ci_diff,angle=90,length=0,col=my_cols_lines[as.factor(f_cat)]))
# dev.off()

# pdf("figures/loss_sim_output2-08march2022.pdf",width=18)
# tiff('figures/loss_sim_output2-08march2022.tiff', units="in", width=20, height=10, res=700, compression = 'lzw')
tiff('figures/loss_sim_output2-14aug2023.tiff', units="in", width=20, height=10, res=100, compression = 'lzw')

cex_lab=1.3
par(mfrow=c(1,5),pch=16,cex=1.4,cex.main=1.3,cex.lab=cex_lab,mar=c(4.1,2.5,2,.5))
for(my_size in unique(output_summ$size_m)){
    plot_df=output_summ %>% filter(size_m==my_size)
    size_plot=size_df %>% filter(size_m==my_size)
    
    
    
    if(my_size==600){
        par(oma=c(.5,1.5,1,1))
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],pch=my_pch[as.factor(f_cat)],
                  main=paste0(size_ha[1],' ha'),ylim=c(50,100),xlim=c(10,200),
                  ylab='',xlab='')
        )
        mtext('% species remaining', side = 2, line = 0, outer = T,cex=1.9,adj=.6)
        legend('bottomright',cex=1,legend=legend_labs2,title = 'forest remaining',col=my_cols[as.factor(legend_labs)],pch=my_pch[as.factor(legend_labs)],bty='n')
    }
    if(my_size==1740){
        
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],pch=my_pch[as.factor(f_cat)],
                  main=paste0(round(size_ha[1]),' ha'),ylim=c(50,100),xlim=c(10,200),
                  ylab='',xlab='edge density (m/ha)')
        )
        
    }
    if(!my_size %in% c(600,1740)){
        
        with(plot_df,
             plot(edge_density,100-median_percent_loss,col=my_cols[as.factor(f_cat)],pch=my_pch[as.factor(f_cat)],
                  main=paste0(round(size_ha[1]),' ha'),ylim=c(50,100),xlim=c(10,200),
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
 dev.off()

