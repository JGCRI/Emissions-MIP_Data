# 30p-60p comparison

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ivmte)

# Run after running summary_plots_diff.R with "global" region
global_data <- summary_long

global_30p <- filter(global_data,experiment=="shp_30p_red")
global_60p <- filter(global_data,experiment=="shp_60p_red")

## Get ratio of 30p difference to 60p difference
df_ratio <- left_join(global_30p,global_60p,by=c("variable","model"))
# Correct srfdms to dms
df_ratio$variable[which(df_ratio$variable == "srfdms")] <- "dms"
# Calculate ratio
#df_ratio <- drop_na(df_ratio) # can't divide by zero
df_ratio$ratio <- df_ratio$value.x/df_ratio$value.y
df_ratio$ratio.sd <- df_ratio$ratio*sqrt((df_ratio$sd.x/df_ratio$value.x)^2+(df_ratio$sd.y/df_ratio$value.y)^2)

# Get combined radiative flux variables
# Net radiative flux (rsut + rlut)
net_rad <- df_ratio %>% filter(variable=="rsut")
net_rad$variable[which(net_rad$variable == "rsut")] <- "net_rad"
net_rad$value.x <- net_rad$value.x + (df_ratio %>% filter(variable=="rlut"))$value.x
net_rad$sd.x <- sqrt(net_rad$sd.x^2 + (df_ratio %>% filter(variable=="rlut"))$sd.x^2)
net_rad$value.y <- net_rad$value.y + (df_ratio %>% filter(variable=="rlut"))$value.y
net_rad$sd.y <- sqrt(net_rad$sd.y^2 + (df_ratio %>% filter(variable=="rlut"))$sd.y^2)
net_rad$ratio <- net_rad$value.x/net_rad$value.y
net_rad$ratio.sd <- net_rad$ratio*sqrt((net_rad$sd.x/net_rad$value.x)^2+(net_rad$sd.y/net_rad$value.y)^2)

# Net clear-sky radiative flux (rsutcs + rlutcs)
net_rad_cs <- df_ratio %>% filter(variable=="rsutcs")
net_rad_cs$variable[which(net_rad_cs$variable == "rsutcs")] <- "net_rad_cs"
net_rad_cs$value.x <- net_rad_cs$value.x + (df_ratio %>% filter(variable=="rlutcs"))$value.x
net_rad_cs$sd.x <- sqrt(net_rad_cs$sd.x^2 + (df_ratio %>% filter(variable=="rlutcs"))$sd.x^2)
net_rad_cs$value.y <- net_rad_cs$value.y + (df_ratio %>% filter(variable=="rlutcs"))$value.y
net_rad_cs$sd.y <- sqrt(net_rad_cs$sd.y^2 + (df_ratio %>% filter(variable=="rlutcs"))$sd.y^2)
net_rad_cs$ratio <- net_rad_cs$value.x/net_rad_cs$value.y
net_rad_cs$ratio.sd <- net_rad_cs$ratio*sqrt((net_rad_cs$sd.x/net_rad_cs$value.x)^2+(net_rad_cs$sd.y/net_rad_cs$value.y)^2)

# Implied cloud response (net_rad - net_rad_cs)
imp_cld <- net_rad
imp_cld$variable[which(imp_cld$variable == "net_rad")] <- "imp_cld"
imp_cld$value.x <- net_rad$value.x - net_rad_cs$value.x
imp_cld$sd.x <- sqrt(net_rad$sd.x^2 + net_rad_cs$sd.x^2)
imp_cld$value.y <- net_rad$value.y - net_rad_cs$value.y
imp_cld$sd.y <- sqrt(net_rad$sd.y^2 + net_rad_cs$sd.y^2)
imp_cld$ratio <- imp_cld$value.x/imp_cld$value.y
imp_cld$ratio.sd <- imp_cld$ratio*sqrt((imp_cld$sd.x/imp_cld$value.x)^2+(imp_cld$sd.y/imp_cld$value.y)^2)

# Add combined variables to df_ratio
df_ratio <- bind_rows(df_ratio,net_rad,net_rad_cs,imp_cld)

# Ratio = 30p value / 60p value. If linear, ratio should be ~0.5

plot_ratios <- function(var) {
  df <- filter(df_ratio,variable==var)
  plot <- ggplot(df,aes(x=model,y=ratio,color=model,shape=model)) +
    geom_point(position=position_dodge(width=0.4), size = 1.5) +
    geom_errorbar(aes(ymin=ratio-ratio.sd,ymax=ratio+ratio.sd),
                  position=position_dodge(width=0.4)) +
    geom_hline(yintercept = 0.5) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_symbols) +
    labs(title=paste0(var),
         x="",y="shp-30p-red/shp-60p-red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  return(plot)
}

plot_scatter <- function(var) {
  df <- filter(df_ratio,variable==var)
  plot <- ggplot(df,aes(x=value.x,y=value.y))
  
  # Adjust error bar height/width
  yrange <- layer_scales(plot)$y$range$range[2]-layer_scales(plot)$y$range$range[1]
  xrange <- layer_scales(plot)$x$range$range[2]-layer_scales(plot)$x$range$range[1]
  
  ywidth <- 0#yrange/20
  xheight <- 0#xrange/20
  
  # Replot with adjusted error bar height/width
  plot <- ggplot(df,aes(x=value.x,y=value.y,color=model,shape=model)) +
    geom_point(size=1.5) +
    geom_errorbar(aes(ymin=value.y-sd.y,ymax=value.y+sd.y,width=ywidth)) +
    geom_errorbarh(aes(xmin=value.x-sd.x,xmax=value.x+sd.x,height=xheight)) +
    geom_abline(intercept=0,slope=2) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_symbols) +
    labs(title=paste0(var),
         x="shp-30p-red",y="shp-60p-red") +
    theme_bw()
  
  return(plot)
}

# Get ratio plots
cl_plot <- plot_ratios("cl")
clivi_plot <- plot_ratios("clivi")
clt_plot <- plot_ratios("clt")
cltc_plot <- plot_ratios("cltc")
dms_plot <- plot_ratios("dms")
drybc_plot <- plot_ratios("drybc")
dryso2_plot <- plot_ratios("dryso2")
dryso4_plot <- plot_ratios("dryso4")
emibc_plot <- plot_ratios("emibc")
emiso2_plot <- plot_ratios("emiso2")
loadbc_plot <- plot_ratios("loadbc")
loadso2_plot <- plot_ratios("loadso2")
loadso4_plot <- plot_ratios("loadso4")
mmrbc_plot <- plot_ratios("mmrbc")
mmrso4_plot <- plot_ratios("mmrso4")
od550aer_plot <- plot_ratios("od550aer")
rlut_plot <- plot_ratios("rlut")
rlutcs_plot <- plot_ratios("rlutcs")
rsut_plot <- plot_ratios("rsut")
rsutcs_plot <- plot_ratios("rsutcs")
rsdt_plot <- plot_ratios("rsdt")
so2_plot <- plot_ratios("so2")
wetbc_plot <- plot_ratios("wetbc")
wetso2_plot <- plot_ratios("wetso2")
wetso4_plot <- plot_ratios("wetso4")
net_rad_plot <- plot_ratios("net_rad")
net_rad_cs_plot <- plot_ratios("net_rad_cs")
imp_cld_plot <- plot_ratios("imp_cld")

# Get scatter plots
cl_scatter <- plot_scatter("cl")
clivi_scatter <- plot_scatter("clivi")
clt_scatter <- plot_scatter("clt")
cltc_scatter <- plot_scatter("cltc")
dms_scatter <- plot_scatter("dms")
drybc_scatter <- plot_scatter("drybc")
dryso2_scatter <- plot_scatter("dryso2")
dryso4_scatter <- plot_scatter("dryso4")
emibc_scatter <- plot_scatter("emibc")
emiso2_scatter <- plot_scatter("emiso2")
loadbc_scatter <- plot_scatter("loadbc")
loadso2_scatter <- plot_scatter("loadso2")
loadso4_scatter <- plot_scatter("loadso4")
mmrbc_scatter <- plot_scatter("mmrbc")
mmrso4_scatter <- plot_scatter("mmrso4")
od550aer_scatter <- plot_scatter("od550aer")
rlut_scatter <- plot_scatter("rlut")
rlutcs_scatter <- plot_scatter("rlutcs")
rsut_scatter <- plot_scatter("rsut")
rsutcs_scatter <- plot_scatter("rsutcs")
rsdt_scatter <- plot_scatter("rsdt")
so2_scatter <- plot_scatter("so2")
wetbc_scatter <- plot_scatter("wetbc")
wetso2_scatter <- plot_scatter("wetso2")
wetso4_scatter <- plot_scatter("wetso4")
net_rad_scatter <- plot_scatter("net_rad")
net_rad_cs_scatter <- plot_scatter("net_rad_cs")
imp_cld_scatter <- plot_scatter("imp_cld")


# Arrange plots
emissions_plot <- grid_arrange_shared_legend(emibc_plot,
                                             emiso2_plot,
                                             mmrbc_plot,
                                             mmrso4_plot,
                                             so2_plot,
                                             dms_plot)

forcing_plot <- grid_arrange_shared_legend(rlut_plot,
                                           rsut_plot,
                                           net_rad_plot,
                                           rlutcs_plot,
                                           rsutcs_plot,
                                           net_rad_cs_plot,
                                           imp_cld_plot)

cloud_plot <- grid_arrange_shared_legend(od550aer_plot,
                                         clt_plot,
                                         cltc_plot,
                                         cl_plot,
                                         clivi_plot)

deposition_plot <- grid_arrange_shared_legend(drybc_plot,
                                              wetbc_plot,
                                              #tot_bc_plot,
                                              dryso2_plot,
                                              wetso2_plot,
                                              dryso4_plot,
                                              wetso4_plot)
                                              #tot_s_plot)

column_plot <- grid_arrange_shared_legend(loadbc_plot,
                                          loadso2_plot,
                                          loadso4_plot)
                                          #so4_lifetime_plot,
                                          #so2_lifetime_plot)

# Print to pdf
# Output plots as pdf
setwd(paste0("C:/Users/done231/OneDrive - PNNL/Desktop/Phase1b_data/output/"))
pdf(paste0('30p-60p_ratio_comparison.pdf'), height = 11, width = 8.5, paper = "letter")

grid.draw(emissions_plot)
grid.newpage()
grid.draw(forcing_plot)
grid.newpage()
grid.draw(cloud_plot)
grid.newpage()
grid.draw(deposition_plot)
grid.newpage()
grid.draw(column_plot)
dev.off()

# Arrange scatter plots
emissions_scatter <- grid_arrange_shared_legend(emibc_scatter,
                                                emiso2_scatter,
                                                mmrbc_scatter,
                                                mmrso4_scatter,
                                                so2_scatter,
                                                dms_scatter)

forcing_scatter <- grid_arrange_shared_legend(rlut_scatter,
                                              rsut_scatter,
                                              net_rad_scatter,
                                              rlutcs_scatter,
                                              rsutcs_scatter,
                                              net_rad_cs_scatter,
                                              imp_cld_scatter)

cloud_scatter <- grid_arrange_shared_legend(od550aer_scatter,
                                            clt_scatter,
                                            cltc_scatter,
                                            cl_scatter,
                                            clivi_scatter)

deposition_scatter <- grid_arrange_shared_legend(drybc_scatter,
                                                 wetbc_scatter,
                                                 #tot_bc_scatter,
                                                 dryso2_scatter,
                                                 wetso2_scatter,
                                                 dryso4_scatter,
                                                 wetso4_scatter)

column_scatter <- grid_arrange_shared_legend(loadbc_scatter,
                                             loadso2_scatter,
                                             loadso4_scatter)

# Print to pdf
# Output plots as pdf
setwd(paste0("C:/Users/done231/OneDrive - PNNL/Desktop/Phase1b_data/output/"))
pdf(paste0('30p-60p_scatter_comparison.pdf'), height = 11, width = 8.5, paper = "letter")

grid.draw(emissions_scatter)
grid.newpage()
grid.draw(forcing_scatter)
grid.newpage()
grid.draw(cloud_scatter)
grid.newpage()
grid.draw(deposition_scatter)
grid.newpage()
grid.draw(column_scatter)
dev.off()