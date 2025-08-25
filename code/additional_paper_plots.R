## Additional plots for paper

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ncdf4)
library(cowplot)

# Specify and navigate to the location of Emissions-MIP directory
emi_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Documents/GitHub/Emissions-MIP_Data/Emissions-MIP_Data/")
setwd(paste0(emi_dir))

# Select sorting method (region or experiment)
sort_by <- "experiment"
exper <- "shp-60p-red"
diff_or_perdiff <- 'diff'

# Source functions from other file (must be after defining sort_by because of if statements)
source("code/plot_functions.R")

# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#c4c4c3", "#4477aa", "#228833", "#66ccee", "#ccbb44","#ee6677", "#aa3377")

model_colors <- c('CESM1' = cbPalette[1], 'GISS-E2.1' = cbPalette[2], 'CAM-ATRAS' = cbPalette[3], 'GEOS' = cbPalette[4], 'NorESM2' = cbPalette[5], 'GFDL-ESM4' = cbPalette[6], 'E3SM' = cbPalette[7])
model_symbols <- c("CESM1" = 15, "GISS-E2.1" = 15, "CAM-ATRAS" = 17,  "GEOS" = 17, "NorESM2" = 17, "GFDL-ESM4" = 19, "E3SM" = 15)

# Read in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models <- excluded_models %>% drop_na() #gets rid of any empty spaces

# Call function to read in data
summary_long <- read_in_data()

# Generate plots
title_font <- 13
axis_font <- 10
axis_title_font <- 12

#filters each species from summary_long
loadso4 <- filter_species(summary_long, "loadso4")
loadso2 <- filter_species(summary_long, "loadso2")
dryso2_diff <- filter_species(summary_long, "dryso2")
wetso2_diff <- filter_species(summary_long, "wetso2")
dryso4_diff <- filter_species(summary_long, "dryso4")
wetso4_diff <- filter_species(summary_long, "wetso4")
rsut <- filter_species(summary_long,"rsut")
rlut <- filter_species(summary_long,"rlut")
rsutcs <- filter_species(summary_long, "rsutcs")
rlutcs <- filter_species(summary_long, "rlutcs")

if (sort_by=="region") {
  #creates plots based on each species using the plot_species function
  loadso4_plot  <- plot_species(loadso4, region, value, 'd) SO4 Absolute Difference', expression(Delta*~loadso4~(kg~m^-2)), region, model_colors, model_symbols, ylimit=c(NA,0))
  loadso2_plot <- plot_species(loadso2, region, value, 'a) SO2 Absolute Difference', expression(Delta*~loadso2~(kg~m^-2)), region, model_colors, model_symbols, ylimit=c(NA,0))
}

if (sort_by=="experiment") {
  #creates plots based on each species using the plot_species function
  loadso4_plot  <- plot_species(loadso4, region, value, 'd) SO4 Absolute Difference', expression(Delta*~loadso4~(kg~m^-2)), exper, model_colors, model_symbols, ylimit=c(NA,0))
  loadso2_plot  <- plot_species(loadso2, region, value, 'a) SO2 Absolute Difference', expression(Delta*~loadso2~(kg~m^-2)), exper, model_colors, model_symbols, ylimit=c(NA,0))
}

# Save loadso2 and loadso4
loadso2_diff <- loadso2_plot
loadso4_diff <- loadso4_plot

loadso2_diff_df <- loadso2
loadso4_diff_df <- loadso4

## Generate loadso4/net rad linear regressions
net_rad <- dplyr::left_join(rlut, rsut, by = c("model", "region"))
net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, value, sd))

net_rad_cs <- dplyr::left_join(rlutcs, rsutcs, by = c("model", "region"))
net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, value, sd))

imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("model", "region"))
imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, value, sd))

loadso4_netrad_combined <- dplyr::left_join(loadso4,net_rad, by = c("model","region"))
loadso4_netrad_combined <- drop_na(loadso4_netrad_combined)
loadso4_netrad_combined <- loadso4_netrad_combined %>% filter(model!="GEOS")
loadso4_netrad_plot <- ggplot(loadso4_netrad_combined, aes(value.x,value.y,color=model)) +
  geom_point() +
  
  # trend lines
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="CESM1"),se=FALSE, linetype = "dashed") +
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="CAM-ATRAS"),se=FALSE, linetype = "dashed") +
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="E3SM"),se=FALSE, linetype = "dashed") +
  #geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="GEOS"),se=FALSE, linetype = "dashed") +
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="GFDL-ESM4"),se=FALSE, linetype = "dashed") +
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="GISS-E2.1"),se=FALSE, linetype = "dashed") +
  geom_smooth(method=lm,data=filter(loadso4_rsut_combined,model=="NorESM2"),se=FALSE, linetype = "dashed") +
  
  # aesthetics
  scale_color_manual(values = model_colors) +
  #facet_wrap(vars(model)) +
  geom_errorbar(aes(ymin=value.y-sd.y,ymax=value.y+sd.y),width=0) +
  geom_errorbarh(aes(xmin=value.x-sd.x,xmax=value.x+sd.x),height=0) +
  #xlim(-10000,10000) +
  #ylim(-1,0.6) +
  ylab(expression(Delta~net~radiative~flux~(W~m^-2))) +
  xlab(expression(Delta~SO4~column~burden~(kg~m^-2)))

# Get values in table
data <- loadso4_netrad_combined
data_name <- "loadso4_netrad_combined"

# Split data by model
data_CAMATRAS <- data[data$model=="CAM-ATRAS",]
data_CESM1 <- data[data$model=="CESM1",]
data_E3SM <- data[data$model=="E3SM",]
data_GFDL <- data[data$model=="GFDL-ESM4",]
data_GISS <- data[data$model=="GISS-E2.1",]
data_NorESM <- data[data$model=="NorESM2",]

# Linear regressions
model_CAMATRAS <- lm(value.y ~ value.x, data = data_CAMATRAS)
r_CAMATRAS <- summary(model_CAMATRAS)$r.squared
intercept_CAMATRAS <- model_CAMATRAS$coefficients[1]
slope_CAMATRAS <- model_CAMATRAS$coefficients[2]

model_CESM1 <- lm(value.y ~ value.x, data = data_CESM1)
r_CESM1 <- summary(model_CESM1)$r.squared
intercept_CESM1 <- model_CESM1$coefficients[1]
slope_CESM1 <- model_CESM1$coefficients[2]

model_E3SM <- lm(value.y ~ value.x, data = data_E3SM)
r_E3SM <- summary(model_E3SM)$r.squared
intercept_E3SM <- model_E3SM$coefficients[1]
slope_E3SM <- model_E3SM$coefficients[2]

model_GFDL <- lm(value.y ~ value.x, data = data_GFDL)
r_GFDL <- summary(model_GFDL)$r.squared
intercept_GFDL <- model_GFDL$coefficients[1]
slope_GFDL <- model_GFDL$coefficients[2]

model_GISS <- lm(value.y ~ value.x, data = data_GISS)
r_GISS <- summary(model_GISS)$r.squared
intercept_GISS <- model_GISS$coefficients[1]
slope_GISS <- model_GISS$coefficients[2]

model_NorESM <- lm(value.y ~ value.x, data = data_NorESM)
r_NorESM <- summary(model_NorESM)$r.squared
intercept_NorESM <- model_NorESM$coefficients[1]
slope_NorESM <- model_NorESM$coefficients[2]

# Put values in table
models <- c("CAM-ATRAS","CESM1","E3SM","GFDL-ESM4","GISS-E2.1","NorESM2")
slopes <- c(slope_CAMATRAS,slope_CESM1,slope_E3SM,slope_GFDL,slope_GISS,slope_NorESM) %>% as.numeric()
intercepts <- c(intercept_CAMATRAS,intercept_CESM1,intercept_E3SM,intercept_GFDL,intercept_GISS,intercept_NorESM) %>% as.numeric()
r_squared <- c(r_CAMATRAS,r_CESM1,r_E3SM,r_GFDL,r_GISS,r_NorESM)

table <- data.frame(models=models,
                    slopes=slopes,
                    intercepts=intercepts,
                    r_squared=r_squared)

# Note: R2 is much less reliable across regions with one experiment, than across experiments with one region
# Makes sense that different regions are impacted differently

if (sort_by=="experiment") {
  #creates plots based on each species using the plot_species function
  net_rad_filtered <- net_rad %>% filter(model!="GEOS")
  rsut_filtered <- rsut %>% filter(model!="GEOS")
  rlut_filtered <- rlut %>% filter(model!="GEOS")
  imp_cld_filtered <- imp_cld %>% filter(model!="GEOS")
  
  rlut_plot  <- plot_species(rlut_filtered, region, value, 'a) Upwelling Longwave Flux', expression(Delta*~rlut~(W~m^-2)), exper, model_colors, model_symbols, c(-0.75,0.75))
  rsut_plot  <- plot_species(rsut_filtered, region, value, 'b) Upwelling Shortwave Flux', expression(Delta*~rsut~(W~m^-2)), exper, model_colors, model_symbols, c(-0.75,0.75))
  net_rad_plot  <- plot_species(net_rad_filtered, region, value, 'c) Net Radiative Flux', expression(Delta*~net_rad~(W~m^-2)), exper, model_colors, model_symbols, c(-0.75,0.75))
  imp_cld_plot <- plot_species(imp_cld_filtered, region, value, 'd) Implied Cloud Response', expression(Delta*~imp_cld~(W~m^-2)), exper, model_colors, model_symbols, c(-0.75,0.75))
}

# Save so2/so4 plots together
flux_plots <- grid_arrange_shared_legend(rlut_plot,rsut_plot,net_rad_plot,imp_cld_plot,ncol=4)
png(paste0(emi_dir,"/output/flux_plots_shp-60p-red.png"))
grid.draw(flux_plots)
dev.off()

## Regenerate plots using per-diff
diff_or_perdiff <- 'per-diff'
setwd(paste0(emi_dir))
source("code/plot_functions.R")

# Call function to read in data
summary_long <- read_in_data()

#filters each species from summary_long
loadso4 <- filter_species(summary_long, "loadso4")
loadso2 <- filter_species(summary_long, "loadso2")

if (sort_by=="region") {
  #creates plots based on each species using the plot_species function
  loadso4_plot  <- plot_species(loadso4, region, value, 'f) SO4 Percent Difference', "%", region, model_colors, model_symbols, ylimit=c(NA,0))
  loadso2_plot <- plot_species(loadso2, region, value, 'c) SO2 Percent Difference', "%", region, model_colors, model_symbols, ylimit=c(NA,0))
}

if (sort_by=="experiment") {
  #creates plots based on each species using the plot_species function
  loadso4_plot  <- plot_species(loadso4, region, value, 'f) SO4 Percent Difference', "%", exper, model_colors, model_symbols, ylimit=c(NA,0))
  loadso2_plot  <- plot_species(loadso2, region, value, 'c) SO2 Percent Difference', "%", exper, model_colors, model_symbols, ylimit=c(NA,0))
}

# Save loadso2 and loadso4
loadso2_perdiff <- loadso2_plot
loadso4_perdiff <- loadso4_plot

# Save so2/so4 plots together
#sox_plots <- grid_arrange_shared_legend(loadso2_diff,loadso4_diff,loadso2_perdiff,loadso4_perdiff)
legend <- get_legend(loadso2_diff)
loadso2_diff2 <- loadso2_diff + theme(legend.position = "none")
loadso4_diff2 <- loadso4_diff + theme(legend.position = "none")
loadso2_perdiff2 <- loadso2_perdiff + theme(legend.position = "none")
loadso4_perdiff2 <- loadso4_perdiff + theme(legend.position = "none")
sox_plots <- plot_grid(loadso2_diff2,loadso4_diff2,loadso2_perdiff2,loadso4_perdiff2,ncol=2,align="hv",axis="tblr")
sox_plots_grid <- plot_grid(sox_plots,legend,ncol=2,rel_heights=c(4,5))
png(paste0(emi_dir,"/output/SOx_plots_shp-60p-red.png"))
grid.draw(sox_plots_grid)
dev.off()

## Generate load plots with absolute amounts (kg over entire basin)
global_area <- get_basin_area("gridarea_global.nc")
sea_area <- get_basin_area("gridarea_sea.nc")
land_area <- global_area - sea_area # no gridarea_land file?
atlantic_area <- get_basin_area("gridarea_NH-atlantic.nc")
pacific_area <- get_basin_area("gridarea_NH-pacific.nc")
indian_area <- get_basin_area("gridarea_NH-indian.nc")
arctic_area <- get_basin_area("gridarea_arctic.nc")

# Filter out unneeded regions
filtered_regions <- c("global","land","sea","NH_atlantic","NH_pacific","NH_indian","arctic")
loadso2 <- loadso2_diff_df %>% filter(region %in% filtered_regions)
loadso4 <- loadso4_diff_df %>% filter(region %in% filtered_regions)

# Multiply by area of basin and 1e-9 (kg -> Tg), to get units in Tg over basin
loadso2$value[loadso2$region=="global"] <- loadso2$value[loadso2$region=="global"]*global_area*1e-9
loadso2$value[loadso2$region=="sea"] <- loadso2$value[loadso2$region=="sea"]*sea_area*1e-9
loadso2$value[loadso2$region=="land"] <- loadso2$value[loadso2$region=="land"]*land_area*1e-9
loadso2$value[loadso2$region=="NH_atlantic"] <- loadso2$value[loadso2$region=="NH_atlantic"]*atlantic_area*1e-9
loadso2$value[loadso2$region=="NH_pacific"] <- loadso2$value[loadso2$region=="NH_pacific"]*pacific_area*1e-9
loadso2$value[loadso2$region=="NH_indian"] <- loadso2$value[loadso2$region=="NH_indian"]*indian_area*1e-9
loadso2$value[loadso2$region=="arctic"] <- loadso2$value[loadso2$region=="arctic"]*arctic_area*1e-9

loadso4$value[loadso4$region=="global"] <- loadso4$value[loadso4$region=="global"]*global_area*1e-9
loadso4$value[loadso4$region=="sea"] <- loadso4$value[loadso4$region=="sea"]*sea_area*1e-9
loadso4$value[loadso4$region=="land"] <- loadso4$value[loadso4$region=="land"]*land_area*1e-9
loadso4$value[loadso4$region=="NH_atlantic"] <- loadso4$value[loadso4$region=="NH_atlantic"]*atlantic_area*1e-9
loadso4$value[loadso4$region=="NH_pacific"] <- loadso4$value[loadso4$region=="NH_pacific"]*pacific_area*1e-9
loadso4$value[loadso4$region=="NH_indian"] <- loadso4$value[loadso4$region=="NH_indian"]*indian_area*1e-9
loadso4$value[loadso4$region=="arctic"] <- loadso4$value[loadso4$region=="arctic"]*arctic_area*1e-9

# Multiply standard deviation by constants
loadso2$sd[loadso2$region=="global"] <- loadso2$sd[loadso2$region=="global"]*global_area*1e-9
loadso2$sd[loadso2$region=="sea"] <- loadso2$sd[loadso2$region=="sea"]*sea_area*1e-9
loadso2$sd[loadso2$region=="land"] <- loadso2$sd[loadso2$region=="land"]*land_area*1e-9
loadso2$sd[loadso2$region=="NH_atlantic"] <- loadso2$sd[loadso2$region=="NH_atlantic"]*atlantic_area*1e-9
loadso2$sd[loadso2$region=="NH_pacific"] <- loadso2$sd[loadso2$region=="NH_pacific"]*pacific_area*1e-9
loadso2$sd[loadso2$region=="NH_indian"] <- loadso2$sd[loadso2$region=="NH_indian"]*indian_area*1e-9
loadso2$sd[loadso2$region=="arctic"] <- loadso2$sd[loadso2$region=="arctic"]*arctic_area*1e-9

loadso4$sd[loadso4$region=="global"] <- loadso4$sd[loadso4$region=="global"]*global_area*1e-9
loadso4$sd[loadso4$region=="sea"] <- loadso4$sd[loadso4$region=="sea"]*sea_area*1e-9
loadso4$sd[loadso4$region=="land"] <- loadso4$sd[loadso4$region=="land"]*land_area*1e-9
loadso4$sd[loadso4$region=="NH_atlantic"] <- loadso4$sd[loadso4$region=="NH_atlantic"]*atlantic_area*1e-9
loadso4$sd[loadso4$region=="NH_pacific"] <- loadso4$sd[loadso4$region=="NH_pacific"]*pacific_area*1e-9
loadso4$sd[loadso4$region=="NH_indian"] <- loadso4$sd[loadso4$region=="NH_indian"]*indian_area*1e-9
loadso4$sd[loadso4$region=="arctic"] <- loadso4$sd[loadso4$region=="arctic"]*arctic_area*1e-9

## Plot
loadso2_basin_plot  <- plot_species(loadso2, region, value, 'b) SO2 Basin Total', expression(Delta*~loadso2~(Tg)), exper, model_colors, model_symbols,ylimit=c(NA,0.005))
loadso4_basin_plot  <- plot_species(loadso4, region, value, 'e) SO4 Basin Total', expression(Delta*~loadso4~(Tg)), exper, model_colors, model_symbols,ylimit=c(NA,0.01))

# Save so2/so4 plots together
sox_basin_plots <- grid_arrange_shared_legend(loadso2_basin_plot,loadso4_basin_plot)
png(paste0(emi_dir,"/output/SOx_basin_plots_shp-60p-red.png"))
grid.draw(sox_basin_plots)
dev.off()

# All so2/so4 plots (absolute diff per m, absolute diff in total kg, % diff)
loadso2_basin_plot2 <- loadso2_basin_plot + theme(legend.position = "none")
loadso4_basin_plot2 <- loadso4_basin_plot + theme(legend.position = "none")
all_sox_plots <- plot_grid(loadso2_diff2,loadso2_basin_plot2,loadso2_perdiff2,loadso4_diff2,
                           loadso4_basin_plot2,loadso4_perdiff2,
                           ncol=3,align="hv",axis="tblr")
all_sox_plots_grid <- plot_grid(all_sox_plots,legend,ncol=3,rel_heights=c(4,5))

png(paste0(emi_dir,"/output/all_SOx_plots_shp-60p-red.png"))
grid.draw(all_sox_plots_grid)
dev.off()

## Total deposition over basins
dryso2 <- dryso2_diff
wetso2 <- wetso2_diff
dryso4 <- dryso4_diff
wetso4 <- wetso4_diff

# Multiply by area of basin and 1e-9 (kg -> Tg), to get units in Tg over basin
dryso2$value[dryso2$region=="global"] <- dryso2$value[dryso2$region=="global"]*global_area*31536000
dryso2$value[dryso2$region=="sea"] <- dryso2$value[dryso2$region=="sea"]*sea_area*31536000
dryso2$value[dryso2$region=="land"] <- dryso2$value[dryso2$region=="land"]*land_area*31536000
dryso2$value[dryso2$region=="NH_atlantic"] <- dryso2$value[dryso2$region=="NH_atlantic"]*atlantic_area*31536000
dryso2$value[dryso2$region=="NH_pacific"] <- dryso2$value[dryso2$region=="NH_pacific"]*pacific_area*31536000
dryso2$value[dryso2$region=="NH_indian"] <- dryso2$value[dryso2$region=="NH_indian"]*indian_area*31536000
dryso2$value[dryso2$region=="arctic"] <- dryso2$value[dryso2$region=="arctic"]*arctic_area*31536000

dryso4$value[dryso4$region=="global"] <- dryso4$value[dryso4$region=="global"]*global_area*31536000
dryso4$value[dryso4$region=="sea"] <- dryso4$value[dryso4$region=="sea"]*sea_area*31536000
dryso4$value[dryso4$region=="land"] <- dryso4$value[dryso4$region=="land"]*land_area*31536000
dryso4$value[dryso4$region=="NH_atlantic"] <- dryso4$value[dryso4$region=="NH_atlantic"]*atlantic_area*31536000
dryso4$value[dryso4$region=="NH_pacific"] <- dryso4$value[dryso4$region=="NH_pacific"]*pacific_area*31536000
dryso4$value[dryso4$region=="NH_indian"] <- dryso4$value[dryso4$region=="NH_indian"]*indian_area*31536000
dryso4$value[dryso4$region=="arctic"] <- dryso4$value[dryso4$region=="arctic"]*arctic_area*31536000

wetso2$value[wetso2$region=="global"] <- wetso2$value[wetso2$region=="global"]*global_area*31536000
wetso2$value[wetso2$region=="sea"] <- wetso2$value[wetso2$region=="sea"]*sea_area*31536000
wetso2$value[wetso2$region=="land"] <- wetso2$value[wetso2$region=="land"]*land_area*31536000
wetso2$value[wetso2$region=="NH_atlantic"] <- wetso2$value[wetso2$region=="NH_atlantic"]*atlantic_area*31536000
wetso2$value[wetso2$region=="NH_pacific"] <- wetso2$value[wetso2$region=="NH_pacific"]*pacific_area*31536000
wetso2$value[wetso2$region=="NH_indian"] <- wetso2$value[wetso2$region=="NH_indian"]*indian_area*31536000
wetso2$value[wetso2$region=="arctic"] <- wetso2$value[wetso2$region=="arctic"]*arctic_area*31536000

wetso4$value[wetso4$region=="global"] <- wetso4$value[wetso4$region=="global"]*global_area*31536000
wetso4$value[wetso4$region=="sea"] <- wetso4$value[wetso4$region=="sea"]*sea_area*31536000
wetso4$value[wetso4$region=="land"] <- wetso4$value[wetso4$region=="land"]*land_area*31536000
wetso4$value[wetso4$region=="NH_atlantic"] <- wetso4$value[wetso4$region=="NH_atlantic"]*atlantic_area*31536000
wetso4$value[wetso4$region=="NH_pacific"] <- wetso4$value[wetso4$region=="NH_pacific"]*pacific_area*31536000
wetso4$value[wetso4$region=="NH_indian"] <- wetso4$value[wetso4$region=="NH_indian"]*indian_area*31536000
wetso4$value[wetso4$region=="arctic"] <- wetso4$value[wetso4$region=="arctic"]*arctic_area*31536000

# Multiply standard deviation by constants
dryso2$sd[dryso2$region=="global"] <- dryso2$sd[dryso2$region=="global"]*global_area*31536000
dryso2$sd[dryso2$region=="sea"] <- dryso2$sd[dryso2$region=="sea"]*sea_area*31536000
dryso2$sd[dryso2$region=="land"] <- dryso2$sd[dryso2$region=="land"]*land_area*31536000
dryso2$sd[dryso2$region=="NH_atlantic"] <- dryso2$sd[dryso2$region=="NH_atlantic"]*atlantic_area*31536000
dryso2$sd[dryso2$region=="NH_pacific"] <- dryso2$sd[dryso2$region=="NH_pacific"]*pacific_area*31536000
dryso2$sd[dryso2$region=="NH_indian"] <- dryso2$sd[dryso2$region=="NH_indian"]*indian_area*31536000
dryso2$sd[dryso2$region=="arctic"] <- dryso2$sd[dryso2$region=="arctic"]*arctic_area*31536000

dryso4$sd[dryso4$region=="global"] <- dryso4$sd[dryso4$region=="global"]*global_area*31536000
dryso4$sd[dryso4$region=="sea"] <- dryso4$sd[dryso4$region=="sea"]*sea_area*31536000
dryso4$sd[dryso4$region=="land"] <- dryso4$sd[dryso4$region=="land"]*land_area*31536000
dryso4$sd[dryso4$region=="NH_atlantic"] <- dryso4$sd[dryso4$region=="NH_atlantic"]*atlantic_area*31536000
dryso4$sd[dryso4$region=="NH_pacific"] <- dryso4$sd[dryso4$region=="NH_pacific"]*pacific_area*31536000
dryso4$sd[dryso4$region=="NH_indian"] <- dryso4$sd[dryso4$region=="NH_indian"]*indian_area*31536000
dryso4$sd[dryso4$region=="arctic"] <- dryso4$sd[dryso4$region=="arctic"]*arctic_area*31536000

wetso2$sd[wetso2$region=="global"] <- wetso2$sd[wetso2$region=="global"]*global_area*31536000
wetso2$sd[wetso2$region=="sea"] <- wetso2$sd[wetso2$region=="sea"]*sea_area*31536000
wetso2$sd[wetso2$region=="land"] <- wetso2$sd[wetso2$region=="land"]*land_area*31536000
wetso2$sd[wetso2$region=="NH_atlantic"] <- wetso2$sd[wetso2$region=="NH_atlantic"]*atlantic_area*31536000
wetso2$sd[wetso2$region=="NH_pacific"] <- wetso2$sd[wetso2$region=="NH_pacific"]*pacific_area*31536000
wetso2$sd[wetso2$region=="NH_indian"] <- wetso2$sd[wetso2$region=="NH_indian"]*indian_area*31536000
wetso2$sd[wetso2$region=="arctic"] <- wetso2$sd[wetso2$region=="arctic"]*arctic_area*31536000

wetso4$sd[wetso4$region=="global"] <- wetso4$sd[wetso4$region=="global"]*global_area*31536000
wetso4$sd[wetso4$region=="sea"] <- wetso4$sd[wetso4$region=="sea"]*sea_area*31536000
wetso4$sd[wetso4$region=="land"] <- wetso4$sd[wetso4$region=="land"]*land_area*31536000
wetso4$sd[wetso4$region=="NH_atlantic"] <- wetso4$sd[wetso4$region=="NH_atlantic"]*atlantic_area*31536000
wetso4$sd[wetso4$region=="NH_pacific"] <- wetso4$sd[wetso4$region=="NH_pacific"]*pacific_area*31536000
wetso4$sd[wetso4$region=="NH_indian"] <- wetso4$sd[wetso4$region=="NH_indian"]*indian_area*31536000
wetso4$sd[wetso4$region=="arctic"] <- wetso4$sd[wetso4$region=="arctic"]*arctic_area*31536000

## Plot
wetso2 <- wetso2 %>% filter(model!="E3SM")
wetso2_basin_plot  <- plot_species(wetso2, region, value, 'SO2 Wet Deposition', expression(Delta*~wetso2~(kg)), exper, model_colors, model_symbols)
dryso2_basin_plot  <- plot_species(dryso2, region, value, 'SO2 Dry Deposition', expression(Delta*~dryso2~(kg)), exper, model_colors, model_symbols)
wetso4_basin_plot  <- plot_species(wetso4, region, value, 'SO4 Wet Deposition', expression(Delta*~wetso4~(kg)), exper, model_colors, model_symbols)
dryso4_basin_plot  <- plot_species(dryso4, region, value, 'SO4 Dry Deposition', expression(Delta*~dryso4~(kg)), exper, model_colors, model_symbols)

# Save so2/so4 plots together
deposition_basin_plots <- grid_arrange_shared_legend(dryso2_basin_plot,wetso2_basin_plot,
                                                     dryso4_basin_plot,wetso4_basin_plot)
png(paste0(emi_dir,"/output/deposition_basin_plots_shp-60p-red.png"))
grid.draw(deposition_basin_plots)
dev.off()

so2_deposition <- left_join(wetso2,dryso2,by=c("model","region"))
so2_deposition$total_deposition <- so2_deposition$value.x + so2_deposition$value.y
so4_deposition <- left_join(wetso4,dryso4,by=c("model","region"))
so4_deposition$total_deposition <- so4_deposition$value.x + so4_deposition$value.y
