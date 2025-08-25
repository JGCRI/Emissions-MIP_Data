# Compare atlantic/pacific sensitivity

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(purrr)
library(ncdf4)

# Specify and navigate to the location of Emissions-MIP directory
emi_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Documents/GitHub/Emissions-MIP_data/Emissions-MIP_data/")
setwd(paste0(emi_dir))

## Load in data using same method from summary plot script
# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#c4c4c3", "#4477aa", "#228833", "#66ccee", "#ccbb44","#ee6677", "#aa3377")

model_colors <- c('CESM1' = cbPalette[1], 'GISS-E2.1' = cbPalette[2], 'CAM-ATRAS' = cbPalette[3], 'GEOS' = cbPalette[4], 'NorESM2' = cbPalette[5], 'GFDL-ESM4' = cbPalette[6], 'E3SM' = cbPalette[7])
model_symbols <- c("CESM1" = 15, "GISS-E2.1" = 15, "CAM-ATRAS" = 17,  "GEOS" = 17, "NorESM2" = 17, "GFDL-ESM4" = 19, "E3SM" = 15)

# ------------------------------------------------------------------------------
# Reads in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models <- excluded_models %>% drop_na() #gets rid of any empty spaces
#-------------------------------------------------------------------------------
#extracts data for each perturbation experiment from csv files
data_accumulation <- function(emi_dir, reg_name, exper){
  
  setwd(paste0(emi_dir,'/input/', reg_name,'/', exper, '/diff'))
  
  # Read in csv files and bind into single data frame
  target_filename <- list.files(getwd(), "*.csv")
  regional_data <- rbind(map(target_filename, read.csv))
  regional_data <- lapply(regional_data, function(x) {x["unit"] <- NULL; x})
  regional_data <- bind_rows(regional_data)
  
  # Extract model from file names (fifth segment) and bind to experiment data frame
  models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[5])
  rep_models <- rep(models, each = 5) # five years
  regional_data$model <- rep_models
  
  # Convert SO2 volume mixing ratio to mass mixing ratio by multiplying by molar
  # mass of SO2 and dividing by molar mass of air, invert sign of forcing variables
  # to be consistent with convention (i.e. positive value denotes a heating effect),
  # then take the average over all years for each variable and calculate std dev
  regional_data_summary <- regional_data %>%
    dplyr::group_by(variable, model) %>%
    within(value <- ifelse(variable == "so2", 64.066 / 28.96, 1) * value) %>%
    within(value <- ifelse(variable %in% c("rlut", "rsut", "rlutcs", "rsutcs"), -1, 1) * value) %>%
    within(value <- ifelse(variable == "dms", 62.13 / 28.96, 1) * value) %>%
    # Convert from NH4HSO4 to SO4 mass
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "E3SM", 96/115, 1) * value) %>%
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "CESM", 96/115, 1) * value) %>%
    # Convert from H2SO4 to SO4 mass
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "NorESM2", 96/98, 1) * value) %>%
    dplyr::summarise(regional_data = mean(value), regional_data_sd = sd(value))
  
  return(regional_data_summary)
}

#read in data for each region
NH_atlantic_atl <- data_accumulation(emi_dir,'NH-atlantic',"shp-atl-shift") %>% 
  mutate(experiment="shp_atl_shift") %>% select(variable,model,experiment,everything())
NH_pacific_atl <- data_accumulation(emi_dir,'NH-pacific',"shp-atl-shift") %>% 
  mutate(experiment="shp_atl_shift") %>% select(variable,model,experiment,everything())
NH_atlantic_atl1950 <- data_accumulation(emi_dir,'NH-atlantic',"shp-atl-shift-1950") %>% 
  mutate(experiment="shp_atl_shift_1950") %>% select(variable,model,experiment,everything())
NH_pacific_atl1950 <- data_accumulation(emi_dir,'NH-pacific',"shp-atl-shift-1950") %>% 
  mutate(experiment="shp_atl_shift_1950") %>% select(variable,model,experiment,everything())

NH_atlantic_30 <- data_accumulation(emi_dir,'NH-atlantic',"shp-30p-red") %>% 
  mutate(experiment="shp_30p_red") %>% select(variable,model,experiment,everything())
NH_pacific_30 <- data_accumulation(emi_dir,'NH-pacific',"shp-30p-red") %>% 
  mutate(experiment="shp_30p_red") %>% select(variable,model,experiment,everything())
NH_atlantic_60 <- data_accumulation(emi_dir,'NH-atlantic',"shp-60p-red") %>% 
  mutate(experiment="shp_60p_red") %>% select(variable,model,experiment,everything())
NH_pacific_60 <- data_accumulation(emi_dir,'NH-pacific',"shp-60p-red") %>% 
  mutate(experiment="shp_60p_red") %>% select(variable,model,experiment,everything())
NH_atlantic_60_1950 <- data_accumulation(emi_dir,'NH-atlantic',"shp-60p-red-1950") %>% 
  mutate(experiment="shp_60p_red_1950") %>% select(variable,model,experiment,everything())
NH_pacific_60_1950 <- data_accumulation(emi_dir,'NH-pacific',"shp-60p-red-1950") %>% 
  mutate(experiment="shp_60p_red_1950") %>% select(variable,model,experiment,everything())

# Bind together by region
#NH_atlantic <- bind_rows(list(NH_atlantic_atl,NH_atlantic_atl1950,NH_atlantic_30,NH_atlantic_60,NH_atlantic_60_1950))
NH_atlantic <- bind_rows(list(NH_atlantic_atl,NH_atlantic_atl1950))
#NH_atlantic <- bind_rows(list(NH_atlantic_30,NH_atlantic_60,NH_atlantic_60_1950))
#NH_pacific <- bind_rows(list(NH_pacific_atl,NH_pacific_atl1950,NH_pacific_30,NH_pacific_60,NH_pacific_60_1950))
NH_pacific <- bind_rows(list(NH_pacific_atl,NH_pacific_atl1950))
#NH_pacific <- bind_rows(list(NH_pacific_30,NH_pacific_60,NH_pacific_60_1950))

#rename the mean and standard deviation results columns in each data frame
NH_atlantic <- rename(NH_atlantic, NH_atlantic = regional_data)
NH_pacific <- rename(NH_pacific, NH_pacific = regional_data)

NH_atlantic <- rename(NH_atlantic, NH_atlantic_sd = regional_data_sd)
NH_pacific <- rename(NH_pacific, NH_pacific_sd = regional_data_sd)

# Bind data together
summary_data <- list(NH_atlantic, NH_pacific) %>% reduce(left_join, by = c("variable", "model","experiment"))

# Correct model names
summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"
summary_data$model[which(summary_data$model == "CAM5")] <- "CAM-ATRAS"
summary_data$model[which(summary_data$model == "GFDL")] <- "GFDL-ESM4"
summary_data$model[which(summary_data$model == "GISS")] <- "GISS-E2.1"

# Change to long format
summary_long_exp <- summary_data %>%
  gather(region, value, -c(model, variable, experiment, NH_atlantic_sd, NH_pacific_sd)) %>%
  select(variable, model, experiment, region, value) %>%
  drop_na()

summary_long_sd <- summary_data %>%
  gather(region, sd, -c(model, variable, experiment, NH_atlantic, NH_pacific)) %>%
  select(variable, model, experiment, region, sd) %>%
  drop_na()

summary_long_sd$region <- gsub("_sd", "", summary_long_sd$region)

summary_long <- dplyr::left_join(summary_long_exp, summary_long_sd)

#runs through each excluded model pair and filters them out of summary_long
if(nrow(excluded_models) != 0) { #only runs if the data frame is not empty
  for (val in 1:nrow(excluded_models)) {
    
    summary_long <- filter(summary_long, experiment != excluded_models$Scenario[val] | model != excluded_models$Model[val] | variable != excluded_models$Variable[val])
  }
}

# Multiply by areas of ocean basins to get values over entire basin
get_basin_area <- function(grid_file) {
  nc_file <- nc_open(paste0(emi_dir,"input/",grid_file))
  area_array <- ncvar_get(nc_file,"cell_area")
  nc_close(nc_file)
  total_area <- sum(area_array,na.rm=TRUE)
  return(total_area)
}

atlantic_area <- get_basin_area("gridarea_NH-atlantic.nc")
pacific_area <- get_basin_area("gridarea_NH-pacific.nc")

summary_long$value[summary_long$region=="NH_atlantic"] <- summary_long$value[summary_long$region=="NH_atlantic"]*atlantic_area
summary_long$value[summary_long$region=="NH_pacific"] <- summary_long$value[summary_long$region=="NH_pacific"]*pacific_area

summary_long$sd[summary_long$region=="NH_atlantic"] <- summary_long$sd[summary_long$region=="NH_atlantic"]*atlantic_area
summary_long$sd[summary_long$region=="NH_pacific"] <- summary_long$sd[summary_long$region=="NH_pacific"]*pacific_area

# Generate plots
title_font <- 9.5
axis_font <- 9
axis_title_font <- 9

#creates a function that filters species out of a database and multiplies by ocean basin area
filter_species <- function(database, species){
  species <- dplyr::filter(database, variable == species) %>%
    #dplyr::mutate(experiment=exper) %>%
    dplyr::select(variable,model,region,experiment,everything())
  return(species)
}

#filters each species from summary_long
rlut <- filter_species(summary_long, "rlut")
rsut <- filter_species(summary_long, "rsut")
rlutcs <- filter_species(summary_long, "rlutcs")
rsutcs <- filter_species(summary_long, "rsutcs")
dryso2 <- filter_species(summary_long, "dryso2")
wetso2 <- filter_species(summary_long, "wetso2")
dryso4 <- filter_species(summary_long, "dryso4")
wetso4 <- filter_species(summary_long, "wetso4")
loadso4 <- filter_species(summary_long, "loadso4")
loadso2 <- filter_species(summary_long, "loadso2")

# Combined variables
net_rad <- dplyr::left_join(rlut, rsut, by = c("model", "region","experiment"))
net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, experiment, value, sd))
net_rad <- dplyr::mutate(net_rad,variable="net_rad") %>% 
  select(variable,everything()) %>%
  group_by(variable)

net_rad_cs <- dplyr::left_join(rlut, rsut, by = c("model", "region","experiment"))
net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, experiment, value, sd))
net_rad_cs <- dplyr::mutate(net_rad_cs,variable="net_rad_cs") %>% 
  select(variable,everything()) %>%
  group_by(variable)

imp_cld <- dplyr::left_join(rlut, rsut, by = c("model", "region","experiment"))
imp_cld <- dplyr::mutate(imp_cld, value = value.x + value.y) %>%
  dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
  dplyr::select(c(model, region, experiment, value, sd))
imp_cld <- dplyr::mutate(imp_cld,variable="imp_cld") %>% 
  select(variable,everything()) %>%
  group_by(variable)

# Get differences between Atlantic and Pacific (absolute values, magnitude of change)
atl_pac_difference <- function(variable) {
  variable %>%
    pivot_wider(names_from = region, values_from = c(value,sd)) %>%
    mutate(
      value_diff = abs(value_NH_atlantic) - abs(value_NH_pacific),
      sd_diff = sqrt(sd_NH_atlantic^2 + sd_NH_pacific^2)
    ) %>%
    select(variable,model,experiment,value_diff,sd_diff)
}

rsut_diff <- atl_pac_difference(rsut)
rlut_diff <- atl_pac_difference(rlut)
rsutcs_diff <- atl_pac_difference(rsutcs)
rlutcs_diff <- atl_pac_difference(rlutcs)
dryso2_diff <- atl_pac_difference(dryso2)
wetso2_diff <- atl_pac_difference(wetso2)
dryso4_diff <- atl_pac_difference(dryso4)
wetso4_diff <- atl_pac_difference(wetso4)
loadso4_diff <- atl_pac_difference(loadso4)
loadso2_diff <- atl_pac_difference(loadso2)
net_rad_diff <- atl_pac_difference(net_rad)
net_rad_cs_diff <- atl_pac_difference(net_rad_cs)
imp_cld_diff <- atl_pac_difference(imp_cld)

# Plot differences
plot_species <- function(variable, title, units, model_colors, model_symbols, ylimit=c(NA,NA)){
  species <- variable
  species_plot <- ggplot(species, aes(x = experiment, y = value_diff, color = model, shape = model))+
    theme_bw()+
    labs(title=paste0(title,' \n|Atlantic| - |Pacific|'), y=units) +
    theme(plot.title = element_text(hjust = 0.5, size = title_font),
          axis.text = element_text(size = axis_font),
          axis.title = element_text(size = axis_title_font),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank()) +
    scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value_diff))-max(abs(species$sd_diff)), max(abs(species$value_diff))+max(abs(species$sd_diff))))+
    scale_colour_manual(values = model_colors) +
    scale_shape_manual(values = model_symbols) +
    geom_point(position=position_dodge(width=0.4), size = 1.5) +
    geom_errorbar(aes(ymin=value_diff-sd_diff, ymax=value_diff+sd_diff), width=0.2, position=position_dodge(0.4), show.legend = F) +
    geom_hline(yintercept = 0) +
    if(missing(ylimit)) {
      
    } else {
      coord_cartesian(xlim=NULL,ylim=c(ylimit[1],ylimit[2]))
    }
  
  return(species_plot)
}

rsut_diff_plot <- plot_species(rsut_diff, "upwelling shortwave flux at TOA -",expression(Delta*~rsut~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
rlut_diff_plot <- plot_species(rlut_diff, "upwelling longwave flux at TOA -",expression(Delta*~rlut~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
rsutcs_diff_plot <- plot_species(rsutcs_diff, "upwelling clear-sky shortwave flux at TOA -",expression(Delta*~rsutcs~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
rlutcs_diff_plot <- plot_species(rlutcs_diff, "upwelling clear-sky longwave flux at TOA -",expression(Delta*~rlutcs~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
net_rad_diff_plot <- plot_species(net_rad_diff, "net radiative flux at TOA -",expression(Delta*~rlut~+~rsut~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
net_rad_cs_diff_plot <- plot_species(net_rad_cs_diff, "clear-sky net radiative flux at TOA -",expression(Delta*~rlutcs~+~rsutcs~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
imp_cld_diff_plot <- plot_species(imp_cld_diff, "implied cloud response at TOA -",expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs~(total~W)),model_colors,model_symbols,c(-2.2e13,2.2e13))
dryso2_diff_plot <- plot_species(dryso2_diff, "dry deposition rate of SO2 -",expression(Delta*~dryso2~(kg~s^-1)),model_colors,model_symbols,c(-30,30))
wetso2_diff_plot <- plot_species(wetso2_diff, "wet deposition rate of SO2 -",expression(Delta*~wetso2~(kg~s^-1)),model_colors,model_symbols,c(-30,30))
dryso4_diff_plot <- plot_species(dryso4_diff, "dry deposition rate of SO4 -",expression(Delta*~dryso4~(kg~s^-1)),model_colors,model_symbols,c(-30,30))
wetso4_diff_plot <- plot_species(wetso4_diff, "wet deposition rate of SO4 -",expression(Delta*~wetso4~(kg~s^-1)),model_colors,model_symbols,c(-30,30))
loadso2_diff_plot <- plot_species(loadso2_diff, "SO2 column burden -",expression(Delta*~loadso2~(kg)),model_colors,model_symbols,c(-4e6,4e6))
loadso4_diff_plot <- plot_species(loadso4_diff, "SO4 column burden -",expression(Delta*~loadso4~(kg)),model_colors,model_symbols,c(-4e6,4e6))

# Arrange plots
# Function from stack exchange to generate a shared legend
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom",
                                     legend.title = element_blank(),
                                     legend.text = element_text(size = 9,
                                                                margin = margin(r = 10, unit = "pt"))))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - 1.5 * lheight, lheight), # the "1.5" adds room for title
    top = textGrob("|Atlantic| - |Pacific| -- ocean basin totals", gp = gpar(fontsize = 12)))
}

forcing_plots <- grid_arrange_shared_legend(rlut_diff_plot,
                                            rsut_diff_plot,
                                            net_rad_diff_plot,
                                            imp_cld_diff_plot,
                                            rlutcs_diff_plot,
                                            rsutcs_diff_plot)

deposition_plots <- grid_arrange_shared_legend(dryso2_diff_plot,
                                               wetso2_diff_plot,
                                               dryso4_diff_plot,
                                               wetso4_diff_plot)

burden_plots <- grid_arrange_shared_legend(loadso2_diff_plot,
                                           loadso4_diff_plot)

# Print plots
setwd("../../../../output/comparison_plots")
pdf("atl-pac_comparison_atl-scenarios_diff.pdf", height = 11, width = 8.5, paper = "letter")

grid.draw(forcing_plots)
grid.newpage()
grid.draw(deposition_plots)
grid.newpage()
grid.draw(burden_plots)
dev.off()