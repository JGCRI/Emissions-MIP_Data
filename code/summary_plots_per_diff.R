# ------------------------------------------------------------------------------
# Program Name: summary_plots_per_diff.R
# Authors: Hamza Ahsan
# Date Last Modified: November 15th, 2023
# Program Purpose: Produces summary plots of the percent difference between the
# perturbations and the reference case averaged over all years
# Input Files: ~Emissions-MIP/input/
# Output Files: ~Emissions-MIP/output/
# TODO:
# ------------------------------------------------------------------------------

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)

# Specify and navigate to the location of Emissions-MIP directory
emi_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Documents/GitHub/Emissions-MIP_Data/Emissions-MIP_Data/")
setwd(paste0(emi_dir))

# Specify what you are sorting by and either the region (i.e., global, land, sea, arctic, NH-land, NH-sea, SH-land, SH-sea) or experiment (i.e., bc-no-season, high-so4, no-so4, reference, so2-at-height, so2-no-season)
#The command line would look like: rscript <rscript>.r <"experiment" or "region"> <specific experiment or region you are sorting by>
sorting <- commandArgs(trailingOnly = TRUE) #pulling region from command line
sort_by <- sorting[1]
if (sort_by == "region"){region <- sorting[2]}
if (sort_by == "experiment"){exper <- sorting[2]}

sort_by <- "experiment"
#region <- "NH-indian"
exper <- "shp-60p-red"

# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
#cbPalette <- c("#999999", "#D55E00", "#117733", "#490092", "#F0E442","#0072B2", "#E69F00")
cbPalette <- c("#c4c4c3", "#4477aa", "#228833", "#66ccee", "#ccbb44","#ee6677", "#aa3377")

model_colors <- c('CESM1' = cbPalette[1], 'GISS-E2.1' = cbPalette[2], 'CAM-ATRAS' = cbPalette[3], 'GEOS' = cbPalette[4], 'NorESM2' = cbPalette[5], 'GFDL-ESM4' = cbPalette[6], 'E3SM' = cbPalette[7])
model_symbols <- c("CESM1" = 15, "GISS-E2.1" = 17, "CAM-ATRAS" = 17, "NorESM2" = 17, "GEOS" = 17, "GFDL-ESM4" = 19, "E3SM" = 15)

# ------------------------------------------------------------------------------
#reads in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models %>% drop_na() #gets rid of any empty spaces
#-----------------------------------------------------------------------------
#extracts data for each perturbation experiment from csv files
data_accumulation <- function(emi_dir, reg_name, exper){
  
  # Setup directory for bc-no-seas percent difference data
  setwd(paste0(emi_dir,'/input/', reg_name,'/', exper, '/per-diff'))
  
  # Read in csv files and bind into single data frame
  target_filename <- list.files(getwd(), "*.csv")
  regional_data <- rbind(map(target_filename, read.csv))
  regional_data <- lapply(regional_data, function(x) {x["unit"] <- NULL; x})
  regional_data <- bind_rows(regional_data)
  
  # Extract model from file names (fifth segment) and bind to experiment data frame
  models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[5])
  rep_models <- rep(models, each = 5) # five years
  regional_data$model <- rep_models
  
  # Take the average over all years for each variable
  regional_data_summary <- regional_data %>% dplyr::group_by(variable, model) %>%
    dplyr::summarise(regional_data = mean(value), regional_data_sd = sd(value))
  
  regional_data_summary <- filter(regional_data_summary,model!="GEOS")
  
  return(regional_data_summary)
}


if (sort_by == "region"){
  
  shp_30p_red_summary <- data_accumulation(emi_dir,region,"shp-30p-red")
  shp_60p_red_summary <- data_accumulation(emi_dir,region,"shp-60p-red")
  shp_60p_red_1950_summary <- data_accumulation(emi_dir,region,"shp-60p-red-1950")
  shp_atl_shift_summary <- data_accumulation(emi_dir,region,"shp-atl-shift")
  shp_atl_shift_1950_summary <- data_accumulation(emi_dir,region,"shp-atl-shift-1950")
  shp_ind_shift_summary <- data_accumulation(emi_dir,region,"shp-ind-shift")
  shp_ind_shift_1950_summary <- data_accumulation(emi_dir,region,"shp-ind-shift-1950")
  
  shp_30p_red_summary <- rename(shp_30p_red_summary, shp_30p_red = regional_data)
  shp_60p_red_summary <- rename(shp_60p_red_summary, shp_60p_red = regional_data)
  shp_60p_red_1950_summary <- rename(shp_60p_red_1950_summary, shp_60p_red_1950 = regional_data)
  shp_atl_shift_summary <- rename(shp_atl_shift_summary, shp_atl_shift = regional_data)
  shp_atl_shift_1950_summary <- rename(shp_atl_shift_1950_summary, shp_atl_shift_1950 = regional_data)
  shp_ind_shift_summary <- rename(shp_ind_shift_summary, shp_ind_shift = regional_data)
  shp_ind_shift_1950_summary <- rename(shp_ind_shift_1950_summary, shp_ind_shift_1950 = regional_data)
  
  shp_30p_red_summary <- rename(shp_30p_red_summary, shp_30p_red_sd = regional_data_sd)
  shp_60p_red_summary <- rename(shp_60p_red_summary, shp_60p_red_sd = regional_data_sd)
  shp_60p_red_1950_summary <- rename(shp_60p_red_1950_summary, shp_60p_red_1950_sd = regional_data_sd)
  shp_atl_shift_summary <- rename(shp_atl_shift_summary, shp_atl_shift_sd = regional_data_sd)
  shp_atl_shift_1950_summary <- rename(shp_atl_shift_1950_summary, shp_atl_shift_1950_sd = regional_data_sd)
  shp_ind_shift_summary <- rename(shp_ind_shift_summary, shp_ind_shift_sd = regional_data_sd)
  shp_ind_shift_1950_summary <- rename(shp_ind_shift_1950_summary, shp_ind_shift_1950_sd = regional_data_sd)
  
  # Bind data together
  summary_data <- list(shp_30p_red_summary, shp_60p_red_summary, shp_60p_red_1950_summary, shp_atl_shift_summary, shp_atl_shift_1950_summary, shp_ind_shift_summary, shp_ind_shift_1950_summary) %>% reduce(left_join, by = c("variable", "model"))
  
  # Correct model names
  summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"
  summary_data$model[which(summary_data$model == "GISS")] <- "GISS-E2.1"
  summary_data$model[which(summary_data$model == "CAM5")] <- "CAM-ATRAS"
  summary_data$model[which(summary_data$model == "GFDL")] <- "GFDL-ESM4"
  
  # Change to long format
  summary_long_exp <- summary_data %>%
    gather(experiment, value, -c(model, variable, shp_30p_red_sd, shp_60p_red_sd, shp_60p_red_1950_sd, shp_atl_shift_sd, shp_ind_shift_sd, shp_atl_shift_1950_sd, shp_ind_shift_1950_sd)) %>%
    select(variable, model, experiment, value) %>%
    drop_na()
  
  summary_long_sd <- summary_data %>%
    gather(experiment, sd, -c(model, variable, shp_30p_red, shp_60p_red, shp_60p_red_1950, shp_atl_shift, shp_atl_shift_1950, shp_ind_shift, shp_ind_shift_1950)) %>%
    select(variable, model, experiment, sd) %>%
    drop_na()
  
  summary_long_sd$experiment <- gsub("_sd", "", summary_long_sd$experiment)
  
  summary_long <- dplyr::left_join(summary_long_exp, summary_long_sd)
  
  #runs through each excluded model pair and filters them out of summary_long
  if(nrow(excluded_models) != 0) { #only runs if the data frame is not empty
    for (val in 1:nrow(excluded_models)) {
      summary_long <- filter(summary_long, experiment != excluded_models$Scenario[val] | model != excluded_models$Model[val] | variable != excluded_models$Variable[val])
    }
  }
  
}

if (sort_by == "experiment"){
  #read in data for each region
  arctic <- data_accumulation(emi_dir,'arctic',exper)
  global <- data_accumulation(emi_dir,'global',exper)
  land <- data_accumulation(emi_dir,'land',exper)
  NH_atlantic <- data_accumulation(emi_dir,'NH-atlantic',exper)
  NH_indian <- data_accumulation(emi_dir,'NH-indian',exper)
  NH_land <- data_accumulation(emi_dir,'NH-land',exper)
  NH_pacific <- data_accumulation(emi_dir,'NH-pacific',exper)
  NH_sea <- data_accumulation(emi_dir,'NH-sea',exper)
  sea <- data_accumulation(emi_dir,'sea',exper)
  SH_land <- data_accumulation(emi_dir,'SH-land',exper)
  SH_sea <- data_accumulation(emi_dir,'SH-sea',exper)
  
  #rename the mean and standard deviation results columns in each data frame
  arctic <- rename(arctic, arctic = regional_data)
  global <- rename(global, global = regional_data)
  land <- rename(land, land = regional_data)
  NH_atlantic <- rename(NH_atlantic, NH_atlantic = regional_data)
  NH_indian <- rename(NH_indian, NH_indian = regional_data)
  NH_land <- rename(NH_land, NH_land = regional_data)
  NH_pacific <- rename(NH_pacific, NH_pacific = regional_data)
  NH_sea <- rename(NH_sea, NH_sea = regional_data)
  sea <- rename(sea, sea = regional_data)
  SH_land <- rename(SH_land, SH_land = regional_data)
  SH_sea <- rename(SH_sea, SH_sea = regional_data)
  
  arctic <- rename(arctic, arctic_sd = regional_data_sd)
  global <- rename(global, global_sd = regional_data_sd)
  land <- rename(land, land_sd = regional_data_sd)
  NH_atlantic <- rename(NH_atlantic, NH_atlantic_sd = regional_data_sd)
  NH_indian <- rename(NH_indian, NH_indian_sd = regional_data_sd)
  NH_land <- rename(NH_land, NH_land_sd = regional_data_sd)
  NH_pacific <- rename(NH_pacific, NH_pacific_sd = regional_data_sd)
  NH_sea <- rename(NH_sea, NH_sea_sd = regional_data_sd)
  sea <- rename(sea, sea_sd = regional_data_sd)
  SH_land <- rename(SH_land, SH_land_sd = regional_data_sd)
  SH_sea <- rename(SH_sea, SH_sea_sd = regional_data_sd)
  
  # Bind data together
  summary_data <- list(arctic, global, land, NH_atlantic, NH_indian, NH_land, NH_pacific, NH_sea, sea, SH_land, SH_sea) %>% reduce(left_join, by = c("variable", "model"))
  
  # Correct model names
  summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"
  summary_data$model[which(summary_data$model == "GISS")] <- "GISS-E2.1"
  summary_data$model[which(summary_data$model == "CAM5")] <- "CAM-ATRAS"
  summary_data$model[which(summary_data$model == "GFDL")] <- "GFDL-ESM4"
  
  # Change to long format
  summary_long_exp <- summary_data %>%
    gather(region, value, -c(model, variable, arctic_sd, global_sd, land_sd, NH_atlantic_sd, NH_indian_sd, NH_land_sd, NH_pacific_sd, NH_sea_sd, sea_sd, SH_land_sd, SH_sea_sd)) %>%
    select(variable, model, region, value) %>%
    drop_na()
  
  summary_long_sd <- summary_data %>%
    gather(region, sd, -c(model, variable, arctic, global, land, NH_atlantic, NH_indian, NH_land, NH_pacific, NH_sea, sea, SH_land, SH_sea)) %>%
    select(variable, model, region, sd) %>%
    drop_na()
  
  summary_long_sd$region <- gsub("_sd", "", summary_long_sd$region)
  
  summary_long <- dplyr::left_join(summary_long_exp, summary_long_sd)
}

# Generate plots
title_font <- 9.5
axis_font <- 9
axis_title_font <- 9

#creates a function that filters species out of a database
filter_species <- function(database, species){
  species <- dplyr::filter(database, variable == species)
  return(species)
}

#filters each species from summary_long
emibc <- filter_species(summary_long, "emibc")
emiso2 <- filter_species(summary_long, "emiso2")
mmrbc <- filter_species(summary_long, "mmrbc")
mmrso4 <- filter_species(summary_long, "mmrso4")
so2 <- filter_species(summary_long, "so2")
rlut <- filter_species(summary_long, "rlut")
rsut <- filter_species(summary_long, "rsut")
rsdt <- filter_species(summary_long, "rsdt")
rlutcs <- filter_species(summary_long, "rlutcs")
rsutcs <- filter_species(summary_long, "rsutcs")
drybc <- filter_species(summary_long, "drybc")
wetbc <- filter_species(summary_long, "wetbc")
dryso2 <- filter_species(summary_long, "dryso2")
wetso2 <- filter_species(summary_long, "wetso2")
dryso4 <- filter_species(summary_long, "dryso4")
wetso4 <- filter_species(summary_long, "wetso4")
od550aer <- filter_species(summary_long, "od550aer")
clt <- filter_species(summary_long, "clt")
cltc <- filter_species(summary_long, "cltc")
cl <- filter_species(summary_long, "cl")
clivi <- filter_species(summary_long, "clivi")
dms <- filter_species(summary_long, "dms")
loadso4 <- filter_species(summary_long, "loadso4")
loadbc <- filter_species(summary_long, "loadbc")
loadso2 <- filter_species(summary_long, "loadso2")

# add E3SM's srfdms to dms
srfdms <- filter_species(summary_long, "srfdms")
dms <- rbind(dms, srfdms) # append srfdms to dms
dms[dms=="srfdms"] <- "dms" #replace srfdms variable name with dms

#Creates a function that creates plots for the data based on each species
if (sort_by == "region"){
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
    species <- variable
    species_plot <- ggplot(species, aes(x = experiment, y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=paste0(title,' - ', region_or_exper), y="Percent") +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_y_continuous(labels = function(x) paste0(format(x,digits=2), "%"), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd)))) +
      scale_colour_manual(values = model_colors) +
      scale_shape_manual(values = model_symbols) +
      geom_point( position=position_dodge(width = 0.4), size = 1.5) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F) +
      if(missing(ylimit)) {
        
      } else {
        ylim(ylimit[1],ylimit[2])
      }
    
    return(species_plot)
  }
  
  #creates plots based on each species using the plot_species function
  emibc_plot <- plot_species(emibc, region, value, 'surface flux of BC', expression(Delta*~emibc), region, model_colors, model_symbols)
  emiso2_plot <- plot_species(emiso2, region, value, 'surface flux of SO2', expression(Delta*~emiso2), region, model_colors, model_symbols)
  mmrbc_plot <- plot_species(mmrbc, region, value, 'surface concentration of BC', expression(Delta*~mmrbc), region, model_colors, model_symbols)
  mmrso4_plot <- plot_species(mmrso4, region, value, 'surface concentration of SO4', expression(Delta*~mmrso4), region, model_colors, model_symbols)
  so2_plot <- plot_species(so2, region, value, 'surface concentration of SO2', expression(Delta*~so2), region, model_colors, model_symbols)
  rlut_plot <- plot_species(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut), region, model_colors, model_symbols)
  rsut_plot <- plot_species(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut), region, model_colors, model_symbols)
  rsdt_plot <- plot_species(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt), region, model_colors, model_symbols)
  rlutcs_plot <- plot_species(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs), region, model_colors, model_symbols)
  rsutcs_plot <- plot_species(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs), region, model_colors, model_symbols)
  drybc_plot <- plot_species(drybc, region, value, 'dry deposition rate \n of BC', expression(Delta*~drybc), region, model_colors, model_symbols)
  wetbc_plot <- plot_species(wetbc, region, value, 'wet deposition rate \n of BC', expression(Delta*~wetbc), region, model_colors, model_symbols)
  dryso2_plot <- plot_species(dryso2, region, value, 'dry deposition rate \n of so2', expression(Delta*~dryso2), region, model_colors, model_symbols)
  wetso2_plot <- plot_species(wetso2, region, value, 'wet deposition rate \n of so2', expression(Delta*~wetso2), region, model_colors, model_symbols)
  dryso4_plot <- plot_species(dryso4, region, value, 'dry deposition rate \n of so4', expression(Delta*~dryso4), region, model_colors, model_symbols)
  wetso4_plot <- plot_species(wetso4, region, value, 'wet deposition rate \n of so4', expression(Delta*~wetso4), region, model_colors, model_symbols)
  od550aer_plot <-  plot_species(od550aer, region, value, 'ambient aerosol optical \n thickness at 550nm', expression(Delta*~od550aer), region, model_colors, model_symbols)
  clt_plot <- plot_species(clt, region, value, 'total cloud cover \n percentage', expression(Delta*~clt), region, model_colors, model_symbols)
  cltc_plot <- plot_species(cltc, region, value, 'convective cloud cover \n percentage', expression(Delta*~cltc), region, model_colors, model_symbols)
  cl_plot <- plot_species(cl, region, value, 'cloud cover \n percentage', "expression cl (%)", region, model_colors, model_symbols)
  clivi_plot <- plot_species(clivi, region, value, 'Ice water path', expression(Delta*~clivi~(kg~m^-2)), region, model_colors, model_symbols)
  dms_plot <- plot_species(dms, region, value, 'Dimethyl sulphide (DMS) mole fraction', expression(Delta*~dms~(mol~mol^-1)), region, model_colors, model_symbols)
  loadso4_plot  <- plot_species(loadso4, region, value, 'load \n of so4', expression(Delta*~loadso4~(kg~m^-2)), region, model_colors, model_symbols)
  loadbc_plot  <- plot_species(loadbc, region, value, 'load \n of bc', expression(Delta*~loadbc~(kg~m^-2)), region, model_colors, model_symbols)
  loadso2_plot <- plot_species(loadso2, region, value, 'load \n of so2', expression(Delta*~loadso2~(kg~m^-2)), region, model_colors, model_symbols)
}
if (sort_by == "experiment"){
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
    species <- variable
    species_plot <- ggplot(species, aes(x = region, y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=paste0(title,' - ', region_or_exper), y="Percent") +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_y_continuous(labels = function(x) paste0(format(x,digits=2), "%"), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd)))) +
      scale_colour_manual(values = model_colors) +
      scale_shape_manual(values = model_symbols) +
      geom_point( position=position_dodge(width = 0.4), size = 1.5) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F) +
      if(missing(ylimit)) {
        
      } else {
        ylim(ylimit[1],ylimit[2])
      }
    return(species_plot)
  }
  #creates plots based on each species using the plot_species function
  emibc_plot <- plot_species(emibc, region, value, 'surface flux of BC', expression(Delta*~emibc), exper, model_colors, model_symbols)
  emiso2_plot <- plot_species(emiso2, region, value, 'surface flux of SO2', expression(Delta*~emiso2), exper, model_colors, model_symbols)
  mmrbc_plot <- plot_species(mmrbc, region, value, 'surface concentration of BC', expression(Delta*~mmrbc), exper, model_colors, model_symbols)
  mmrso4_plot <- plot_species(mmrso4, region, value, 'surface concentration of SO4', expression(Delta*~mmrso4), exper, model_colors, model_symbols)
  so2_plot <- plot_species(so2, region, value, 'surface concentration of SO2', expression(Delta*~so2), exper, model_colors, model_symbols)
  rlut_plot <- plot_species(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut), exper, model_colors, model_symbols)
  rsut_plot <- plot_species(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut), exper, model_colors, model_symbols)
  rsdt_plot <- plot_species(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt), exper, model_colors, model_symbols)
  rlutcs_plot <- plot_species(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs), exper, model_colors, model_symbols)
  rsutcs_plot <- plot_species(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs), exper, model_colors, model_symbols)
  drybc_plot <- plot_species(drybc, region, value, 'dry deposition rate \n of BC', expression(Delta*~drybc), exper, model_colors, model_symbols)
  wetbc_plot <- plot_species(wetbc, region, value, 'wet deposition rate \n of BC', expression(Delta*~wetbc), exper, model_colors, model_symbols)
  dryso2_plot <- plot_species(dryso2, region, value, 'dry deposition rate \n of so2', expression(Delta*~dryso2), exper, model_colors, model_symbols)
  wetso2_plot <- plot_species(wetso2, region, value, 'wet deposition rate \n of so2', expression(Delta*~wetso2), exper, model_colors, model_symbols)
  dryso4_plot <- plot_species(dryso4, region, value, 'dry deposition rate \n of so4', expression(Delta*~dryso4), exper, model_colors, model_symbols)
  wetso4_plot <- plot_species(wetso4, region, value, 'wet deposition rate \n of so4', expression(Delta*~wetso4), exper, model_colors, model_symbols)
  od550aer_plot <-  plot_species(od550aer, region, value, 'ambient aerosol optical \n thickness at 550nm', expression(Delta*~od550aer), exper, model_colors, model_symbols)
  clt_plot <- plot_species(clt, region, value, 'total cloud cover \n percentage', expression(Delta*~clt), exper, model_colors, model_symbols)
  cltc_plot <- plot_species(cltc, region, value, 'convective cloud cover \n percentage', expression(Delta*~cltc), exper, model_colors, model_symbols)
  cl_plot <- plot_species(cl, region, value, 'cloud cover \n percentage', "expression cl (%)", exper, model_colors, model_symbols)
  clivi_plot <- plot_species(clivi, region, value, 'Ice water path', expression(Delta*~clivi~(kg~m^-2)), exper, model_colors, model_symbols)
  dms_plot <- plot_species(dms, region, value, 'Dimethyl sulphide (DMS) mole fraction', expression(Delta*~dms~(mol~mol^-1)), exper, model_colors, model_symbols)
  loadso4_plot  <- plot_species(loadso4, region, value, 'load \n of so4', expression(Delta*~loadso4~(kg~m^-2)), exper, model_colors, model_symbols,c(-10,1))
  loadbc_plot  <- plot_species(loadbc, region, value, 'load \n of bc', expression(Delta*~loadbc~(kg~m^-2)), exper, model_colors, model_symbols)
  loadso2_plot  <- plot_species(loadso2, region, value, 'load \n of so2', expression(Delta*~loadso2~(kg~m^-2)), exper, model_colors, model_symbols,c(-10,1))
  
}

# Define normal and clear-sky net radiative flux (sum of longwave and shortwave radiation)
if (sort_by == "region"){
  net_rad <- dplyr::left_join(rlut, rsut, by = c("model", "experiment"))
  net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  net_rad_cs <- dplyr::left_join(rlutcs, rsutcs, by = c("model", "experiment"))
  net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  #plots normal and clear sky net radiative flux using the plot_species function
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut), region, model_colors, model_symbols,c(-1,1))
  net_rad_cs_plot <- plot_species(net_rad_cs, region, value, 'clear-sky net radiative flux \n at TOA', expression(Delta*~rlutcs~+~rsutcs), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  net_rad <- dplyr::left_join(rlut, rsut, by = c("model", "region"))
  net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  net_rad_cs <- dplyr::left_join(rlutcs, rsutcs, by = c("model", "region"))
  net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  #plots normal and clear sky net radiative flux using the plot_species function
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut), exper, model_colors, model_symbols)
  net_rad_cs_plot <- plot_species(net_rad_cs, region, value, 'clear-sky net radiative flux \n at TOA', expression(Delta*~rlutcs~+~rsutcs), exper, model_colors, model_symbols)
}

# Define total BC deposition rate (sum of dry and wet BC )
if (sort_by == "region"){
  tot_bc <- dplyr::left_join(drybc, wetbc, by = c("model", "experiment"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  tot_bc_plot <- plot_species(tot_bc, region, value, 'total deposition rate \n of BC', expression(Delta*~drybc~+~wetbc), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  tot_bc <- dplyr::left_join(drybc, wetbc, by = c("model", "region"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  tot_bc_plot <- plot_species(tot_bc, region, value, 'total deposition rate \n of BC', expression(Delta*~drybc~+~wetbc), exper, model_colors, model_symbols)
}

# Define total S deposition rate (sum of dry and wet SO2/SO4 )
if (sort_by == "region"){
  dry_s <- dplyr::left_join(dryso2, dryso4, by = c("model", "experiment"))
  dry_s <- dplyr::mutate(dry_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  wet_s <- dplyr::left_join(wetso2, wetso4, by = c("model", "experiment"))
  wet_s <- dplyr::mutate(wet_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  tot_s <- dplyr::left_join(dry_s, wet_s, by = c("model", "experiment"))
  tot_s <- dplyr::mutate(tot_s, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  tot_s_plot <- plot_species(tot_s, region, value, 'total deposition rate \n of S', expression(Delta*(~dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  dry_s <- dplyr::left_join(dryso2, dryso4, by = c("model", "region"))
  dry_s <- dplyr::mutate(dry_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  wet_s <- dplyr::left_join(wetso2, wetso4, by = c("model", "region"))
  wet_s <- dplyr::mutate(wet_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  tot_s <- dplyr::left_join(dry_s, wet_s, by = c("model", "region"))
  tot_s <- dplyr::mutate(tot_s, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  tot_s_plot <- plot_species(tot_s, region, value, 'total deposition rate \n of S', expression(Delta*(~dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3), exper, model_colors, model_symbols)
}

if (sort_by == "region"){
  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("model", "experiment"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  imp_cld_plot <- plot_species(imp_cld, region, value, 'implied cloud response at TOA - \n', expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("model", "region"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  imp_cld_plot <- plot_species(imp_cld, region, value, 'implied cloud response at TOA - \n', expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs),exper, model_colors, model_symbols)
}

if (sort_by == "region"){
  # Define total SO4 (sum of dry and wet SO4)
  tot_so4 <- dplyr::left_join(dryso4, wetso4, by = c("model","experiment"))
  tot_so4 <- dplyr::mutate(tot_so4, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model,experiment,value,sd))
  
  tot_so4_plot <- plot_species(tot_so4, region, value, 'total SO4 - \n', expression(Delta*~dryso4~+~wetso4), region, model_colors, model_symbols)
}

if (sort_by == "region"){
  # define so4 lifetime (loadso4/tot_so4), convert from seconds to days
  so4_lifetime <- dplyr::left_join(loadso4, tot_so4, by = c("model","experiment"))
  so4_lifetime <- dplyr::mutate(so4_lifetime, value = (value.x/value.y)/86400) %>%
    dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
    dplyr::filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
    dplyr::filter(value >= 0) %>%
    dplyr::select(c(model,experiment,value,sd))
  
  so4_lifetime_plot <- plot_species(so4_lifetime, region, value, 'SO4 lifetime - \n ', expression(Delta*~loadso4/(dryso4~+~wetso4)~(days)), region, model_colors, model_symbols)
}

# so2 lifetime (loadso2/emiso2), convert from seconds to days
if (sort_by == "region"){
  # define so2 lifetime (loadso2/emiso2), convert seconds to days
  so2_lifetime <- dplyr::left_join(loadso2, emiso2, by = c("model","experiment"))
  so2_lifetime <- dplyr::mutate(so2_lifetime, value = (value.x/value.y)/86400) %>%
    dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
    dplyr::filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
    dplyr::filter(value >= 0) %>%
    dplyr::select(c(model,experiment,value,sd))
  
  so2_lifetime_plot <- plot_species(so2_lifetime, region, value, 'SO2 lifetime - \n', expression(Delta*~loadso2/emiso2~(days)), region, model_colors, model_symbols)
}

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
    top = textGrob("Summary - percent difference", gp = gpar(fontsize = 12)))
}

emissions_plot <- grid_arrange_shared_legend(emibc_plot,
                                             emiso2_plot,
                                             mmrbc_plot,
                                             mmrso4_plot,
                                             so2_plot,
                                             dms_plot)

forcing_plot <- grid_arrange_shared_legend(rlut_plot,
                                           rsut_plot,
                                           net_rad_plot,
                                           rsdt_plot,
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
                                              tot_bc_plot,
                                              dryso2_plot,
                                              wetso2_plot,
                                              dryso4_plot,
                                              wetso4_plot,
                                              tot_s_plot)

column_plot <- grid_arrange_shared_legend(loadbc_plot,
                                          loadso2_plot,
                                          loadso4_plot,
                                          so4_lifetime_plot,
                                          so2_lifetime_plot)

# Print plots
if (sort_by == 'region'){
  setwd(paste0('../../../../output/', region, '/summary'))
  
  pdf(paste0(region, '_summary_plots_per_diff-noGEOS.pdf'), height = 11, width = 8.5, paper = "letter")
}

if (sort_by == 'experiment'){
  setwd(paste0('../../../../output/', exper, '/summary'))
  
  pdf(paste0(exper, '_summary_plots_per_diff-noGEOS.pdf'), height = 11, width = 8.5, paper = "letter")
}

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