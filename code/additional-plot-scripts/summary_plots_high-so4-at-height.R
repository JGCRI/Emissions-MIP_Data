# ------------------------------------------------------------------------------
# Program Name: summary_plots_diff_high-so4-at-height.R
# Authors: Hamza Ahsan
# Date Last Modified: July 1, 2024
# Program Purpose: Produces summary plots of the difference between the
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
emi_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Desktop/Phase1b_data")
setwd(paste0(emi_dir))

# Specify what you are sorting by and either the region (i.e., global, land, sea, arctic, NH-land, NH-sea, SH-land, SH-sea) or experiment (i.e., bc-no-season, high-so4, no-so4, reference, so2-at-height, so2-no-season)
#The command line would look like: rscript <rscript>.r <"experiment" or "region"> <specific experiment or region you are sorting by>
sorting <- commandArgs(trailingOnly = TRUE) #pulling region from command line
sort_by <- sorting[1]
if (sort_by == "region") {region <- sorting[2]}
if (sort_by == "experiment") {exper <- sorting[2]}

sort_by <- "region"
region <- "SH-land"
#exper <- "shp-60p-red"


# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#c4c4c3", "#4477aa", "#228833", "#66ccee", "#ccbb44","#ee6677", "#aa3377")

model_colors <- c('CESM1' = cbPalette[1], 'GISS-E2.1' = cbPalette[2], 'CAM-ATRAS' = cbPalette[3], 'GEOS' = cbPalette[4], 'NorESM2' = cbPalette[5], 'GFDL-ESM4' = cbPalette[6], 'E3SM' = cbPalette[7])
model_symbols <- c("CESM1" = 15, "GISS-E2.1" = 15, "CAM-ATRAS" = 17,  "GEOS" = 17, "NorESM2" = 17, "GFDL-ESM4" = 19, "E3SM" = 15)

# ------------------------------------------------------------------------------
# Reads in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models %>% drop_na() #gets rid of any empty spaces
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

#Creates a function to be used if accumulating data from a reference case
data_accumulation_reference <- function(emi_dir, reg_name, exper){
  
  setwd(paste0(emi_dir,'/input/', reg_name,'/', exper))
  
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
    dplyr::summarise(regional_data = mean(value), regional_data_sd = sd(value))
  
  return(regional_data_summary)
}

#-----------------------------------------------------------------------------
if (sort_by == "region"){
  high_so4_at_height_summary <- data_accumulation(emi_dir,region,"high-so4-at-height")
  high_so4_summary <- data_accumulation(emi_dir,region,"high-so4")
  so2_at_height_summary <- data_accumulation(emi_dir,region,"so2-at-height")
  
  high_so4_at_height_summary <- rename(high_so4_at_height_summary, high_so4_at_height = regional_data)
  high_so4_summary <- rename(high_so4_summary, high_so4 = regional_data)
  so2_at_height_summary <- rename(so2_at_height_summary, so2_at_height = regional_data)
  
  high_so4_at_height_summary <- rename(high_so4_at_height_summary, high_so4_at_height_sd = regional_data_sd)
  high_so4_summary <- rename(high_so4_summary, high_so4_sd = regional_data_sd)
  so2_at_height_summary <- rename(so2_at_height_summary, so2_at_height_sd = regional_data_sd)
  
  # Rename srfdms
  high_so4_at_height_summary$variable[which(high_so4_at_height_summary$variable == "srfdms")] <- "dms"
  
  # Bind data together
  summary_data <- list(high_so4_at_height_summary, high_so4_summary, so2_at_height_summary) %>% reduce(left_join, by = c("variable", "model"))
  
  # Correct model names for CESM, GISS, CAM-ATRAS, GFDL
  summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"
  summary_data$model[which(summary_data$model == "GISS")] <- "GISS-E2.1"
  summary_data$model[which(summary_data$model == "CAM5")] <- "CAM-ATRAS"
  summary_data$model[which(summary_data$model == "GFDL")] <- "GFDL-ESM4"
  
  # Change to long format
  summary_long_exp <- summary_data %>%
    gather(experiment, value, -c(model, variable, high_so4_at_height_sd, high_so4_sd, so2_at_height_sd)) %>%
    select(variable, model, experiment, value) %>%
    drop_na()
  
  summary_long_sd <- summary_data %>%
    gather(experiment, sd, -c(model, variable, high_so4_at_height, high_so4, so2_at_height)) %>%
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
  if (exper == "reference"){
    arctic <- data_accumulation_reference(emi_dir,'arctic',exper)
    global <- data_accumulation_reference(emi_dir,'global',exper)
    land <- data_accumulation_reference(emi_dir,'land',exper)
    NH_atlantic <- data_accumulation_reference(emi_dir,'NH-atlantic',exper)
    NH_indian <- data_accumulation_reference(emi_dir,'NH-indian',exper)
    NH_land <- data_accumulation_reference(emi_dir,'NH-land',exper)
    NH_pacific <- data_accumulation_reference(emi_dir,'NH-pacific',exper)
    NH_sea <- data_accumulation_reference(emi_dir,'NH-sea',exper)
    sea <- data_accumulation_reference(emi_dir,'sea',exper)
    SH_land <- data_accumulation_reference(emi_dir,'SH-land',exper)
    SH_sea <- data_accumulation_reference(emi_dir,'SH-sea',exper)
  }
  
  else {
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
  }
  
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
  summary_data$model[which(summary_data$model == "CAM5")] <- "CAM-ATRAS"
  summary_data$model[which(summary_data$model == "GFDL")] <- "GFDL-ESM4"
  summary_data$model[which(summary_data$model == "GISS")] <- "GISS-E2.1"
  
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
# srfdms <- filter_species(summary_long, "srfdms")
# dms <- rbind(dms, srfdms) # append srfdms to dms
# dms[dms=="srfdms"] <- "dms" #replace srfdms variable name with dms

#Creates a function that creates plots for the data based on each species
if (sort_by == "region"){
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
    species <- variable
    species_plot <- ggplot(variable, aes(x = experiment, y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=paste0(title,' - ', region_or_exper), y=units) +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd))))+
      scale_colour_manual(values = model_colors) +
      scale_shape_manual(values = model_symbols) +
      geom_point( position=position_dodge(width=0.4), size = 1.5) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F) +
      if(missing(ylimit)) {
        
      } else {
        ylim(ylimit[1],ylimit[2])
      }
    return(species_plot)
  }
  
  #creates plots based on each species using the plot_species function
  emibc_plot <- plot_species(emibc, region, value, 'surface flux of BC', expression(Delta*~emibc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  emiso2_plot <- plot_species(emiso2, region, value, 'surface flux of SO2', expression(Delta*~emiso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-10e-13,10e-13))
  mmrbc_plot <- plot_species(mmrbc, region, value, 'surface concentration of BC', expression(Delta*~mmrbc~(kg~kg-1)), region, model_colors, model_symbols)
  mmrso4_plot <- plot_species(mmrso4, region, value, 'surface concentration of SO4', expression(Delta*~mmrso4~(kg~kg-1)), region, model_colors, model_symbols)#,c(-2e-11,2e-11))
  so2_plot <- plot_species(so2, region, value, 'surface concentration of SO2', expression(Delta*~so2~(kg~kg-1)), region, model_colors, model_symbols)
  rlut_plot <- plot_species(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  rsut_plot <- plot_species(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  rsdt_plot <- plot_species(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  rlutcs_plot <- plot_species(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  rsutcs_plot <- plot_species(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  drybc_plot <- plot_species(drybc, region, value, 'dry deposition rate \n of BC', expression(Delta*~drybc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  wetbc_plot <- plot_species(wetbc, region, value, 'wet deposition rate \n of BC', expression(Delta*~wetbc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  dryso2_plot <- plot_species(dryso2, region, value, 'dry deposition rate \n of so2', expression(Delta*~dryso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-2e-13,2e-13))
  wetso2_plot <- plot_species(wetso2, region, value, 'wet deposition rate \n of so2', expression(Delta*~wetso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-1e-3,1e-3))
  dryso4_plot <- plot_species(dryso4, region, value, 'dry deposition rate \n of so4', expression(Delta*~dryso4~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-1e-13,1e-13))
  wetso4_plot <- plot_species(wetso4, region, value, 'wet deposition rate \n of so4', expression(Delta*~wetso4~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  od550aer_plot <-  plot_species(od550aer, region, value, 'ambient aerosol optical \n thickness at 550nm', expression(Delta*~od550aer), region, model_colors, model_symbols)
  clt_plot <- plot_species(clt, region, value, 'total cloud cover \n percentage',"expression clt (%)", region, model_colors, model_symbols)
  cltc_plot <- plot_species(cltc, region, value, 'convective cloud cover \n percentage', "expression cltc (%)", region, model_colors, model_symbols)#,c(-0.1,0.1))
  cl_plot <- plot_species(cl, region, value, 'cloud cover \n percentage', "expression cl (%)", region, model_colors, model_symbols)#,c(-0.2,0.2))
  clivi_plot <- plot_species(clivi, region, value, 'Ice water path', expression(Delta*~clivi~(kg~m^-2)), region, model_colors, model_symbols)
  dms_plot <- plot_species(dms, region, value, 'Dimethyl sulphide (DMS) mole fraction', expression(Delta*~dms~(mol~mol^-1)), region, model_colors, model_symbols)
  loadso4_plot  <- plot_species(loadso4, region, value, 'load \n of so4', expression(Delta*~loadso4~(kg~m^-2)), region, model_colors, model_symbols)
  loadbc_plot  <- plot_species(loadbc, region, value, 'load \n of bc', expression(Delta*~loadbc~(kg~m^-2)), region, model_colors, model_symbols)
  loadso2_plot <- plot_species(loadso2, region, value, 'load \n of so2', expression(Delta*~loadso2~(kg~m^-2)), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols){
    species <- variable
    species_plot <- ggplot(species, aes(x = region, y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=paste0(title,' - ', region_or_exper), y=units) +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd))))+
      scale_colour_manual(values = model_colors) +
      scale_shape_manual(values = model_symbols) +
      geom_point( position=position_dodge(width=0.4), size = 1.5) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F)
    
    return(species_plot)
  }
  #creates plots based on each species using the plot_species function
  emibc_plot <- plot_species(emibc, region, value, 'surface flux of BC', expression(Delta*~emibc~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  emiso2_plot <- plot_species(emiso2, region, value, 'surface flux of SO2', expression(Delta*~emiso2~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  mmrbc_plot <- plot_species(mmrbc, region, value, 'surface concentration of BC', expression(Delta*~mmrbc~(kg~kg-1)), exper, model_colors, model_symbols)
  mmrso4_plot <- plot_species(mmrso4, region, value, 'surface concentration of SO4', expression(Delta*~mmrso4~(kg~kg-1)), exper, model_colors, model_symbols)
  so2_plot <- plot_species(so2, region, value, 'surface concentration of SO2', expression(Delta*~so2~(kg~kg-1)), exper, model_colors, model_symbols)
  rlut_plot <- plot_species(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut~(W~m-2)), exper, model_colors, model_symbols)
  rsut_plot <- plot_species(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut~(W~m-2)), exper, model_colors, model_symbols)
  rsdt_plot <- plot_species(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt~(W~m-2)), exper, model_colors, model_symbols)
  rlutcs_plot <- plot_species(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs~(W~m-2)), exper, model_colors, model_symbols)
  rsutcs_plot <- plot_species(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs~(W~m-2)), exper, model_colors, model_symbols)
  drybc_plot <- plot_species(drybc, region, value, 'dry deposition rate \n of BC', expression(Delta*~drybc~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  wetbc_plot <- plot_species(wetbc, region, value, 'wet deposition rate \n of BC', expression(Delta*~wetbc~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  dryso2_plot <- plot_species(dryso2, region, value, 'dry deposition rate \n of so2', expression(Delta*~dryso2~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  wetso2_plot <- plot_species(wetso2, region, value, 'wet deposition rate \n of so2', expression(Delta*~wetso2~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  dryso4_plot <- plot_species(dryso4, region, value, 'dry deposition rate \n of so4', expression(Delta*~dryso4~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  wetso4_plot <- plot_species(wetso4, region, value, 'wet deposition rate \n of so4', expression(Delta*~wetso4~(kg~m^-2~s^-1)), exper, model_colors, model_symbols)
  od550aer_plot <-  plot_species(od550aer, region, value, 'ambient aerosol optical \n thickness at 550nm', expression(Delta*~od550aer), exper, model_colors, model_symbols)
  clt_plot <- plot_species(clt, region, value, 'total cloud cover \n percentage',"expression clt (%)", exper, model_colors, model_symbols)
  cltc_plot <- plot_species(cltc, region, value, 'convective cloud cover \n percentage', "expression cltc (%)", exper, model_colors, model_symbols)
  cl_plot <- plot_species(cl, region, value, 'cloud cover \n percentage', "expression cl (%)", exper, model_colors, model_symbols)
  clivi_plot <- plot_species(clivi, region, value, 'Ice water path', expression(Delta*~clivi~(kg~m^-2)), exper, model_colors, model_symbols)
  dms_plot <- plot_species(dms, region, value, 'Dimethyl sulphide (DMS) mole fraction', expression(Delta*~dms~(mol~mol^-1)), exper, model_colors, model_symbols)
  loadso4_plot  <- plot_species(loadso4, region, value, 'load \n of so4', expression(Delta*~loadso4~(kg~m^-2)), exper, model_colors, model_symbols)
  loadbc_plot  <- plot_species(loadbc, region, value, 'load \n of bc', expression(Delta*~loadbc~(kg~m^-2)), exper, model_colors, model_symbols)
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
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
  net_rad_cs_plot <- plot_species(net_rad_cs, region, value, 'clear-sky net radiative flux \n at TOA', expression(Delta*~rlutcs~+~rsutcs~(W~m-2)), region, model_colors, model_symbols,c(-1.25,1.25))
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
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut~(W~m-2)), exper, model_colors, model_symbols)
  net_rad_cs_plot <- plot_species(net_rad_cs, region, value, 'clear-sky net radiative flux \n at TOA', expression(Delta*~rlutcs~+~rsutcs~(W~m-2)), exper, model_colors, model_symbols)
}

# Define total BC deposition rate (sum of dry and wet BC )
if (sort_by == "region"){
  tot_bc <- dplyr::left_join(drybc, wetbc, by = c("model", "experiment"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  tot_bc_plot <- plot_species(tot_bc, region, value, 'total deposition rate \n of BC', expression(Delta*~drybc~+~wetbc~(kg~m-2~s-1)), region, model_colors, model_symbols)
}

if (sort_by == "experiment"){
  tot_bc <- dplyr::left_join(drybc, wetbc, by = c("model", "region"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  tot_bc_plot <- plot_species(tot_bc, region, value, 'total deposition rate \n of BC', expression(Delta*~drybc~+~wetbc~(kg~m-2~s-1)), exper, model_colors, model_symbols)
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
  
  tot_s_plot <- plot_species(tot_s, region, value, 'total deposition rate \n of S', expression(atop((dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3, (kg~m^-2~s^-1))), region, model_colors, model_symbols)
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
  
  tot_s_plot <- plot_species(tot_s, region, value, 'total deposition rate \n of S', expression(atop((dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3, (kg~m^-2~s^-1))), exper, model_colors, model_symbols)
}

if (sort_by == "region"){
  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("model", "experiment"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, experiment, value, sd))
  
  imp_cld_plot <- plot_species(imp_cld, region, value, 'implied cloud response at TOA - \n', expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)), region, model_colors, model_symbols,c(-1.25,1.25))
}

if (sort_by == "experiment"){
  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("model", "region"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(model, region, value, sd))
  
  imp_cld_plot <- plot_species(imp_cld, region, value, 'implied cloud response at TOA - \n', expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)),exper, model_colors, model_symbols)
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

# Plot with high_so4 + so2_at_height
if (sort_by == "region"){
  combo_plot <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
    species <- variable
    species_so4 <- filter(species,experiment=="high_so4")
    species_height <- filter(species,experiment=="so2_at_height")
    species_so4_at_height <- filter(species,experiment=="high_so4_at_height")
    
    species_sum <- dplyr::left_join(species_so4, species_height, by = c("variable", "model"))
    species_sum <- dplyr::mutate(species_sum, value = value.x + value.y) %>%
      dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
      dplyr::select(c(variable, model, value, sd))
    species_sum$experiment <- "combined-high-so4-and-height"
    
    species <- bind_rows(species_sum,species_so4_at_height)
    
    species_plot <- ggplot(species, aes(x = experiment, y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=paste0(title,' - ', region_or_exper), y=units) +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd))))+
      scale_colour_manual(values = model_colors) +
      scale_shape_manual(values = model_symbols) +
      geom_point( position=position_dodge(width=0.4), size = 1.5) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F) +
      if(missing(ylimit)) {
        
      } else {
        ylim(ylimit[1],ylimit[2])
      }
    return(species_plot)
  }
  
  #creates plots with high-so4 and so2-at-height summed together, compared to high-so4-at-height
  emibc_combo <- combo_plot(emibc, region, value, 'surface flux of BC', expression(Delta*~emibc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  emiso2_combo <- combo_plot(emiso2, region, value, 'surface flux of SO2', expression(Delta*~emiso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-10e-13,10e-13))
  mmrbc_combo <- combo_plot(mmrbc, region, value, 'surface concentration of BC', expression(Delta*~mmrbc~(kg~kg-1)), region, model_colors, model_symbols)
  mmrso4_combo <- combo_plot(mmrso4, region, value, 'surface concentration of SO4', expression(Delta*~mmrso4~(kg~kg-1)), region, model_colors, model_symbols)#,c(-2e-11,2e-11))
  so2_combo <- combo_plot(so2, region, value, 'surface concentration of SO2', expression(Delta*~so2~(kg~kg-1)), region, model_colors, model_symbols)
  rlut_combo <- combo_plot(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut~(W~m-2)), region, model_colors, model_symbols)
  rsut_combo <- combo_plot(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut~(W~m-2)), region, model_colors, model_symbols)
  rsdt_combo <- combo_plot(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt~(W~m-2)), region, model_colors, model_symbols)
  rlutcs_combo <- combo_plot(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs~(W~m-2)), region, model_colors, model_symbols)
  rsutcs_combo <- combo_plot(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs~(W~m-2)), region, model_colors, model_symbols)
  drybc_combo <- combo_plot(drybc, region, value, 'dry deposition rate \n of BC', expression(Delta*~drybc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  wetbc_combo <- combo_plot(wetbc, region, value, 'wet deposition rate \n of BC', expression(Delta*~wetbc~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  dryso2_combo <- combo_plot(dryso2, region, value, 'dry deposition rate \n of so2', expression(Delta*~dryso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-2e-13,2e-13))
  wetso2_combo <- combo_plot(wetso2, region, value, 'wet deposition rate \n of so2', expression(Delta*~wetso2~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-1e-3,1e-3))
  dryso4_combo <- combo_plot(dryso4, region, value, 'dry deposition rate \n of so4', expression(Delta*~dryso4~(kg~m^-2~s^-1)), region, model_colors, model_symbols)#,c(-1e-13,1e-13))
  wetso4_combo <- combo_plot(wetso4, region, value, 'wet deposition rate \n of so4', expression(Delta*~wetso4~(kg~m^-2~s^-1)), region, model_colors, model_symbols)
  od550aer_combo <-  combo_plot(od550aer, region, value, 'ambient aerosol optical \n thickness at 550nm', expression(Delta*~od550aer), region, model_colors, model_symbols)
  clt_combo <- combo_plot(clt, region, value, 'total cloud cover \n percentage',"expression clt (%)", region, model_colors, model_symbols)
  cltc_combo <- combo_plot(cltc, region, value, 'convective cloud cover \n percentage', "expression cltc (%)", region, model_colors, model_symbols)#,c(-0.1,0.1))
  cl_combo <- combo_plot(cl, region, value, 'cloud cover \n percentage', "expression cl (%)", region, model_colors, model_symbols)#,c(-0.2,0.2))
  clivi_combo <- combo_plot(clivi, region, value, 'Ice water path', expression(Delta*~clivi~(kg~m^-2)), region, model_colors, model_symbols)
  dms_combo <- combo_plot(dms, region, value, 'Dimethyl sulphide (DMS) mole fraction', expression(Delta*~dms~(mol~mol^-1)), region, model_colors, model_symbols)
  loadso4_combo  <- combo_plot(loadso4, region, value, 'load \n of so4', expression(Delta*~loadso4~(kg~m^-2)), region, model_colors, model_symbols)
  loadbc_combo  <- combo_plot(loadbc, region, value, 'load \n of bc', expression(Delta*~loadbc~(kg~m^-2)), region, model_colors, model_symbols)
  #loadso2_combo <- combo_plot(loadso2, region, value, 'load \n of so2', expression(Delta*~loadso2~(kg~m^-2)), region, model_colors, model_symbols)
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
    top = textGrob("Summary - absolute difference", gp = gpar(fontsize = 12)))
}

emissions_combo <- grid_arrange_shared_legend(emibc_combo,
                                              emiso2_combo,
                                              mmrbc_combo,
                                              mmrso4_combo,
                                              so2_combo,
                                              dms_combo)

forcing_combo <- grid_arrange_shared_legend(rlut_combo,
                                            rsut_combo,
                                            #net_rad_combo,
                                            #net_rad_cs_combo,
                                            #imp_cld_combo,
                                            rsdt_combo,
                                            rlutcs_combo,
                                            rsutcs_combo)

cloud_combo <- grid_arrange_shared_legend(od550aer_combo,
                                          clt_combo,
                                          cltc_combo,
                                          cl_combo,
                                          clivi_combo)

deposition_combo <- grid_arrange_shared_legend(drybc_combo,
                                               wetbc_combo,
                                               #tot_bc_combo,
                                               dryso2_combo,
                                               wetso2_combo,
                                               dryso4_combo,
                                               wetso4_combo)
                                               #tot_s_combo)

column_combo <- grid_arrange_shared_legend(loadbc_combo,
                                           #loadso2_combo,
                                           loadso4_combo)
                                           #so4_lifetime_combo,
                                           #so2_lifetime_combo)


# Print plots
if (sort_by == 'region'){
  setwd(paste0('../../../../output/', region, '/summary'))
  pdf(paste0(region, '_summary_plots_diff_high-so4-at-height.pdf'), height = 11, width = 8.5, paper = "letter")
}

if (sort_by == 'experiment'){
  setwd(paste0('../../../../output/', exper, '/summary'))
  
  pdf(paste0(exper, '_summary_plots_diff_high-so4-at-height.pdf'), height = 11, width = 8.5, paper = "letter")
}

grid.draw(emissions_combo)
grid.newpage()
grid.draw(forcing_combo)
grid.newpage()
grid.draw(cloud_combo)
grid.newpage()
grid.draw(deposition_combo)
grid.newpage()
grid.draw(column_combo)
#grid.newpage()
#grid.draw(other_plot)
#grid.newpage()
#grid.draw(other_plot_2)
dev.off()
