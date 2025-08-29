## Get average values for basins for shp-60p-red
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

sort_by <- "experiment"
exper <- "shp-60p-red"


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
  
  #regional_data_summary <- filter(regional_data_summary,model!="GEOS")
  
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
  summary_data <- list(shp_30p_red_summary, shp_60p_red_summary, shp_60p_red_1950_summary, shp_atl_shift_summary, shp_ind_shift_summary, shp_atl_shift_1950_summary, shp_ind_shift_1950_summary) %>% reduce(left_join, by = c("variable", "model"))
  
  # Correct model names for CESM, GISS, CAM-ATRAS, GFDL
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
    gather(experiment, sd, -c(model, variable, shp_30p_red, shp_60p_red, shp_60p_red_1950, shp_atl_shift, shp_ind_shift, shp_atl_shift_1950, shp_ind_shift_1950)) %>%
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
srfdms <- filter_species(summary_long, "srfdms")
dms <- rbind(dms, srfdms) # append srfdms to dms
dms[dms=="srfdms"] <- "dms" #replace srfdms variable name with dms

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
  rlut_plot <- plot_species(rlut, region, value, 'upwelling longwave flux \n at TOA', expression(Delta*~rlut~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
  rsut_plot <- plot_species(rsut, region, value, 'upwelling shortwave flux \n at TOA', expression(Delta*~rsut~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
  rsdt_plot <- plot_species(rsdt, region, value, 'incident shortwave flux \n at TOA', expression(Delta*~rsdt~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
  rlutcs_plot <- plot_species(rlutcs, region, value, 'upwelling clear-sky longwave \n flux at TOA', expression(Delta*~rlutcs~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
  rsutcs_plot <- plot_species(rsutcs, region, value, 'upwelling clear-sky shortwave \n flux at TOA', expression(Delta*~rsutcs~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
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
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
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
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F) +
      if(missing(ylimit)) {
        
      } else {
        ylim(ylimit[1],ylimit[2])
      }
    
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
  loadso2_plot  <- plot_species(loadso2, region, value, 'load \n of so2', expression(Delta*~loadso2~(kg~m^-2)), exper, model_colors, model_symbols)
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
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut~(W~m-2)), region, model_colors, model_symbols,c(-0.5,0.5))
  net_rad_cs_plot <- plot_species(net_rad_cs, region, value, 'clear-sky net radiative flux \n at TOA', expression(Delta*~rlutcs~+~rsutcs~(W~m-2)), region, model_colors, model_symbols,c(-0.25,0.25))
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
  net_rad_plot <- plot_species(net_rad, region, value, 'net radiative flux \n at TOA', expression(Delta*~rlut~+~rsut~(W~m-2)), exper, model_colors, model_symbols,c(-0.8,0.8))
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
  
  imp_cld_plot <- plot_species(imp_cld, region, value, 'implied cloud response at TOA - \n', expression(Delta*~rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)), region, model_colors, model_symbols,c(-0.25,0.25))
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

## Get average net rad values
net_rad_global <- net_rad %>% filter(region=="global")
net_rad_nh_atlantic <- net_rad %>% filter(region=="NH_atlantic")
net_rad_nh_pacific <- net_rad %>% filter(region=="NH_pacific")
net_rad_nh_indian <- net_rad %>% filter(region=="NH_indian")
net_rad_land <- net_rad %>% filter(region=="land")
net_rad_sea <- net_rad %>% filter(region=="sea")

# Average values
net_rad_global_avg <- mean(net_rad_global$value)
net_rad_nh_atlantic_avg <- mean(net_rad_nh_atlantic$value)
net_rad_nh_pacific_avg <- mean(net_rad_nh_pacific$value)
net_rad_nh_indian_avg <- mean(net_rad_nh_indian$value)
net_rad_land_avg <- mean(net_rad_land$value)
net_rad_sea_avg <- mean(net_rad_sea$value)

# Average SDs
net_rad_global_sd <- sqrt(sum(net_rad_global$sd^2))
net_rad_nh_atlantic_sd <- sqrt(sum(net_rad_nh_atlantic$sd^2))
net_rad_nh_pacific_sd <- sqrt(sum(net_rad_nh_pacific$sd^2))
net_rad_nh_indian_sd <- sqrt(sum(net_rad_nh_indian$sd^2))
net_rad_land_sd <- sqrt(sum(net_rad_land$sd^2))
net_rad_sea_sd <- sqrt(sum(net_rad_sea$sd^2))

# Combine into dataset
net_rad_basin_averages <- data.frame(
  region=c("global","NH_atlantic","NH_pacific","NH_indian","land","sea"),
  mean=c(net_rad_global_avg,net_rad_nh_atlantic_avg,net_rad_nh_pacific_avg,
         net_rad_nh_indian_avg,net_rad_land_avg,net_rad_sea_avg),
  sd=c(net_rad_global_sd,net_rad_nh_atlantic_sd,net_rad_nh_pacific_sd,
       net_rad_nh_indian_sd,net_rad_land_sd,net_rad_sea_sd)
)

write.csv(net_rad_basin_averages,"../../../../output/basin_averages/net_rad_basin_averages.csv")

## Same process but without GEOS (outlier)
## Get average net rad values
net_rad_global <- net_rad %>% filter(region=="global",model!="GEOS")
net_rad_nh_atlantic <- net_rad %>% filter(region=="NH_atlantic",model!="GEOS")
net_rad_nh_pacific <- net_rad %>% filter(region=="NH_pacific",model!="GEOS")
net_rad_nh_indian <- net_rad %>% filter(region=="NH_indian",model!="GEOS")
net_rad_land <- net_rad %>% filter(region=="land",model!="GEOS")
net_rad_sea <- net_rad %>% filter(region=="sea",model!="GEOS")

# Average values
net_rad_global_avg <- mean(net_rad_global$value)
net_rad_nh_atlantic_avg <- mean(net_rad_nh_atlantic$value)
net_rad_nh_pacific_avg <- mean(net_rad_nh_pacific$value)
net_rad_nh_indian_avg <- mean(net_rad_nh_indian$value)
net_rad_land_avg <- mean(net_rad_land$value)
net_rad_sea_avg <- mean(net_rad_sea$value)

# Average SDs
net_rad_global_sd <- sqrt(sum(net_rad_global$sd^2))
net_rad_nh_atlantic_sd <- sqrt(sum(net_rad_nh_atlantic$sd^2))
net_rad_nh_pacific_sd <- sqrt(sum(net_rad_nh_pacific$sd^2))
net_rad_nh_indian_sd <- sqrt(sum(net_rad_nh_indian$sd^2))
net_rad_land_sd <- sqrt(sum(net_rad_land$sd^2))
net_rad_sea_sd <- sqrt(sum(net_rad_sea$sd^2))

# Combine into dataset
net_rad_basin_averages <- data.frame(
  region=c("global","NH_atlantic","NH_pacific","NH_indian","land","sea"),
  mean=c(net_rad_global_avg,net_rad_nh_atlantic_avg,net_rad_nh_pacific_avg,
         net_rad_nh_indian_avg,net_rad_land_avg,net_rad_sea_avg),
  sd=c(net_rad_global_sd,net_rad_nh_atlantic_sd,net_rad_nh_pacific_sd,
       net_rad_nh_indian_sd,net_rad_land_sd,net_rad_sea_sd)
)

write.csv(net_rad_basin_averages,"../../../../output/basin_averages/net_rad_basin_averages-noGEOS.csv")