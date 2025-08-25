## Functions to be used in additional_paper_plots.R

#extracts data for each perturbation experiment from csv files
data_accumulation <- function(emi_dir, reg_name, exper){
  
  setwd(paste0(emi_dir,'/input/', reg_name,'/', exper, '/',diff_or_perdiff))
  
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

read_in_data <- function() {
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
  return(summary_long)
}

#creates a function that filters species out of a database
filter_species <- function(database, species){
  filtered_regions <- c("global","land","sea","NH_atlantic","NH_pacific","NH_indian","arctic")
  species <- dplyr::filter(database, variable == species, region %in% filtered_regions)
  return(species)
}

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
}

if (sort_by == "experiment"){
  plot_species <- function(variable, x, y, title, units, region_or_exper, model_colors, model_symbols, ylimit=c(NA,NA)){
    filtered_regions <- c("global","land","sea","NH_atlantic","NH_pacific","NH_indian","arctic")
    species <- variable
    species_plot <- ggplot(species, aes(x = factor(region, level=filtered_regions), y = value, color = model, shape = model))+
      theme_bw()+
      labs(title=title, y=units) +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(labels=c("global"="Global","land"="Land","sea"="Ocean","NH_atlantic"="NH Atlantic","NH_pacific"="NH Pacific","NH_indian"="NH Indian","arctic"="Arctic")) +
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
    #widths = 5,
    heights = unit.c(unit(1, "npc") - 1.5 * lheight, lheight)) # the "1.5" adds room for title
}
