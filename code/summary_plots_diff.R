# ------------------------------------------------------------------------------
# Program Name: summary_plots_diff.R
# Authors: Hamza Ahsan, Harrison Suchyta
# Date Last Modified: April 4, 2022
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

# Specify location of Emissions-MIP directory
emi_dir <- paste0('C:/Users/such559/Documents/Emissions-MIP_Data')

setwd(paste0(emi_dir))

#determines whether the script sorts by region or experiment if nothing is put into command line
sort_by <- 'region'

#determines which region or experiment is sorted by
region <- 'global'
exper <- 'high-SO4-at-height'

# Specify what you are sorting by and either the region (i.e., global, land, sea, arctic, NH-land, NH-sea, SH-land, SH-sea) or experiment (i.e., bc-no-season, high-so4, no-so4, reference, so2-at-height, so2-no-season)
#The command line would look like: rscript <rscript>.r <"experiment" or "region"> <specific experiment or region you are sorting by>
sorting <- commandArgs(trailingOnly = TRUE) #pulling region from command line
if (!is_empty(sorting)){sort_by <- sorting[1]}
if (!is_empty(sorting) & sort_by == "region"){region <- sorting[2]}
if (!is_empty(sorting) & sort_by == "experiment"){experiment <- sorting[2]}


#-------------------------------------------------------------------------------
#Read in variable, region, and experiment names and sort them into their own lists
var_master_list <- read.csv(file = paste0(emi_dir, '/input/var_master_list.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
master_vars <- var_master_list$vars
master_com_var <- var_master_list$com_vars
master_region <- var_master_list$reg
master_exp <- var_master_list$exper
#remove blank strings
master_vars <- master_vars[master_vars != ""]
master_com_var<- master_com_var[master_com_var != ""]
master_region <- master_region[master_region != ""]
master_exp <- master_exp[master_exp != ""]
#-------------------------------------------------------------------------------
#Read in csv responsible for organizing combined variable outputs
combined_vars <- read.csv(file=paste0(emi_dir,'/input/combined_variables.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
#-------------------------------------------------------------------------------
#Read in fixed data to determine whether the axes should be fixed or grouped
fixed_data <- read.csv(file = paste0(emi_dir, '/input/fixed_axes.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
#-------------------------------------------------------------------------------
#Read in the variable min/max data
variables <- read.csv(file = paste0(emi_dir, '/input/variables.csv'), fileEncoding="UTF-8-BOM")
rownames(variables) <- variables$Variable
variables <- subset(variables, select = -c(Variable))
#creates a list of all the variables as strings
list_of_variable_strings <- rownames(variables)
#-------------------------------------------------------------------------------
#Checks to see if data was entered into csv files correctly

#Checks if total variables are consistent between variables.csv and var_master_list.csv
#Add 2 for the wet and dry s variables that would not appear in the master lists
if(nrow(variables) != length(master_vars) + length(master_com_var) + 3){
  stop('discrepancy between variables.csv and var_master_list.csv files')
}
#checks if the number of combined vars in combined_variables.csv added to the
#master list of variables is equivalent to the number of vars in variables
if(nrow(variables) != length(master_vars) + ncol(combined_vars)){
  stop('missing combined variables in combined_variables.csv')
}

# ------------------------------------------------------------------------------
# Reads in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models %>% drop_na() #gets rid of any empty spaces
#-------------------------------------------------------------------------------
# Define default ggplot colors and associate with models (in case a plot is
# missing a model, the color scheme will remain consistent)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(10)
# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#920000",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#490092",
               "#117733")

model_colors <- c(CESM1 = cbPalette[1], E3SM = cbPalette[2], GISS = cbPalette[3])
                  #,
                  # CESM2 = cbPalette[4], MIROC = cbPalette[5], NorESM2 = cbPalette[6],
                  # GFDL = cbPalette[7], OsloCTM3 = cbPalette[8], UKESM = cbPalette[9],
                  # GEOS = cbPalette[10], CAM5 = cbPalette[11])

model_symbols <- c(CESM1 = 15, E3SM = 15, GISS = 17) #, CESM2 = 19, MIROC = 15,
                   #NorESM2 = 17, GFDL = 19, OsloCTM3 = 19, UKESM = 15, GEOS = 17,
                   #CAM5 = 17)

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
  rep_models <- rep(models, each = 5) # four years
  regional_data$model <- rep_models

  # Convert SO2 volume mixing ratio to mass mixing ratio by multiplying by molar
  # mass of SO2 and dividing by molar mass of air, invert sign of forcing variables
  # to be consistent with convention (i.e. positive value denotes a heating effect),
  # then take the average over all years for each variable and calculate std dev
  regional_data_summary <- regional_data %>%
    within(value <- ifelse(variable == "so2", 64.066 / 28.96, 1) * value) %>%
    within(value <- ifelse(variable %in% c("rlut", "rsut", "rlutcs", "rsutcs"), -1, 1) * value) %>%
    within(value <- ifelse(variable %in% c("wetbc", "wetso2", "wetso4") & model == "CESM2", -1, 1) * value) %>%
    within(value <- ifelse(variable == "dms", 62.13 / 28.96, 1) * value) %>%
    dplyr::group_by(variable, model) %>%
    dplyr::summarise(regional_data = mean(value), regional_data_sd = sd(value))

  return(regional_data_summary)
}
#-------------------------------------------------------------------------------
#create a function that finds the max value for a group
group_max <- function(dataframe_column, variables){

  true_max <- 0
  true_min <- 0

  for (var in dataframe_column){
    cur_max <- variables[var, 'Max']
    cur_min <- variables[var, 'Min']

    if(cur_max > true_max){true_max <- cur_max}
    if(cur_min < true_min){true_min <- cur_min}
  }

  for (var in dataframe_column){
    variables[var, 'Max'] <- true_max
    variables[var, 'Min'] <- true_min
  }
  return(variables)
}
#-------------------------------------------------------------------------------
#Creates a function that combines singular variables into a combined variable
if(sort_by == 'region'){
  make_combined_var <- function(var1,var2,add_or_subtract){
    if(add_or_subtract == "add"){
      output <- dplyr::left_join(var1, var2, by = c("model", "experiment"))
      output <- dplyr::mutate(output, value = value.x + value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, experiment, value, sd))
      return(output)
    }
    if (add_or_subtract == "so2_and_so4_add"){
      output <- dplyr::left_join(var1, var2, by = c("model", "region"))
      output <- dplyr::mutate(output, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, region, value, sd))
      return(output)
    }
    if(add_or_subtract == "subtract"){
      output <- dplyr::left_join(var1, var2, by = c("model", "experiment"))
      output <- dplyr::mutate(output, value = value.x - value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, experiment, value, sd))
      return(output)
    }

    if(add_or_subtract == "divide"){
      output <- dplyr::left_join(var1, var2, by = c("model", "experiment"))
      output <- dplyr::mutate(output, value = value.x / value.y) %>%
        dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
        dplyr::select(c(model, experiment, value, sd))
      return(output)
    }
  }


}

if(sort_by == 'experiment'){
  make_combined_var <- function(var1,var2, add_or_subtract){
    if (add_or_subtract == "add"){
      output <- dplyr::left_join(var1, var2, by = c("model", "region"))
      output <- dplyr::mutate(output, value = value.x + value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, region, value, sd))
      return(output)
    }
    if (add_or_subtract == "so2_and_so4_add"){
      output <- dplyr::left_join(var1, var2, by = c("model", "region"))
      output <- dplyr::mutate(output, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, region, value, sd))
      return(output)
    }
    if (add_or_subtract == "subtract"){
      output <- dplyr::left_join(var1, var2, by = c("model", "region"))
      output <- dplyr::mutate(output, value = value.x - value.y) %>%
        dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
        dplyr::select(c(model, region, value, sd))
      return(output)
    }
    if(add_or_subtract == "divide"){
      output <- dplyr::left_join(var1, var2, by = c("model", "region"))
      output <- dplyr::mutate(output, value = value.x / value.y) %>%
        dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
        dplyr::select(c(model, experiment, value, sd))
      return(output)
    }

  }
}
#-------------------------------------------------------------------------------
#Creates a function that accumulates (using data_accumulation) data and appropriately
#renames specific columns (depending on what you are sorting by)
data_accum_and_rename_region <- function(direc, reg, experiment){
  experiment_subbed <- gsub('_','-',experiment)
  output <- data_accumulation(direc, reg, experiment_subbed)
  names(output)[names(output) == 'regional_data'] <- experiment
  exper_sd <- paste0(experiment, '_sd')
  names(output)[names(output) == 'regional_data_sd'] <- exper_sd
  return(output)
}

data_accum_and_rename_experiment <- function(direc, reg, experiment){
  reg_subbed <- gsub('_','-',reg)
  output <- data_accumulation(direc, reg_subbed, experiment)
  names(output)[names(output) == 'regional_data'] <- reg
  reg_sd <- paste0(reg, '_sd')
  names(output)[names(output) == 'regional_data_sd'] <- reg_sd
  return(output)
}

#create a function that sorts data for each region
if(sort_by == "region"){
  create_summary_long <- function(region, directory){

    summary_data_list <- list()

    for (i in 1:length(master_exp)){
      summary_data_list[[i]] <- assign(paste0(master_exp[i], '_summary'), data_accum_and_rename_region(directory,region,master_exp[i]))
    }

    # Bind data together
    #summary_data <- list(bc_no_season_summary, high_so4_summary, no_so4_summary, so2_at_height_summary, so2_no_season_summary) %>% reduce(left_join, by = c("variable", "model"))
    summary_data <- summary_data_list %>% reduce(left_join, by = c("variable", "model"))

    # Correct model names for CESM and CESM2
    summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"

    # Change to long format
    summary_long_exp <- summary_data %>%
      gather(experiment, value, -c(model, variable)) %>%
      filter(!grepl("_sd", experiment)) %>%
      select(variable, model, experiment, value) %>%
      drop_na()

    summary_long_sd <- summary_data %>%
      gather(experiment, sd, -c(model, variable)) %>%
      filter(grepl("_sd", experiment)) %>%
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
    return(summary_long)
  }
}

#create a function that sorts data for each experiment
if(sort_by == "experiment"){
  create_summary_long <- function(experiment, directory){

    summary_data_list <- list()

    for (i in 1:length(master_region)){
      summary_data_list[[i]] <- assign(paste0(master_region[i]), data_accum_and_rename_experiment(directory,master_region[i],experiment))
    }

    # Bind data together
    summary_data <- summary_data_list %>% reduce(left_join, by = c("variable", "model"))

    # Correct model names for CESM and CESM2
    summary_data$model[which(summary_data$model == "CESM")] <- "CESM1"

    # Change to long format
    summary_long_exp <- summary_data %>%
      gather(region, value, -c(model, variable)) %>%
      filter(!grepl("_sd", region)) %>%
      select(variable, model, region, value) %>%
      drop_na()

    summary_long_sd <- summary_data %>%
      gather(region, sd, -c(model, variable)) %>%
      filter(grepl("_sd", region)) %>%
      select(variable, model, region, sd) %>%
      drop_na()

    summary_long_sd$region <- gsub("_sd", "", summary_long_sd$region)

    summary_long <- dplyr::left_join(summary_long_exp, summary_long_sd)

    return(summary_long)
  }
}
#-------------------------------------------------------------------------------
#Creates a function that finds min and max values
find_max_min <- function(variable_data, variable, varname){

  #replace max and min if current variable data is more or less than the current max/min
  if(max(variable$value) > variable_data[varname, 'Max']){
    variable_data[varname, 'Max'] <- max(variable$value) + max(variable$sd)
  }
  if(min(variable$value) < variable_data[varname, 'Min']){
    variable_data[varname, 'Min'] <- min(variable$value) - max(variable$sd)
  }

  return(variable_data)
}
#-------------------------------------------------------------------------------
#Creates a function that filters species
#creates a function that filters species out of a database
filter_species <- function(database, species){
  species <- dplyr::filter(database, variable == species)
  return(species)
}
#-------------------------------------------------------------------------------
#Creates functions that process all data if fixing the y-axis
if (sort_by == "region"){

  if(fixed_data[1, 'fixed_by'] == 'fixed'){
    #first determine min and max values
    #finds min and max for each variable in all scenarios
    for (current_region in master_region){
      #substitutes all _ for - so that R can trace through the directory
      current_region <- gsub('_','-',current_region)
      #run the create summary long function for each region
      region_summary_long <- create_summary_long(current_region, emi_dir)

      list_vars <- list()
      #filter out each variable and determine max and min
      for (i in 1:length(master_vars)){
        list_vars[[i]] <- assign(master_vars[i], filter_species(region_summary_long,master_vars[i]))
      }

      for(i in 1:ncol(combined_vars)){
        name <- colnames(combined_vars)[i]
        var1 <- combined_vars[1,i]
        var2 <- combined_vars[2,i]
        operator <- combined_vars[3,i]
        #Adds these variables to total vars directly under where it left off
        list_vars[[i + length(master_vars)]] <- assign(name, make_combined_var(eval(parse(text = var1)),eval(parse(text = var2)),operator))
        intermediate_df <- list_vars[[i + length(master_vars)]]
        list_vars[[i + length(master_vars)]] <- intermediate_df %>%  dplyr::filter(is.infinite(value) == FALSE) %>% na.omit()
      }
      #creates an arbitrary counter variable
      j <- 1
      #runs the find_max_min function for each variable included in the csv file
      for (var in list_vars){
        variables <- find_max_min(variables, var, list_of_variable_strings[j])
        j <- j + 1
      }
    }
  }
  #---------------------------------------------------------------------------
  #Now can get data for the region being studied

  #creates a summary_long with the data being focused on
  summary_long <- create_summary_long(region, emi_dir)

  # source("C:/Users/such559/Documents/Emissions-MIP_Data/code/cobine_so2_so4.R")
  #
  # summed_so2_so4_adjusted <- summed_so2_so4[[paste0('so4_so2_sum_',region)]] %>%
  #   dplyr::rename(value = regional_data) %>%
  #   dplyr::rename(sd = regional_data_sd)
  # summed_so2_so4_adjusted <- summed_so2_so4_adjusted[,c('variable','model','experiment','value','sd')]
  #
  # summary_long <- summary_long %>%
  #   rbind(summed_so2_so4_adjusted)

  if(fixed_data[1, 'fixed_by'] == 'group'){
    #filter out each variable and determine max and min

    #Creates a list of singular variables by filtering from summary_long
    list_vars <- list()
    for (i in 1:length(master_vars)){
      list_vars[[i]] <- assign(master_vars[i], filter_species(summary_long,master_vars[i]))
    }
    #Adds the combined variables to the list
    for(i in 1:ncol(combined_vars)){
      name <- colnames(combined_vars)[i]
      var1 <- combined_vars[1,i]
      var2 <- combined_vars[2,i]
      operator <- combined_vars[3,i]
      #Adds these variables to total vars directly under where it left off
      list_vars[[i + length(master_vars)]] <- assign(name, make_combined_var(eval(parse(text = var1)),eval(parse(text = var2)),operator))
      #creates an intermediate dataframe for editing meant to be overwritten each loop
      intermediate_df <- list_vars[[i + length(master_vars)]]
      #Gets rid of infinite, NA, and NaN values
      list_vars[[i + length(master_vars)]] <- dplyr::filter(intermediate_df,is.infinite(intermediate_df$value) == FALSE) %>% na.omit()
    }
    #establishes a j term used for indexing later on
    j <- 1

    #runs the find_max_min function for each variable included in the csv file
    for (var in list_vars){
      variables <- find_max_min(variables, var, list_of_variable_strings[j])
      j <- j + 1
    }

    #creates a vector with each fixed_data column name as a string, remove the first column
    fixed_colnames <- colnames(fixed_data)
    fixed_colnames <- fixed_colnames[-1]

    #only runs if there is at least one group in the csv file
    if (ncol(fixed_data > 1)){
      #runs through each existing group and groups the min/maxes for each var together
      for (col in fixed_colnames){
        current_vector <- fixed_data[,col]
        variables <- group_max(current_vector, variables)
      }
    }
  }
}


if (sort_by == "experiment"){

  if(fixed_data[1, 'fixed_by'] == 'fixed'){
    #-----------------------------------------------------------------------------
    #find min and max for each variable in all experiments
    all_experiments <- master_exp
    variables <- read.csv(file = paste0(emi_dir, '/input/variables.csv'), fileEncoding="UTF-8-BOM")
    rownames(variables) <- variables$Variable
    variables <- subset(variables, select = -c(Variable))

    #creates a list of all the variables as strings
    list_of_variable_strings <- rownames(variables)

    for (current_experiment in all_experiments){
      #substitutes all _ for - so that R can trace through the directory
      current_experiment <- gsub('_','-',current_experiment)
      #run the create summary long function for each region
      experiment_summary_long <- create_summary_long(current_experiment, emi_dir)
      #Create an empty variable dataframe list
      list_vars <- list()
      #Creates a list of dataframes filtered from summary_long
      for (i in 1:length(master_vars)){
        list_vars[[i]] <- assign(master_vars[i], filter_species(experiment_summary_long,master_vars[i]))
      }
      #Adds combined variables to the list
      for(i in 1:ncol(combined_vars)){
        name <- colnames(combined_vars)[i]
        var1 <- combined_vars[1,i]
        var2 <- combined_vars[2,i]
        operator <- combined_vars[3,i]
        #Adds these variables to total vars directly under where it left off
        list_vars[[i + length(master_vars)]] <- assign(name, make_combined_var(eval(parse(text = var1)),eval(parse(text = var2)),operator))
        #creates an intermediate dataframe for editing meant to be overwritten each loop
        intermediate_df <- list_vars[[i + length(master_vars)]]
        #Gets rid of infinite, NA, and NaN values
        list_vars[[i + length(master_vars)]] <- dplyr::filter(intermediate_df,is.infinite(intermediate_df$value) == FALSE) %>% na.omit()
      }
      #establishes a j term used for indexing later on
      j <- 1

      for (var in list_vars){

        variables <- find_max_min(variables, var, list_of_variable_strings[j])

        j <- j + 1
      }
    }
  }
  #-----------------------------------------------------------------------------
  #Now can get data for experiment being studied
  #read in data for each region
  summary_long <- create_summary_long(exper, emi_dir)

  #source("C:/Users/such559/Documents/Emissions-MIP_Data/code/cobine_so2_so4.R")
#
#   summed_so2_so4_adjusted <- summed_so2_so4[[paste0('so4_so2_sum_',region)]] %>%
#     dplyr::rename(value = regional_data) %>%
#     dplyr::rename(sd = regional_data_sd)
#   summed_so2_so4_adjusted <- summed_so2_so4_adjusted[,c('variable','model','experiment','value','sd')]
#
#   summary_long <- summary_long %>%
#     rbind(summed_so2_so4_adjusted)

  if(fixed_data[1, 'fixed_by'] == 'group'){

    #Create an empty variable dataframe list
    list_vars <- list()
    #filter out each variable and determine max and min
    for (i in 1:length(master_vars)){
      list_vars[[i]] <- assign(master_vars[i], filter_species(summary_long,master_vars[i]))
    }
    #Adds combined variables to the list
    for(i in 1:ncol(combined_vars)){
      name <- colnames(combined_vars)[i]
      var1 <- combined_vars[1,i]
      var2 <- combined_vars[2,i]
      operator <- combined_vars[3,1]
      #Adds these variables to total vars directly under where it left off
      list_vars[[i + length(master_vars)]] <- assign(name, make_combined_var(eval(parse(text = var1)),eval(parse(text = var2)),operator))
      #creates an intermediate dataframe for editing meant to be overwritten each loop
      intermediate_df <- list_vars[[i + length(master_vars)]]
      #Gets rid of infinite, NA, and NaN values
      list_vars[[i + length(master_vars)]] <- dplyr::filter(intermediate_df,is.infinite(intermediate_df$value) == FALSE) %>% na.omit()
    }

    #establishes a j term used for indexing later on
    j <- 1

    #runs the find_max_min function for each variable included in the csv file
    for (var in list_vars){
      variables <- find_max_min(variables, var, list_of_variable_strings[j])
      j <- j + 1
    }

    #creates a vector with each fixed_data column name as a string, remove the first column
    fixed_colnames <- colnames(fixed_data)
    fixed_colnames <- fixed_colnames[-1]

    #only runs if there is at least one group in the csv file
    if (ncol(fixed_data > 1)){
      #runs through each existing group and groups the min/maxes for each var together
      for (col in fixed_colnames){
        current_vector <- fixed_data[,col]
        variables <- group_max(current_vector, variables)
      }
    }

  }
}
#-------------------------------------------------------------------------------

# Generate plots
title_font <- 9
axis_font <- 9
axis_title_font <- 9

if(fixed_data[1, 'fixed_by'] == 'group'){ total_vars <- list_vars}

if(fixed_data[1, 'fixed_by'] != 'group'){
  #create a list of variable dataframes
  total_vars <- list()
  #filters each species from summary_long
  for (i in 1:length(master_vars)){
    total_vars[[i]] <- assign(master_vars[i], filter_species(summary_long,master_vars[i]))
  }

  #creates the combined variables
  for(i in 1:ncol(combined_vars)){
    name <- colnames(combined_vars)[i]
    var1 <- combined_vars[1,i]
    var2 <- combined_vars[2,i]
    operator <- combined_vars[3,i]
    #Adds these variables to total vars directly under where it left off
    total_vars[[i + length(master_vars)]] <- assign(name, make_combined_var(eval(parse(text = var1)),eval(parse(text = var2)),operator))
    #creates an intermediate dataframe for editing meant to be overwritten each loop
    intermediate_df <- total_vars[[i + length(master_vars)]]
    total_vars[[i + length(master_vars)]] <- dplyr::filter(intermediate_df,is.infinite(intermediate_df$value) == FALSE) %>% na.omit()
  }
}

#filter out excluded data again (previous filter step did not include combined variables)
all_vars <- rownames(variables)
if(nrow(excluded_models) !=0){
  for(Model_name in 1:nrow(excluded_models)){
    #locates index number where excluded variable is located in total_vars
    model_index <- which(all_vars == excluded_models[Model_name,3])
    #filters out the experiment and model associated with the filtered variable
    total_vars[[model_index]] <- filter(total_vars[[model_index]], experiment != excluded_models$Scenario[Model_name] | model != excluded_models$Model[Model_name]  )
  }
}
#convert so2_timescale and so4_lifetime units from seconds to days
so2_timescale_index <- which(all_vars == 'so2_timescale')
so4_lifetime_index <- which(all_vars == 'so4_lifetime')
so2_timescale_to_fix <- total_vars[[so2_timescale_index]]
so4_lifetime_to_fix <- total_vars[[so4_lifetime_index]]
total_vars[[so2_timescale_index]] <- mutate(so2_timescale_to_fix, value = value / 86400, sd = sd / 86400)
total_vars[[so4_lifetime_index]] <- mutate(so4_lifetime_to_fix, value = value / 86400, sd = sd / 86400)
#Convert the max min data from seconds to days as well
variables['so2_timescale','Max'] <- variables['so2_timescale','Max'] / 86400
variables['so2_timescale','Min'] <- variables['so2_timescale','Min'] / 86400
variables['so4_lifetime','Max'] <- variables['so4_lifetime','Max'] / 86400
variables['so4_lifetime','Min'] <- variables['so4_lifetime','Min'] / 86400

test <- as.character(variables[2,3])

#Creates a function that creates plots for the data based on each species
if (sort_by == "region"){
  #plots only take max and min into account if fixing or grouping
  if (fixed_data[1, 'fixed_by'] == 'fixed' || fixed_data[1, 'fixed_by'] == 'group'){
    plot_species <- function(variable, x, y, title1,title2, units, region_or_exper, model_colors, model_symbols, ymin, ymax){
      species <- variable
      species_plot <- ggplot(species, aes(x = experiment, y = value, color = model, shape = model))+
        theme_bw()+
        labs(title=paste0(title1,"\n", title2,' - ', region_or_exper), y=units) +
        theme(plot.title = element_text(hjust = 0.5, size = title_font),
              axis.text = element_text(size = axis_font),
              axis.title = element_text(size = axis_title_font),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank()) +
        scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd))))+
        scale_colour_manual(values = model_colors) +
        scale_shape_manual(values = model_symbols) +
        geom_point( position=position_dodge(width=0.5), size = 1.3) +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), size=0.1, width = 0, position=position_dodge(0.5), show.legend = F) +
        ylim(ymax, ymin)

      return(species_plot)
    }
  }else{
    plot_species <- function(variable, x, y, title1,title2, units, region_or_exper, model_colors, model_symbols, ymin, ymax){
      species <- variable
      species_plot <- ggplot(species, aes(x = experiment, y = value, color = model, shape = model))+
        theme_bw()+
        labs(title=paste0(title1, "\n",title2,' - ', region_or_exper), y=units) +
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
  }
  #creates plots based on each species using the plot_species function
  k <- 1
  for(var in total_vars){
    #number of the row in the variables column to access
    #rownum <- which(rownames(variables) == master_vars[k])
    assign(paste0(all_vars[k],'_plot'), plot_species(var,region,value,variables[k,'plot_names'],variables[k,'plot_names_line_2'],eval(parse( text = paste(variables[k,'description_diff']))),region,model_colors,model_symbols,variables[k,'Max'],variables[k,'Min']))

    k <- k + 1
  }


}

if (sort_by == "experiment"){
  #plots only take min or max into account if grouping
  if (fixed_data[1, 'fixed_by'] == 'fixed' || fixed_data[1, 'fixed_by'] == 'group'){
    plot_species <- function(variable, x, y, title1,title2, units, region_or_exper, model_colors, model_symbols, ymin, ymax){
      species <- variable
      species_plot <- ggplot(species, aes(x = region, y = value, color = model, shape = model))+
        theme_bw()+
        labs(title=paste0(title1,"\n",title2,' - ', region_or_exper), y=units) +
        theme(plot.title = element_text(hjust = 0.5, size = title_font),
              axis.text = element_text(size = axis_font),
              axis.title = element_text(size = axis_title_font),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank()) +
        scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(-max(abs(species$value))-max(abs(species$sd)), max(abs(species$value))+max(abs(species$sd))))+
        scale_colour_manual(values = model_colors) +
        scale_shape_manual(values = model_symbols) +
        geom_point( position=position_dodge(width=0.4), size = 1.5) +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.2, position=position_dodge(0.4), show.legend = F)+
        ylim(ymax, ymin)

      return(species_plot)
    }
  } else{
    plot_species <- function(variable, x, y, title,title2, units, region_or_exper, model_colors, model_symbols, ymin, ymax){
      species <- variable
      species_plot <- ggplot(species, aes(x = region, y = value, color = model, shape = model))+
        theme_bw()+
        labs(title=paste0(title,"\n", title2,' - ', region_or_exper), y=units) +
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
  }

  #creates plots based on each species using the plot_species function
  k <- 1
  for(var in total_vars){
    #number of the row in the variables column to access
    #rownum <- which(rownames(variables) == master_vars[k])
    all_vars <- rownames(variables)
    assign(paste0(all_vars[k],'_plot'), plot_species(var,region,value,variables[k,'plot_names'],variables[k,'plot_names_line_2'],eval(parse(text = paste(variables[k,'description_diff']))),exper,model_colors,model_symbols,variables[k,'Max'],variables[k,'Min'] ))

    k <- k + 1
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
    heights = unit.c(unit(1, "npc") - 1.5 * lheight, lheight), # the "1.5" adds room for title
    top = textGrob("Summary - absolute difference", gp = gpar(fontsize = 12)))
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
                                          so2_timescale_plot)

# Print plots
if (sort_by == 'region'){
  setwd(paste0('../../../../output/', region, '/summary'))

  pdf(paste0(region, '_summary_plots_diff.pdf'), height = 11, width = 10)
}

if (sort_by == 'experiment'){
  setwd(paste0('../../../../output/', exper, '/summary'))

  pdf(paste0(exper, '_summary_plots_diff.pdf'), height = 11, width = 12)
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


