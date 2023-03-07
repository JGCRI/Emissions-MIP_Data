# ------------------------------------------------------------------------------
# Program Name: combine_so2_so4.R
# Authors: Harrison Suchyta
# Date Last Modified: September 15, 2022
# Program Purpose: this script combines the high SO2 and SO4 at height values in order to compare
# that value with the high-SO4-at-height experiment
# Input Files: ~Emissions-MIP/input/
# Output Files: ~Emissions-MIP/output/
# TODO: ADD GISS functionality
# ------------------------------------------------------------------------------

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)




#set emissions directory
emi_dir <- paste0('C:/Users/such559/Documents/Emissions-MIP_Data')
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
#create a function to accumulate data from csv files
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
#Make a function to combine experiment outputs
combine_expers <- function(region_name,curr_names,summary_list){

    #select both dataframes associated with the current region
    curr_reg <- curr_names[grepl(region_name,curr_names) ==TRUE]

    curr_df1 <- summary_list[[curr_reg[1]]]
    curr_df2 <- summary_list[[curr_reg[2]]]

    #add the values and standard deviations
    sum_df <- curr_df1 %>%
        left_join(curr_df2,by = c('variable','model')) %>%
        mutate(regional_data = regional_data.x + regional_data.y,
               regional_data_sd = sqrt(regional_data_sd.x^2 + regional_data_sd.y^2),
               experiment = 'so2_so4_sum') %>%
        subset(select = -c(regional_data.x,regional_data.y,regional_data_sd.x,regional_data_sd.y))

    return(sum_df)

}
#-------------------------------------------------------------------------------
#read in high-So2 and SO4 at height values

experiments <- c('high-so4','so2-at-height')

summary_data_list <- list()

#loop through each region and extract values. Filter to just E3SM
for (exper in experiments){
    for (reg_name in master_region){
        summary_data_list[[paste0(exper,'_summary_',reg_name)]] <- assign(paste0(exper,'_summary_',reg_name),data_accumulation(emi_dir,reg_name,exper)) %>%
            filter(model %in% c('E3SM','GISS'))
    }
}

#create a list of the names in summary_data_list
data_names <- names(summary_data_list)

summed_so2_so4 <- list()
#loop through each region and add the results together
for (reg_name in master_region){

    summed_so2_so4[[paste0('so4_so2_sum_',reg_name)]] <- assign(paste0('so4_so2_sum_',reg_name),combine_expers(reg_name,data_names,summary_data_list))

}
