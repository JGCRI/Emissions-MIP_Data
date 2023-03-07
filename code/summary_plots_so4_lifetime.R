# ------------------------------------------------------------------------------
# Program Name: summary_plots_so4_lifetime.R
# Authors: Hamza Ahsan
# Date Last Modified: March 14, 2022
# Program Purpose: Produces summary plots for global sulfate lifetime
# Input Files: ~Emissions-MIP/input/
# Output Files: ~Emissions-MIP/output/
# TODO:
# ------------------------------------------------------------------------------
# 1. Setup directories and parameters

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# Specify location of Emissions-MIP directory
emi_dir <- paste0('C:/Users/such559/Documents/Emissions-MIP_Data')
setwd(paste0(emi_dir))

output_dir <- 'C:/Users/such559/Documents/Emissions-MIP_Data/output'

# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#920000", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#490092", "#117733")

model_colors <- c(CESM1 = cbPalette[1], E3SM = cbPalette[2], GISS = cbPalette[3],
                  CESM2 = cbPalette[4], MIROC = cbPalette[5], NorESM2 = cbPalette[6],
                  GFDL = cbPalette[7], OsloCTM3 = cbPalette[8], UKESM = cbPalette[9],
                  GEOS = cbPalette[10], CAM5ATRAS = cbPalette[11])

model_symbols <- c(CESM1 = 15, E3SM = 15, GISS = 17, CESM2 = 19, MIROC = 15,
                   NorESM2 = 17, GFDL = 19, OsloCTM3 = 19, UKESM = 15, GEOS = 17,
                   CAM5ATRAS = 17)

#-------------------------------------------------------------------------------
# 2. Generate global reference case plot for sulfate lifetime

region <- "global"

# Setup directory for reference case data
setwd(paste0(emi_dir, '/input/', region, '/reference/diff'))

# Read in csv files and bind into single data frame
target_filename <- list.files(getwd(), "*.csv")
ref_data <- rbind(map(target_filename, read.csv))
ref_data <- lapply(ref_data, function(x) {x["unit"] <- NULL; x})
ref_data <- bind_rows(ref_data)

# Extract model from file names (fifth segment) and bind to experiment data frame
models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[4])
rep_models <- rep(models, each = 5) # five years
ref_data$model <- rep_models

# Correct model names
ref_data$model[which(ref_data$model == "CMIP6_AerChemMIP")] <- "GISS"
ref_data$model[which(ref_data$model == "CMIP6_CMIP_CESM")] <- "CESM1"
ref_data$model[which(ref_data$model == "CMIP6_CMIP_E3SM")] <- "E3SM"
ref_data$model[which(ref_data$model == "CMIP6_CMIP_CESM2")] <- "CESM2"
ref_data$model[which(ref_data$model == "CAM5")] <- "CAM5ATRAS"

# Convert SO2 volume mixing ratio to mass mixing ratio by multiplying by molar
# mass of SO2 and dividing by molar mass of air, invert sign of forcing variables
# to be consistent with convention (i.e. positive value denotes a heating effect),
# then take the average over all years for each variable and calculate std dev
ref_data_summary <- ref_data %>%
  dplyr::group_by(variable, model) %>%
  within(value <- ifelse(variable == "so2", 64.066 / 28.96, 1) * value) %>%
  within(value <- ifelse(variable == "dms", 62.13 / 28.96, 1) * value) %>%
  within(value <- ifelse(variable %in% c("rlut", "rsut", "rlutcs", "rsutcs"), -1, 1) * value) %>%
  within(value <- ifelse(variable %in% c("wetbc", "wetso2", "wetso4") & model == "CESM2", -1, 1) * value) %>%
  dplyr::summarise(ref_data = mean(value), ref_data_sd = sd(value))

# Define total sulfate deposition (dryso4 + wetso4) for sulfate lifetime
dryso4 <- dplyr::filter(ref_data_summary, variable == "dryso4")
wetso4 <- dplyr::filter(ref_data_summary, variable == "wetso4")
tot_so4 <- dplyr::left_join(dryso4, wetso4, by = "model")
tot_so4 <- dplyr::mutate(tot_so4, ref_data = ref_data.x + ref_data.y) %>%
  dplyr::mutate(ref_data_sd = sqrt(ref_data_sd.x^2 + ref_data_sd.y^2)) %>%
  dplyr::select(c(model, ref_data, ref_data_sd))

# Define sulfate lifetime (loadso4/(dryso4+wetso4)) and convert from seconds to days
loadso4 <- dplyr::filter(ref_data_summary, variable == "loadso4")
so4_lifetime <- dplyr::left_join(loadso4, tot_so4, by = "model")
so4_lifetime <- dplyr::mutate(so4_lifetime, ref_data = (ref_data.x/ref_data.y)/86400) %>%
  dplyr::mutate(ref_data_sd = ref_data*sqrt((ref_data_sd.x/ref_data.x)^2 + (ref_data_sd.y/ref_data.y)^2)) %>%
  dplyr::filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  dplyr::filter(ref_data >= 0) %>%
  dplyr::select(c(model, ref_data, ref_data_sd)) %>%
  dplyr::mutate(variable = "so4_lifetime")

# Reference plots
so4_lifetime_ref_plot <- ggplot(so4_lifetime, aes(x = model, y = ref_data, color = model, shape = model))+
  theme_bw()+
  labs(title = "Global Reference", y="sulfate lifetime (days)") +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1)))+
  scale_colour_manual(values = model_colors) +
  scale_shape_manual(values = model_symbols) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = ref_data - ref_data_sd, ymax = ref_data + ref_data_sd), size = 0.3, width = 0, position=position_dodge(0.5), show.legend = F) +
  annotate(geom = 'text', size = 6, label = '(a)', x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5)

#-------------------------------------------------------------------------------
# 3. Generate global difference plot for sulfate lifetime

pert_data_gather <- function(emi_dir, region, exper){

  setwd(paste0(emi_dir,'/input/', region,'/', exper, '/reference/diff'))

  # Read in csv files and bind into single data frame
  target_filename <- list.files(getwd(), "*.csv")
  pert_data <- rbind(map(target_filename, read.csv))
  pert_data <- lapply(pert_data, function(x) {x["unit"] <- NULL; x})
  pert_data <- bind_rows(pert_data)

  # Extract model from file names (fifth segment) and bind to experiment data frame
  models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[4])
  rep_models <- rep(models, each = 5) # five years
  pert_data$model <- rep_models

  # Correct model names
  pert_data$model[which(pert_data$model == "CMIP6_AerChemMIP")] <- "GISS"
  pert_data$model[which(pert_data$model == "CMIP6_CMIP_CESM")] <- "CESM1"
  pert_data$model[which(pert_data$model == "CMIP6_CMIP_E3SM")] <- "E3SM"
  pert_data$model[which(pert_data$model == "CMIP6_CMIP_CESM2")] <- "CESM2"
  pert_data$model[which(pert_data$model == "CAM5")] <- "CAM5ATRAS"

  # Calculate time scale variables for perturbations
  pert_data_wide <- pert_data %>%
    spread(variable, value) %>%
    within(wetso4 <- ifelse(model == "CESM2", -1, 1) * wetso4) %>%
    mutate(so4_lifetime = (loadso4/(dryso4+wetso4))/86400,
           so2_lifetime = (loadso2/emiso2)/86400)

  # Convert back to long format
  pert_data_long <- gather(pert_data_wide, variable, value, -c(year, model))

  return(pert_data_long)
}

shp_30p_red_data <- pert_data_gather(emi_dir,region,"shp_30p_red")
shp_60p_red_data <- pert_data_gather(emi_dir,region,"shp-60p-red")
shp_60p_red_1950_data <- pert_data_gather(emi_dir,region,"shp-60p-red-1950")
atl_shift_1950_data <- pert_data_gather(emi_dir,region,"atl-shift-1950")
ind_shift_1950_data <- pert_data_gather(emi_dir,region,"ind-shift-1950")

# Calculate time scale variables for reference
ref_data_wide <- ref_data %>%
  spread(variable, value) %>%
  within(wetso4 <- ifelse(model == "CESM2", -1, 1) * wetso4) %>%
  mutate(so4_lifetime = (loadso4/(dryso4+wetso4))/86400,
         so2_lifetime = (loadso2/emiso2)/86400)

# Convert back to long format
ref_data_long <- gather(ref_data_wide, variable, value, -c(year, model))

# Join reference and perturbation and calculate difference and std dev
summary_data <- left_join(ref_data_long, shp_30p_red_data, by = c("variable", "year", "model")) %>%
  left_join(shp_60p_red_data, by = c("variable", "year", "model")) %>%
  left_join(shp_60p_red_1950_data, by = c("variable", "year", "model")) %>%
  left_join(atl_shift_1950_data, by = c("variable", "year", "model")) %>%
  left_join(ind_shift_1950_data, by = c("variable", "year", "model")) %>%
  rename(ref = value.x, shp_30p_red = value.y, shp_60p_red = value.x.x, shp_60p_red_1950 = value.y.y, atl_shift_1950 = value.x.x.x, ind_shift_1950 = value.y.y.y) %>%
  drop_na() %>%
  mutate(shp_30p_red_diff = shp_30p_red - ref) %>%
  mutate(shp_60p_red_diff = shp_60p_red - ref) %>%
  mutate(shp_60p_red_1950_diff = shp_60p_red_1950 - ref) %>%
  mutate(atl_shift_1950_diff = atl_shift_1950 - ref) %>%
  mutate(ind_shift_1950_diff = ind_shift_1950 - ref) %>%
  select(-ref) %>%
  group_by(variable, model) %>%
  dplyr::summarise(shp_30p_red = mean(shp_30p_red_diff), shp_30p_red_sd = sd(shp_30p_red_diff),
                   shp_60p_red = mean(shp_60p_red_diff), shp_60p_red_sd = sd(shp_60p_red_diff),
                   shp_60p_red_1950 = mean(shp_60p_red_1950_diff), shp_60p_red_1950_sd = sd(shp_60p_red_1950_diff),
                   ind_shift_1950 = mean(ind_shift_1950_diff), ind_shift_1950_sd = sd(ind_shift_1950_diff),
                   atl_shift_1950 = mean(atl_shift_1950_diff), atl_shift_1950_sd = sd(atl_shift_1950_diff))

# Change to long format
data_long <- summary_data %>%
  gather(experiment, value, -c(model, variable, shp_30p_red_sd, shp_60p_red_sd, shp_60p_red_1950_sd, ind_shift_1950_sd, atl_shift_1950_sd)) %>%
  select(variable, model, experiment, value) %>%
  drop_na()

sd_long <- summary_data %>%
  gather(experiment, sd, -c(model, variable, shp_30p_red, ind_shift_1950, shp_60p_red, shp_60p_red_1950, atl_shift_1950)) %>%
  select(variable, model, experiment, sd) %>%
  drop_na()

sd_long$experiment <- gsub("_sd", "", sd_long$experiment)

summary_long <- dplyr::left_join(data_long, sd_long, by = c("variable", "model", "experiment"))

# Set order of experiment
summary_long$experiment <- factor(summary_long$experiment,levels = c("shp_30p_red", "ind_shift_1950", "shp_60p_red_1950", "shp_60p_red", "atl_shift_1950"))

# Difference plot
so4_lifetime <- dplyr::filter(summary_long, variable == "so4_lifetime")

so4_lifetime_diff_plot <- ggplot(so4_lifetime, aes(x = experiment, y = value, color = model, shape = model))+
  theme_bw()+
  labs(title = "Global Perturbation - Reference" , y=expression(Delta*~sulfate~lifetime~(days))) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-max(abs(so4_lifetime$value))-max(abs(so4_lifetime$sd)), max(abs(so4_lifetime$value))+max(abs(so4_lifetime$sd))))+
  scale_colour_manual(values = model_colors) +
  scale_shape_manual(values = model_symbols) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), size = 0.3, width = 0, position=position_dodge(0.5), show.legend = F) +
  annotate(geom = 'text', size = 6, label = '(b)', x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5)


