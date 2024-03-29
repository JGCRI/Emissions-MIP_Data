# ------------------------------------------------------------------------------
# Program Name: timeseries_diff.R
# Authors: Hamza Ahsan
# Date Last Modified: September 22, 2023
# Program Purpose: Produces time series line plots of the difference between
# the perturbations and the reference case
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

#set the working directory to the code directory
setwd('C:/Users/ahsa361/Documents/Emissions-MIP_Data')

# Specify location of Emissions-MIP directory
emi_dir <- paste0('C:/Users/ahsa361/Documents/Emissions-MIP_Data')

#determines whether the script sorts by region or experiment if nothing is put into command line
sort_by <- 'experiment'

#determines which region or experiment is sorted by
region <- 'global'
pert <- 'so2-no-season'

# Specify what you are sorting by and either the region (i.e., global, land, sea, arctic, NH-land, NH-sea, SH-land, SH-sea) or experiment (i.e., bc-no-season, high-so4, no-so4, reference, so2-at-height, so2-no-season)
#The command line would look like: rscript <rscript>.r <"experiment" or "region"> <specific experiment or region you are sorting by>
sorting <- commandArgs(trailingOnly = TRUE) #pulling region from command line
if (!is_empty(sorting)){sort_by <- sorting[1]}
if (!is_empty(sorting) & sort_by == "region"){region <- sorting[2]}
if (!is_empty(sorting) & sort_by == "experiment"){experiment <- sorting[2]}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#920000", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#490092", "#117733")

model_colors <- c('CESM' = cbPalette[1], 'E3SM' = cbPalette[2], 'GISS modelE' = cbPalette[3],
                  'CESM2' = cbPalette[4], 'MIROC-SPRINTARS' = cbPalette[5], 'NorESM2' = cbPalette[6],
                  'GFDL-ESM4' = cbPalette[7], 'OsloCTM3' = cbPalette[8], 'UKESM1' = cbPalette[9],
                  'GEOS' = cbPalette[10], 'CAM-ATRAS' = cbPalette[11])


# ------------------------------------------------------------------------------
# Iterate over the different perturbation/regional experiments
if (sort_by == "region"){scenarios <- c('bc-no-season', 'high-so4', 'no-so4', 'so2-at-height', 'so2-no-season')}
if (sort_by == "experiment"){scenarios <- c("arctic", "global", "land", "NH-atlantic", "NH-indian", "NH-land", "NH-pacific", "NH-sea", "sea", "SH-land", "SH-sea")}
#-------------------------------------------------------------------------------
#reads in csv file specifying which models to exclude from the data
excluded_models <- read.csv(file = paste0(emi_dir, '/input', '/excluded_data.csv'), fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)
excluded_models %>% drop_na() #gets rid of any empty spaces
#-------------------------------------------------------------------------------

for(scenario in scenarios){
  # Specify location of difference data
  if (sort_by == "region"){
    setwd(paste0(emi_dir, '/input/', region, '/', scenario, '/diff'))
  }
  if (sort_by == "experiment"){
    setwd(paste0(emi_dir, '/input/', scenario, '/', pert, '/diff'))
  }
  
  # Read in csv files and bind into single data frame
  target_filename <- list.files(getwd(), "*.csv")
  experiment <- rbind(map(target_filename, read.csv))
  experiment <- lapply(experiment, function(x) {x["unit"] <- NULL; x})
  experiment <- bind_rows(experiment)
  
  # Extract model from file names (fifth segment) and bind to experiment data frame
  models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[5])
  rep_models <- rep(models, each = 5) # five years
  experiment$model <- rep_models
  
  # Correct model names
  experiment$model[which(experiment$model == "GISS")] <- "GISS modelE"
  experiment$model[which(experiment$model == "CAM5")] <- "CAM-ATRAS"
  experiment$model[which(experiment$model == "UKESM")] <- "UKESM1"
  experiment$model[which(experiment$model == "MIROC")] <- "MIROC-SPRINTARS"
  experiment$model[which(experiment$model == "GFDL")] <- "GFDL-ESM4"
  
  # Convert DMS from volume mixing ratio to mass mixing ratio
  experiment <- within(experiment, value <- ifelse(variable == "dms", 62.13 / 28.96, 1) * value)  %>%
    
    # Invert sign of forcing variables to be consistent with convention (i.e. positive
    # value denotes a heating effect)
    within(value <- ifelse(variable %in% c("rlut", "rsut", "rlutcs", "rsutcs"), -1, 1) * value)  %>%
    
    # Invert sign of CESM2 wet deposition variables (i.e., CESM2 wetbc, wetso2, wetso4)
    within(value <- ifelse(variable %in% c("wetbc", "wetso2", "wetso4") & model == "CESM2", -1, 1) * value) %>%
    
    # Convert from NH4HSO4 to SO4 mass
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "E3SM", 96/115, 1) * value) %>%
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "CESM", 96/115, 1) * value) %>%
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "MIROC-SPRINTARS", 96/115, 1) * value) %>%
    
    # Convert from H2SO4 to SO4 mass
    within(value <- ifelse(variable %in% c("dryso4", "loadso4", "mmrso4", "wetso4") & model == "NorESM2", 96/98, 1) * value)
  
  # Rearrange data frame by years descending
  experiment <- dplyr::arrange(experiment, year)
  
  # Convert volume mixing ratio to mass mixing ratio by multiplying by molar mass
  # of SO2 and dividing by molar mass of air
  so2_experiment <- dplyr::filter(experiment, variable == 'so2') %>%
    dplyr::mutate(new_value = value * 64.066 / 28.96)
  
  # Change units from mol/mol to kg/kg
  so2_experiment$unit <- 'kg kg-1'
  
  # Define remaining experiments
  emibc_experiment    <- dplyr::filter(experiment, variable == 'emibc')
  emiso2_experiment   <- dplyr::filter(experiment, variable == 'emiso2')
  mmrbc_experiment    <- dplyr::filter(experiment, variable == 'mmrbc')
  mmrso4_experiment   <- dplyr::filter(experiment, variable == 'mmrso4')
  rlut_experiment     <- dplyr::filter(experiment, variable == 'rlut')
  rlutcs_experiment   <- dplyr::filter(experiment, variable == 'rlutcs')
  rsdt_experiment     <- dplyr::filter(experiment, variable == 'rsdt')
  rsut_experiment     <- dplyr::filter(experiment, variable == 'rsut')
  rsutcs_experiment   <- dplyr::filter(experiment, variable == 'rsutcs')
  drybc_experiment    <- dplyr::filter(experiment, variable == 'drybc')
  wetbc_experiment    <- dplyr::filter(experiment, variable == 'wetbc')
  dryso2_experiment   <- dplyr::filter(experiment, variable == 'dryso2')
  wetso2_experiment   <- dplyr::filter(experiment, variable == 'wetso2')
  dryso4_experiment   <- dplyr::filter(experiment, variable == 'dryso4')
  wetso4_experiment   <- dplyr::filter(experiment, variable == 'wetso4')
  od550aer_experiment <- dplyr::filter(experiment, variable == 'od550aer')
  clt_experiment      <- dplyr::filter(experiment, variable == 'clt')
  cltc_experiment     <- dplyr::filter(experiment, variable == 'cltc')
  
  # Define normal and clear-sky net radiative flux and  (sum of longwave and shortwave radiation)
  net_rad <- dplyr::left_join(rlut_experiment, rsut_experiment, by = c("year", "model"))
  net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
    dplyr::select(c(year, model, value))
  
  net_rad_cs <- dplyr::left_join(rlutcs_experiment, rsutcs_experiment, by = c("year", "model"))
  net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
    dplyr::select(c(year, model, value))
  
  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("year", "model"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::select(c(year, model, value))
  
  # Define total BC deposition rate (sum of dry BC and wet BC deposition)
  tot_bc <- dplyr::left_join(drybc_experiment, wetbc_experiment, by = c("year", "model"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::select(c(year, model, value))
  
  #Define total S deposition rate (sum of dry SO2/SO4 and wet SO2/SO4 deposition)
  dry_s <- dplyr::left_join(dryso2_experiment, dryso4_experiment, by = c("year", "model"))
  dry_s <- dplyr::mutate(dry_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::select(c(year, model, value))
  
  wet_s <- dplyr::left_join(wetso2_experiment, wetso4_experiment, by = c("year", "model"))
  wet_s <- dplyr::mutate(wet_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::select(c(year, model, value))
  
  tot_s <- dplyr::left_join(dry_s, wet_s, by = c("year", "model"))
  tot_s <- dplyr::mutate(tot_s, value = value.x + value.y) %>%
    dplyr::select(c(year, model, value))
  
  # Pre-define plot font sizes
  title_font <- 7
  axis_font <- 6
  axis_title_font <- 7
  
  plot_species <- function(experiment, reg_or_exper,title, units){
    species_plot <- ggplot(experiment, aes(x = year, y = value, color = model)) +
      labs(title=paste0(title, ' - ', reg_or_exper), y=units, x="Year") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = title_font),
            axis.text = element_text(size = axis_font),
            axis.title = element_text(size = axis_title_font)) +
      scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
      scale_colour_manual(values = model_colors) +
      geom_line()
    
    return(species_plot)
  }
  if (sort_by == "region"){
    #Generate plots
    emibc_plot <- plot_species(emibc_experiment, region, "surface flux \n of BC", expression(emibc~(kg~m^-2~s^-1)))
    emiso2_plot <- plot_species(emiso2_experiment, region, "surface flux \n of SO2",expression(emiso2~(kg~m^-2~s^-1)))
    mmrbc_plot <- plot_species(mmrbc_experiment, region, "surface concentration \n of BC",expression(mmrbc~(kg~kg-1)))
    mmrso4_plot <- plot_species(mmrso4_experiment, region, "surface concentration \n of SO4",expression(mmrso4~(kg~kg-1)))
    rlut_plot <- plot_species(rlut_experiment, region, "upwelling longwave flux \n at TOA",expression(rlut~(W~m-2)))
    rlutcs_plot <- plot_species(rlutcs_experiment, region, "upwelling clear-sky longwave \n flux at TOA",expression(rlutcs~(W~m-2)))
    rsut_plot <- plot_species(rsut_experiment, region, "upwelling shortwave flux \n at TOA",expression(rsut~(W~m-2)))
    rsutcs_plot <- plot_species(rsutcs_experiment, region, "upwelling clear-sky shortwave \n flux at TOA",expression(rsutcs~(W~m-2)))
    rsdt_plot <- plot_species(rsdt_experiment, region, "incident shortwave flux \n at TOA",expression(rsdt~(W~m-2)))
    net_rad_plot <- plot_species(net_rad, region, "net radiative flux \n at TOA",expression(rlut~+~rsut~(W~m^-2)))
    net_rad_cs_plot <- plot_species(net_rad_cs, region, "clear-sky net radiative \n flux at TOA",expression(rlutcs~+~rsutcs~(W~m^-2)))
    so2_plot <- plot_species(so2_experiment, region, "surface concentration \n of SO2",expression(so2~(kg~kg-1)))
    drybc_plot <- plot_species(drybc_experiment, region, "dry deposition rate \n of BC",expression(drybc~(kg~m^-2~s^-1)))
    wetbc_plot <- plot_species(wetbc_experiment, region, "wet deposition rate \n of BC",expression(wetbc~(kg~m^-2~s^-1)))
    tot_bc_plot <- plot_species(tot_bc, region, "total deposition rate \n of BC", expression(drybc~+~wetbc~(kg~m^-2~s^-1)))
    dryso2_plot <- plot_species(dryso2_experiment, region, "dry deposition rate \n of SO2", expression(dryso2~(kg~m^-2~s^-1)))
    wetso2_plot <- plot_species(wetso2_experiment, region, "wet deposition rate \n of SO2", expression(wetso2~(kg~m^-2~s^-1)))
    dryso4_plot <- plot_species(dryso4_experiment, region, "dry deposition rate \n of SO4", expression(dryso4~(kg~m^-2~s^-1)))
    wetso4_plot <- plot_species(wetso4_experiment, region, "wet deposition rate \n of SO4", expression(wetso4~(kg~m^-2~s^-1)))
    tot_s_plot <- plot_species(tot_s, region, "total deposition rate \n of S", expression(atop((dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3, (kg~m^-2~s^-1))))
    od550aer_plot <- plot_species(od550aer_experiment, region, "ambient aerosol optical \n thickness at 550nm", expression(od550aer))
    clt_plot <- plot_species(clt_experiment, region, "total cloud cover \n percentage", "expression clt (%)")
    cltc_plot <- plot_species(cltc_experiment, region, "convective cloud cover \n percentage", "expression cltc (%)")
    imp_cld_plot <- plot_species(imp_cld, region, "implied cloud response \n at TOA", expression(rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)))
  }
  
  if (sort_by == "experiment"){
    #Generate plots
    emibc_plot <- plot_species(emibc_experiment, pert, "surface flux \n of BC", expression(emibc~(kg~m^-2~s^-1)))
    emiso2_plot <- plot_species(emiso2_experiment, pert, "surface flux \n of SO2",expression(emiso2~(kg~m^-2~s^-1)))
    mmrbc_plot <- plot_species(mmrbc_experiment, pert, "surface concentration \n of BC",expression(mmrbc~(kg~kg-1)))
    mmrso4_plot <- plot_species(mmrso4_experiment, pert, "surface concentration \n of SO4",expression(mmrso4~(kg~kg-1)))
    rlut_plot <- plot_species(rlut_experiment, pert, "upwelling longwave flux \n at TOA",expression(rlut~(W~m-2)))
    rlutcs_plot <- plot_species(rlutcs_experiment, pert, "upwelling clear-sky longwave \n flux at TOA",expression(rlutcs~(W~m-2)))
    rsut_plot <- plot_species(rsut_experiment, pert, "upwelling shortwave flux \n at TOA",expression(rsut~(W~m-2)))
    rsutcs_plot <- plot_species(rsutcs_experiment, pert, "upwelling clear-sky shortwave \n flux at TOA",expression(rsutcs~(W~m-2)))
    rsdt_plot <- plot_species(rsdt_experiment, pert, "incident shortwave flux \n at TOA",expression(rsdt~(W~m-2)))
    net_rad_plot <- plot_species(net_rad, pert, "net radiative flux \n at TOA",expression(rlut~+~rsut~(W~m^-2)))
    net_rad_cs_plot <- plot_species(net_rad_cs, pert, "clear-sky net radiative \n flux at TOA",expression(rlutcs~+~rsutcs~(W~m^-2)))
    so2_plot <- plot_species(so2_experiment, pert, "surface concentration \n of SO2",expression(so2~(kg~kg-1)))
    drybc_plot <- plot_species(drybc_experiment, pert, "dry deposition rate \n of BC",expression(drybc~(kg~m^-2~s^-1)))
    wetbc_plot <- plot_species(wetbc_experiment, pert, "wet deposition rate \n of BC",expression(wetbc~(kg~m^-2~s^-1)))
    tot_bc_plot <- plot_species(tot_bc, pert, "total deposition rate \n of BC", expression(drybc~+~wetbc~(kg~m^-2~s^-1)))
    dryso2_plot <- plot_species(dryso2_experiment, pert, "dry deposition rate \n of SO2", expression(dryso2~(kg~m^-2~s^-1)))
    wetso2_plot <- plot_species(wetso2_experiment, pert, "wet deposition rate \n of SO2", expression(wetso2~(kg~m^-2~s^-1)))
    dryso4_plot <- plot_species(dryso4_experiment, pert, "dry deposition rate \n of SO4", expression(dryso4~(kg~m^-2~s^-1)))
    wetso4_plot <- plot_species(wetso4_experiment, pert, "wet deposition rate \n of SO4", expression(wetso4~(kg~m^-2~s^-1)))
    tot_s_plot <- plot_species(tot_s, pert, "total deposition rate \n of S", expression(atop((dryso2~+~wetso2)/2~+~(dryso4~+~wetso4)/3, (kg~m^-2~s^-1))))
    od550aer_plot <- plot_species(od550aer_experiment, pert, "ambient aerosol optical \n thickness at 550nm", expression(od550aer))
    clt_plot <- plot_species(clt_experiment, pert, "total cloud cover \n percentage", "expression(clt~(%)")
    cltc_plot <- plot_species(cltc_experiment, pert, "convective cloud cover \n percentage", "cltc~(%)")
    imp_cld_plot <- plot_species(imp_cld, pert, "implied cloud response \n at TOA", expression(rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)))
  }
  
  # Function from stack exchange to generate a shared legend
  grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom",
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 7,
                                                                  margin = margin(r = 10, unit = "pt"))))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
      do.call(arrangeGrob, lapply(plots, function(x)
        x + theme(legend.position="none"))),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - 1.5 * lheight, lheight), # the "1.5" adds room for title
      top = textGrob(paste0(scenario, ': absolute difference'), gp = gpar(fontsize = 12)))
  }
  
  final_plot <- grid_arrange_shared_legend(emibc_plot,
                                           emiso2_plot,
                                           mmrbc_plot,
                                           mmrso4_plot,
                                           so2_plot,
                                           rlut_plot,
                                           rsut_plot,
                                           net_rad_plot,
                                           rsdt_plot,
                                           rlutcs_plot,
                                           rsutcs_plot,
                                           net_rad_cs_plot,
                                           imp_cld_plot,
                                           drybc_plot,
                                           wetbc_plot,
                                           tot_bc_plot,
                                           dryso2_plot,
                                           wetso2_plot,
                                           dryso4_plot,
                                           wetso4_plot,
                                           tot_s_plot,
                                           od550aer_plot,
                                           clt_plot,
                                           cltc_plot)
  
  # Print plots
  
  if (sort_by == "region"){
    setwd(paste0('../../../../output/', region, '/timeseries'))
    # To save to file on A4 paper
    ggsave(paste0(region, '_', scenario ,'_diff.pdf'), final_plot, width = 21, height = 29.7, units = "cm")
  }
  if (sort_by == "experiment"){
    setwd(paste0('../../../../output/', pert, '/timeseries'))
    # To save to file on A4 paper
    ggsave(paste0(scenario, '_', pert ,'_diff.pdf'), final_plot, width = 21, height = 29.7, units = "cm")
  }
}
