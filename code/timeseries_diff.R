# ------------------------------------------------------------------------------
# Program Name: timeseries_diff.R
# Authors: Hamza Ahsan
# Date Last Modified: November 15th, 2023
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

# Specify and navigate to the location of Emissions-MIP directory
emi_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Documents/GitHub/Emissions-MIP_Data/Emissions-MIP_Data/")
setwd(paste0(emi_dir))

# Specify what you are sorting by and either the region (i.e., global, land, sea, arctic, NH-land, NH-sea, SH-land, SH-sea) or experiment (i.e., bc-no-season, high-so4, no-so4, reference, so2-at-height, so2-no-season)
#The command line would look like: rscript <rscript>.r <"experiment" or "region"> <specific experiment or region you are sorting by>
sorting <- commandArgs(trailingOnly = TRUE) #pulling region from command line
sort_by <- sorting[1]
if (sort_by == "region"){region <- sorting[2]}
if (sort_by == "experiment"){pert <- sorting[2]}

sort_by <- 'experiment'
pert <- 'shp-60p-red'

# Define colorblind-friendly palette colors and associate with models (in case a
# plot is missing a model, the color scheme will remain consistent)
cbPalette <- c("#c4c4c3", "#4477aa", "#228833", "#66ccee", "#ccbb44","#ee6677", "#aa3377")

model_colors <- c('CESM1' = cbPalette[1], 'GISS-E2.1' = cbPalette[2], 'CAM-ATRAS' = cbPalette[3], 'GEOS' = cbPalette[4], 'NorESM2' = cbPalette[5], 'GFDL-ESM4' = cbPalette[6], 'E3SM' = cbPalette[7])
model_symbols <- c("CESM1" = 15, "GISS-E2.1" = 17, "CAM-ATRAS" = 17, "NorESM2" = 17, "GEOS" = 17, "GFDL-ESM4" = 19, "E3SM" = 15)

# ------------------------------------------------------------------------------
# Iterate over the different perturbation/regional experiments
if (sort_by == "region"){scenarios <- c('shp-30p-red', 'shp-60p-red', 'shp-60p-red-1950',
                                        'shp-atl-shift', 'shp-atl-shift-1950', 'shp-ind-shift',
                                        'shp-ind-shift-1950')}
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
  experiment <- map_dfr(.x = target_filename, .f = ~read.csv(.x) %>% mutate(unit = as.character(unit)))

  # Extract model from file names (fifth segment) and bind to experiment data frame
  models <- sapply(strsplit(target_filename, "[-.]+"),function(x) x[5])
  rep_models <- rep(models, each = 5) # five years
  experiment$model <- rep_models

  # Correct model names
  experiment$model[which(experiment$model == "CESM")] <- "CESM1"
  experiment$model[which(experiment$model == "GISS")] <- "GISS-E2.1"
  experiment$model[which(experiment$model == "CAM5")] <- "CAM-ATRAS"
  experiment$model[which(experiment$model == "GFDL")] <- "GFDL-ESM4"


  # Rearrange data frame by years descending
  experiment <- dplyr::arrange(experiment, year)

  #only runs excluded models if sorting by region
  if (sort_by == "region"){
    #runs through each excluded model pair and filters them out of experiment
    if(nrow(excluded_models) != 0) { #only runs if the data frame is not empty
      for (val in 1:nrow(excluded_models)) {
        summary_long <- filter(summary_long, experiment != excluded_models$Scenario[val] | model != excluded_models$Model[val] | variable != excluded_models$Variable[val])
      }
    }
  }

  # Convert volume mixing ratio to mass mixing ratio by multiplying by molar mass
  # of SO2 and dividing by molar mass of air
  so2_experiment <- dplyr::filter(experiment, variable == 'so2') %>%
    dplyr::mutate(new_value = value * 64.066 / 28.96)

  # Change units from mol/mol to kg/kg
  so2_experiment$unit <- 'kg kg-1'

  # Define other experiments
  emibc_experiment <- dplyr::filter(experiment, variable == 'emibc')
  emiso2_experiment <- dplyr::filter(experiment, variable == 'emiso2')
  mmrbc_experiment <- dplyr::filter(experiment, variable == 'mmrbc')
  mmrso4_experiment <- dplyr::filter(experiment, variable == 'mmrso4')
  rlut_experiment <- dplyr::filter(experiment, variable == 'rlut')
  rlutcs_experiment <- dplyr::filter(experiment, variable == 'rlutcs')
  rsdt_experiment <- dplyr::filter(experiment, variable == 'rsdt')
  rsut_experiment <- dplyr::filter(experiment, variable == 'rsut')
  rsutcs_experiment <- dplyr::filter(experiment, variable == 'rsutcs')
  drybc_experiment <- dplyr::filter(experiment, variable == 'drybc')
  wetbc_experiment <- dplyr::filter(experiment, variable == 'wetbc')
  dryso2_experiment <- dplyr::filter(experiment, variable == 'dryso2')
  wetso2_experiment <- dplyr::filter(experiment, variable == 'wetso2')
  dryso4_experiment <- dplyr::filter(experiment, variable == 'dryso4')
  wetso4_experiment <- dplyr::filter(experiment, variable == 'wetso4')
  clt_experiment <- dplyr::filter(experiment, variable == 'clt')
  cltc_experiment <- dplyr::filter(experiment, variable == 'cltc')
  cl_experiment <- dplyr::filter(experiment, variable == 'cl')
  clivi_experiment <- dplyr::filter(experiment, variable == 'clivi')
  dms_experiment <- dplyr::filter(experiment, variable == 'dms')
  loadso4_experiment <- dplyr::filter(experiment, variable == 'loadso4')
  loadbc_experiment <- dplyr::filter(experiment, variable == 'loadbc')
  od550aer_experiment <- dplyr::filter(experiment, variable == 'od550aer')
  loadso2_experiment <- dplyr::filter(experiment, variable == 'loadso2')

  # Define normal and clear-sky net radiative flux and  (sum of longwave and shortwave radiation)
  net_rad <- dplyr::left_join(rlut_experiment, rsut_experiment, by = c("year", "unit", "model"))
  net_rad <- dplyr::mutate(net_rad, value = value.x + value.y) %>%
    dplyr::select(c(year, unit, model, value))

  net_rad_cs <- dplyr::left_join(rlutcs_experiment, rsutcs_experiment, by = c("year", "unit", "model"))
  net_rad_cs <- dplyr::mutate(net_rad_cs, value = value.x + value.y) %>%
    dplyr::select(c(year, unit, model, value))

  # Define implied cloud response (net - clearsky) as a new variable to plot
  imp_cld <- dplyr::left_join(net_rad, net_rad_cs, by = c("year", "model"))
  imp_cld <- dplyr::mutate(imp_cld, value = value.x - value.y) %>%
    dplyr::select(c(year, model, value))

  # Define total BC deposition rate (sum of dry BC and wet BC deposition)
  tot_bc <- dplyr::left_join(drybc_experiment, wetbc_experiment, by = c("year", "unit", "model"))
  tot_bc <- dplyr::mutate(tot_bc, value = value.x + value.y) %>%
    dplyr::select(c(year, unit, model, value))

  #Define total S deposition rate (sum of dry SO2/SO4 and wet SO2/SO4 deposition)
  dry_s <- dplyr::left_join(dryso2_experiment, dryso4_experiment, by = c("year", "unit", "model"))
  dry_s <- dplyr::mutate(dry_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::select(c(year, unit, model, value))

  wet_s <- dplyr::left_join(wetso2_experiment, wetso4_experiment, by = c("year", "unit", "model"))
  wet_s <- dplyr::mutate(wet_s, value = (32.065/64.066)*value.x + (32.065/96.06)*value.y) %>%
    dplyr::select(c(year, unit, model, value))

  tot_s <- dplyr::left_join(dry_s, wet_s, by = c("year", "unit", "model"))
  tot_s <- dplyr::mutate(tot_s, value = value.x + value.y) %>%
    dplyr::select(c(year, unit, model, value))
  
  # Define total so4 (sum of dry and wet so4)
  tot_so4 <- dplyr::left_join(dryso4_experiment, wetso4_experiment, by = c("year","unit","model"))
  tot_so4 <- dplyr::mutate(tot_so4, value = value.x + value.y) %>%
    #dplyr::mutate(sd = sqrt(sd.x^2 + sd.y^2)) %>%
    dplyr::select(c(year, unit, model, value))
  
  # Define so4 lifetime (loadso4/tot_so4), convert from seconds to days
  so4_lifetime <- dplyr::left_join(loadso4_experiment, tot_so4, by = c("year","unit","model"))
  so4_lifetime <- dplyr::mutate(so4_lifetime, value = (value.x/value.y)/86400) %>%
    #dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
    #dplyr::filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
    #dplyr::filter(value >= 0) %>%
    dplyr::select(c(year, unit, model, value))
  
  # define so2 lifetime (loadso2/emiso2), convert seconds to days
  so2_lifetime <- dplyr::left_join(loadso2_experiment, emiso2_experiment, by = c("year","unit","model"))
  so2_lifetime <- dplyr::mutate(so2_lifetime, value = (value.x/value.y)/86400) %>%
    #dplyr::mutate(sd = value*sqrt((sd.x/value.x)^2 + (sd.y/value.y)^2)) %>%
    #dplyr::filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
    #dplyr::filter(value >= 0) %>%
    dplyr::select(c(year, unit, model, value))


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
    cl_plot <- plot_species(cl_experiment, region, "cloud cover \n percentage", "expression cltc (%)")
    clivi_plot <- plot_species(clivi_experiment, region, "Ice water path", expression(clivi~(kg~m^-2)))
    dms_plot <- plot_species(dms_experiment, region, 'Dimethyl sulphide (DMS) mole fraction', expression(dms~(mol~mol^-1)))
    loadso4_plot <- plot_species(loadso4_experiment, region, "load \n of so4", expression(loadso4~(kg~m^-2)))
    loadbc_plot <- plot_species(loadbc_experiment, region, "load \n of bc", expression(loadbc~(kg~m^-2)))
    imp_cld_plot <- plot_species(imp_cld, region, "implied cloud response \n at TOA", expression(rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)))
    loadso2_plot <- plot_species(loadso2_experiment, region, "load of so2 \n", expression(loadso2~(kg~m^-2)))
    tot_so4_plot <- plot_species(tot_so4, region, 'total SO4 - \n', expression(dryso4~+~wetso4))
    so4_lifetime <- plot_species(so4_lifetime, region, 'SO4 lifetime - \n ', expression(loadso4/(dryso4~+~wetso4)~(days)))
    so2_lifetime <- plot_species(so2_lifetime, region, 'SO2 lifetime - \n', expression(loadso2/emiso2~(days)))
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
    cl_plot <- plot_species(cl_experiment, pert, "cloud cover \n percentage", "expression cltc (%)")
    clivi_plot <- plot_species(clivi_experiment, pert, "Ice water path", expression(clivi~(kg~m^-2)))
    dms_plot <- plot_species(dms_experiment, pert, 'Dimethyl sulphide (DMS) mole fraction', expression(dms~(mol~mol^-1)))
    loadso4_plot <- plot_species(loadso4_experiment, pert, "load \n of so4", expression(loadso4~(kg~m^-2)))
    loadbc_plot <- plot_species(loadbc_experiment, pert, "load \n of bc", expression(loadbc~(kg~m^-2)))
    imp_cld_plot <- plot_species(imp_cld, pert, "implied cloud response \n at TOA", expression(rlut~+~rsut~-~rlutcs~-~rsutcs~(W~m^-2)))
    loadso2_plot <- plot_species(loadso2_experiment, pert, "load of so2 \n", expression(loadso2~(kg~m^-2)))
    #tot_so4_plot <- plot_species(tot_so4, pert, value, 'total SO4 - \n', expression(dryso4~+~wetso4~(kg~m^-2~s^-1)))
    #so4_lifetime_plot <- plot_species(so4_lifetime, pert, value, 'SO4 lifetime - \n ', expression(loadso4/(dryso4~+~wetso4)~(days)))
    #so2_lifetime_plot <- plot_species(so2_lifetime, pert, value, 'SO2 lifetime - \n', expression(loadso2/emiso2~(days)))
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
                                           cl_plot,
                                           clivi_plot,
                                           dms_plot,
                                           od550aer_plot,
                                           loadso4_plot,
                                           loadbc_plot,
                                           loadso2_plot)#,
                                           #so4_lifetime_plot,
                                           #so2_lifetime_plot)

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
