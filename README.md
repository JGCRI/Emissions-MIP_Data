# Emissions-MIP_Data
Emissions Model Intercomparison Project

Data: [![DOI](https://zenodo.org/badge/362623387.svg)](https://zenodo.org/badge/latestdoi/362623387)

Publication: Ahsan, H., Wang, H., Wu, J., Wu, M., Smith, S. J., Bauer, S., Suchyta, H., Olivié, D., Myhre, G., Matsui, H., Bian, H., Lamarque, J.-F., Carslaw, K., Horowitz, L., Regayre, L., Chin, M., Schulz, M., Skeie, R. B., Takemura, T., and Naik, V.: The Emissions Model Intercomparison Project (Emissions-MIP): quantifying model sensitivity to emission characteristics, Atmos. Chem. Phys., 23, 14779–14799, https://doi.org/10.5194/acp-23-14779-2023, 2023.

# Description
This repository contains the timeseries data generated by ESMValTool for each climate model as well as the scripts used to create the timeseries and summary plots. A current copy of the plots is stored here as well.

The climate models currently included are CAM5-ATRAS, CESM1, CESM2-WACCM, E3SMv1, GEOS-i33p2, GFDL-AM4.1, GISS ModelE2.1, HadGEM-UKCA, MIROC-SPRINTARS, NorESM, OsloCTM3. This repository will be updated in the future as more models are made available.

## Directories
### code
This directory contains the R scripts used to generate the plots. Before running these scripts, remember to specify the location of the "Emissions-MIP" directory to its location on your computer (*i.e.*, variable `emi_dir`). Also, specify the region to generate plots for.

### input
The perturbation and reference data are included for each region (*e.g.*, global, land, sea, etc.). These are csv files of annual, averaged (over respective region) data generated by ESMValTool.

### output
Two types of plots are generated - timeseries and summary plots. The timeseries plots include the reference (base) case, the absolute difference between the perturbations and the reference case, and the percent difference between the perturbations and the reference case. The summary plots exhibit the absolute and percent differences averaged over all years for each perturbation.

## Instructions

### summary/timeseries plot scripts

  These are run from the command line so you need to run Rscript <name of script> <experiment or region (no quotes)> <the specific experiment or region you are running (also no quotes)>
  
  Make sure that the experiments you are running are consistent with what is in the var_master_list.csv file
  
  You can also run straight from the R script file, just update the sort_by, region, and exper variables

### fixed_axes.csv
	
For the first column, enter in the first row which method you want to use in order to determine y axes 
	in the ouput file. The choices are:

	1. group
	If set to group then you must create a group name in the first row of a column and add the variables in
	that group below it. The script will then make sure each variable in a group is plotted with the same y 
	axes. NOTE: do not put a variable in more than one group as it will default to the axes of the last group
	it is listed in. You can have as many groups as you want

	2. fixed
	If set to fixed then the script will determine the min and max values across all regions or experiments
	(depending on what you are sorting each run by) and set the output y axes to those values

	3. neither
	The script runs as normal

**Note:** make sure the order of the variables under 'vars' in var_master_list.csv and 'Variable' in variables.csv are are in the same order. All stand-alone variables (ie mmrbc, emibc...) should be listed before any combined variables (ie net_rad, tot_s...) in variables.csv

### combined_variables.csv
	
The column headers are the name of the ocmbined variable and the values below are the variables joined in order to make the combined variable. You can only join two at a time, so if you wanted to combine four variables you would need some intermediate steps and put them to the left of the final column (see for example tot_s, which is a combination of dry_s and wet_s. Make sure to add these intermediates to variables.csv too).


*Please alert us if you come across any issues or if you would like to make recommendations for improvements.*
