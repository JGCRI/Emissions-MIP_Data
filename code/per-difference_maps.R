# ------------------------------------------------------------------------------
# Program Name: per-difference_maps.R
# Authors: Ciara Donegan
# Date Last Modified: April 15, 2024
# Program Purpose: Produces maps of the difference between perturbation datasets
# and reference datasets for Emissions-MIP.
# Input Files: ~Emissions-MIP/input/
# Output Files: ~Emissions-MIP/output/
# TODO:
# ------------------------------------------------------------------------------

# Load required libraries
library(ncdf4)
library(lattice)
library(RColorBrewer)
library(raster)
library(ggplot2)
library(terra)
library(rasterVis)
library(classInt)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

# set path to netCDF files - replace with your path to files
file_dir <- paste0("C:/Users/done231/OneDrive - PNNL/Desktop/difference_plot_files")
setwd(file_dir)

# select scenario to generate maps for
scenario <- "shp-60p-red"

# read in shapefile for coastline
coast_outline <- shapefile("ne_110m_coastline.shp")

# get data from ncdf files
get_ncdf_data <- function(variable,model) {
  filepath <- paste0(file_dir,"/input/",scenario,"/",model)
  
  # Get control file
  if (model == "CAM5ATRAS" | model == "GEOS" | model == "GFDL" | model == "NorESM2") {
    file_control <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                               pattern=paste0("base.+",variable,".nc"))
  }
  
  if (model == "CESM1" | model == "E3SM") {
    file_control <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                               pattern=paste0("ref.+",variable,".nc"))
  }
  
  if (model == "GISS" & scenario != "shp-60p-red-1950" & scenario != "shp-atl-shift-1950" & scenario != "shp-ind-shift-1950") {
    file_control <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                               pattern=paste0("ref.+",variable,".nc"))
  }
  
  if (model == "GISS" & (scenario == "shp-60p-red-1950" | scenario == "shp-atl-shift-1950" | scenario == "shp-ind-shift-1950")) {
    file_control <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                               pattern=paste0("BW1950.+",variable,".nc"))
  }
  
  # Get experiment file
  if (model != "GISS" & scenario == "shp-30p-red") {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                                           pattern=paste0("30.+",variable,".nc"))
  }
  
  if (model != "GISS" & scenario == "shp-60p-red" | scenario == "shp-60p-red-1950") {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                                           pattern=paste0("60.+",variable,".nc"))
  }
  
  
  if (model != "GISS" & (scenario == "shp-atl-shift" | scenario == "shp-atl-shift-1950")) {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                                           pattern=paste0("atl.+",variable,".nc"))
  }
  
  if (model != "GISS" & (scenario == "shp-ind-shift" | scenario == "shp-ind-shift-1950")) {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                                           pattern=paste0("ind.+",variable,".nc"))
  }
  
  if (model == "GISS" & (scenario == "shp-30p-red" | scenario == "shp-60p-red" | scenario == "shp-60p-red-1950")) {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                           pattern=paste0("dot.+",variable,".nc"))
  }
  
  if (model == "GISS" & scenario != "shp-30p-red" & scenario != "shp-60p-red" & scenario != "shp-60p-red-1950") {
    file_exp <- list.files(path=paste0(file_dir,"/input/",scenario,"/",model),
                           pattern=paste0("SHFT.+",variable,".nc"))
  }
  
  ## get control dataset
  nc_control <- nc_open(paste0(filepath,"/",file_control))
  
  # get coordinate vars
  lon <- ncvar_get(nc_control,"lon") # 0 to 360, rather than -180 to 180
  lon <- lon - 180
  nlon <- dim(lon)
  
  lat <- ncvar_get(nc_control,"lat") # -90 to 90
  nlat <- dim(lat)
  
  year <- ncvar_get(nc_control,"year") # values are averaged across time range
  
  # read in data from chosen variable
  control.array <- ncvar_get(nc_control,variable)
  fillvalue <- ncatt_get(nc_control,variable,"_FillValue") # get value used for missing data
  
  # close netCDF file
  nc_close(nc_control)
  
  # get ncdf coordinates to match shapefile (centered on 0 deg longitude)
  east <- control.array[181:360,]
  west <- control.array[1:180,]
  control.shifted <- rbind(east,west)
  
  # raster
  r <- raster(t(control.shifted),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  r <- flip(r,direction='y')
  
  #  convert to a format usable by ggplot
  r.control <- as.data.frame(r,xy=TRUE)
  
  # remove fill value
  r.control[r.control>1e36] <- NA
  
  ## get experiment dataset
  nc_exp <- nc_open(paste0(filepath,"/",file_exp))
  
  # get coordinate vars
  lon <- ncvar_get(nc_exp,"lon") # 0 to 360, rather than -180 to 180
  lon <- lon - 180
  nlon <- dim(lon)
  
  lat <- ncvar_get(nc_exp,"lat") # -90 to 90
  nlat <- dim(lat)
  
  year <- ncvar_get(nc_exp,"year") # values are averaged across time range
  
  # read in data from chosen variable
  exp.array <- ncvar_get(nc_exp,variable)
  fillvalue <- ncatt_get(nc_exp,variable,"_FillValue") # get value used for missing data
  
  # close netCDF file
  nc_close(nc_exp)
  
  # get ncdf coordinates to match shapefile (centered on 0 deg longitude)
  east <- exp.array[181:360,]
  west <- exp.array[1:180,]
  exp.shifted <- rbind(east,west)
  
  # raster
  r <- raster(t(exp.shifted),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  r <- flip(r,direction='y')
  
  #  convert to a format usable by ggplot
  r.exp <- as.data.frame(r,xy=TRUE)
  
  # remove fill value
  r.exp[r.exp>1e36] <- NA
  
  # Get percent difference dataset
  r.df <- (r.exp - r.control)/r.control*100
  r.df$x <- r.exp$x
  r.df$y <- r.exp$y
  
  ## Invert sign of forcing variables to match convention (positive = warming effect)
  #if (variable == "rlut" | variable == "rsut" | variable == "rlutcs" | variable == "rsutcs") {
  #  r.df$layer <- r.df$layer*-1
  #}
  
  return(r.df)
}

# get plot from dataframe
get_plot <- function(variable,model,df,at,title) {
  
  # levelplot from lattice
  myTheme <- rasterTheme(region = rev(brewer.pal(11,'RdBu')))
  jpeg_map <- levelplot(layer~x*y, data=df,
                        par.settings=myTheme,
                        at=at,
                        margin=F,
                        xlab=NULL,
                        ylab=NULL,
                        zlab="test",
                        main=title) +
    latticeExtra::layer(sp.lines(coast_outline, col="black", lwd=0.5))
  
  return(jpeg_map)
}

## Get dataframes

# CAM5ATRAS
clt_CAM5ATRAS <- get_ncdf_data("clt","CAM5ATRAS")
loadbc_CAM5ATRAS <- get_ncdf_data("loadbc","CAM5ATRAS")
loadso2_CAM5ATRAS <- get_ncdf_data("loadso2","CAM5ATRAS")
loadso4_CAM5ATRAS <- get_ncdf_data("loadso4","CAM5ATRAS")
rlut_CAM5ATRAS <- get_ncdf_data("rlut","CAM5ATRAS")
rlutcs_CAM5ATRAS <- get_ncdf_data("rlutcs","CAM5ATRAS")
rsut_CAM5ATRAS <- get_ncdf_data("rsut","CAM5ATRAS")
rsutcs_CAM5ATRAS <- get_ncdf_data("rsutcs","CAM5ATRAS")

# CESM1
clt_CESM1 <- get_ncdf_data("clt","CESM1")
loadbc_CESM1 <- get_ncdf_data("loadbc","CESM1")
loadso2_CESM1 <- get_ncdf_data("loadso2","CESM1")
loadso4_CESM1 <- get_ncdf_data("loadso4","CESM1")
rlut_CESM1 <- get_ncdf_data("rlut","CESM1")
rlutcs_CESM1 <- get_ncdf_data("rlutcs","CESM1")
rsut_CESM1 <- get_ncdf_data("rsut","CESM1")
rsutcs_CESM1 <- get_ncdf_data("rsutcs","CESM1")

# E3SM
clt_E3SM <- get_ncdf_data("clt","E3SM")
if (scenario != "shp-atl-shift") {
  loadbc_E3SM <- get_ncdf_data("loadbc","E3SM")
  loadso2_E3SM <- get_ncdf_data("loadso2","E3SM")
  loadso4_E3SM <- get_ncdf_data("loadso4","E3SM")
}
rlut_E3SM <- get_ncdf_data("rlut","E3SM")
rlutcs_E3SM <- get_ncdf_data("rlutcs","E3SM")
rsut_E3SM <- get_ncdf_data("rsut","E3SM")
rsutcs_E3SM <- get_ncdf_data("rsutcs","E3SM")

# GEOS
clt_GEOS <- get_ncdf_data("clt","GEOS")
loadbc_GEOS <- get_ncdf_data("loadbc","GEOS")
loadso2_GEOS <- get_ncdf_data("loadso2","GEOS")
loadso4_GEOS <- get_ncdf_data("loadso4","GEOS")
rlut_GEOS <- get_ncdf_data("rlut","GEOS")
rlutcs_GEOS <- get_ncdf_data("rlutcs","GEOS")
rsut_GEOS <- get_ncdf_data("rsut","GEOS")
rsutcs_GEOS <- get_ncdf_data("rsutcs","GEOS")

# GFDL
if (scenario == "shp-60p-red" || scenario == "shp-atl-shift" || scenario == "shp-ind-shift") {
  loadso2_GFDL <- get_ncdf_data("loadso2","GFDL")

}
if (scenario != "shp-ind-shift-1950" & scenario != "shp-30p-red") {
  loadbc_GFDL <- get_ncdf_data("loadbc","GFDL")
  loadso4_GFDL <- get_ncdf_data("loadso4","GFDL")
}
rlut_GFDL <- get_ncdf_data("rlut","GFDL")
rlutcs_GFDL <- get_ncdf_data("rlutcs","GFDL")
rsut_GFDL <- get_ncdf_data("rsut","GFDL")
rsutcs_GFDL <- get_ncdf_data("rsutcs","GFDL")

# GISS
clt_GISS <- get_ncdf_data("clt","GISS")
loadbc_GISS <- get_ncdf_data("loadbc","GISS")
loadso2_GISS <- get_ncdf_data("loadso2","GISS")
loadso4_GISS <- get_ncdf_data("loadso4","GISS")
rlut_GISS <- get_ncdf_data("rlut","GISS")
rlutcs_GISS <- get_ncdf_data("rlutcs","GISS")
rsut_GISS <- get_ncdf_data("rsut","GISS")
rsutcs_GISS <- get_ncdf_data("rsutcs","GISS")

# NorESM2
clt_NorESM2 <- get_ncdf_data("clt","NorESM2")
loadbc_NorESM2 <- get_ncdf_data("loadbc","NorESM2")
loadso2_NorESM2 <- get_ncdf_data("loadso2","NorESM2")
loadso4_NorESM2 <- get_ncdf_data("loadso4","NorESM2")
rlut_NorESM2 <- get_ncdf_data("rlut","NorESM2")
rlutcs_NorESM2 <- get_ncdf_data("rlutcs","NorESM2")
rsut_NorESM2 <- get_ncdf_data("rsut","NorESM2")
rsutcs_NorESM2 <- get_ncdf_data("rsutcs","NorESM2")

## Get breaks for plots, standard across all models
# clt
all_clt <- rbind(clt_CAM5ATRAS,clt_CESM1,clt_E3SM,clt_GEOS,clt_GISS,clt_NorESM2)
all_clt <- drop_na(all_clt)
bound_clt <- max(abs(all_clt$layer))
clt.at <- seq(-bound_clt,bound_clt,length.out=25)

# loadbc
if (scenario == "shp-30p-red" || scenario == "shp-ind-shift-1950") {
  all_loadbc <- rbind(loadbc_CAM5ATRAS,loadbc_CESM1,loadbc_E3SM,loadbc_GEOS,loadbc_GISS,loadbc_NorESM2)
  all_loadbc <- drop_na(all_loadbc)
  bound_loadbc <- max(abs(all_loadbc$layer))
  loadbc.at <- seq(-bound_loadbc,bound_loadbc,length.out=25)
} else if (scenario == "shp-atl-shift") {
  all_loadbc <- rbind(loadbc_CAM5ATRAS,loadbc_CESM1,loadbc_GEOS,loadbc_GFDL,loadbc_GISS,loadbc_NorESM2)
  all_loadbc <- drop_na(all_loadbc)
  bound_loadbc <- max(abs(all_loadbc$layer))
  loadbc.at <- seq(-bound_loadbc,bound_loadbc,length.out=25)
} else {
  all_loadbc <- rbind(loadbc_CAM5ATRAS,loadbc_CESM1,loadbc_E3SM,loadbc_GEOS,loadbc_GFDL,loadbc_GISS,loadbc_NorESM2)
  all_loadbc <- drop_na(all_loadbc)
  bound_loadbc <- max(abs(all_loadbc$layer))
  loadbc.at <- seq(-bound_loadbc,bound_loadbc,length.out=25)
}

# loadso2
if (scenario == "shp-30p-red" || scenario == "shp-60p-red-1950" || scenario == "shp-atl-shift-1950" || scenario == "shp-ind-shift-1950") {
  all_loadso2 <- rbind(loadso2_CAM5ATRAS,loadso2_CESM1,loadso2_E3SM,loadso2_GEOS,loadso2_GISS,loadso2_NorESM2)
  all_loadso2 <- drop_na(all_loadso2)
  bound_loadso2 <- max(abs(all_loadso2$layer))
  loadso2.at <- seq(-bound_loadso2,bound_loadso2,length.out=25)
} else if (scenario == "shp-atl-shift") {
  all_loadso2 <- rbind(loadso2_CAM5ATRAS,loadso2_CESM1,loadso2_GEOS,loadso2_GFDL,loadso2_GISS,loadso2_NorESM2)
  all_loadso2 <- drop_na(all_loadso2)
  bound_loadso2 <- max(abs(all_loadso2$layer))
  loadso2.at <- seq(-bound_loadso2,bound_loadso2,length.out=25)
} else {
  all_loadso2 <- rbind(loadso2_CAM5ATRAS,loadso2_CESM1,loadso2_E3SM,loadso2_GEOS,loadso2_GFDL,loadso2_GISS,loadso2_NorESM2)
  all_loadso2 <- drop_na(all_loadso2)
  bound_loadso2 <- max(abs(all_loadso2$layer))
  loadso2.at <- seq(-bound_loadso2,bound_loadso2,length.out=25)
}

# loadso4
if (scenario == "shp-30p-red" || scenario == "shp-ind-shift-1950") {
  all_loadso4 <- rbind(loadso4_CAM5ATRAS,loadso4_CESM1,loadso4_E3SM,loadso4_GEOS,loadso4_GISS,loadso4_NorESM2)
  all_loadso4 <- drop_na(all_loadso4)
  bound_loadso4 <- max(abs(all_loadso4$layer))
  loadso4.at <- seq(-bound_loadso4,bound_loadso4,length.out=25)
} else if (scenario == "shp-atl-shift") {
  all_loadso4 <- rbind(loadso4_CAM5ATRAS,loadso4_CESM1,loadso4_GEOS,loadso4_GFDL,loadso4_GISS,loadso4_NorESM2)
  all_loadso4 <- drop_na(all_loadso4)
  bound_loadso4 <- max(abs(all_loadso4$layer))
  loadso4.at <- seq(-bound_loadso4,bound_loadso4,length.out=25)
} else {
  all_loadso4 <- rbind(loadso4_CAM5ATRAS,loadso4_CESM1,loadso4_E3SM,loadso4_GEOS,loadso4_GFDL,loadso4_GISS,loadso4_NorESM2)
  all_loadso4 <- drop_na(all_loadso4)
  bound_loadso4 <- max(abs(all_loadso4$layer))
  loadso4.at <- seq(-bound_loadso4,bound_loadso4,length.out=25)
}

# rlut
all_rlut <- rbind(rlut_CAM5ATRAS,rlut_CESM1,rlut_E3SM,rlut_GEOS,rlut_GFDL,rlut_GISS,rlut_NorESM2)
all_rlut <- drop_na(all_rlut)
bound_rlut <- max(abs(all_rlut$layer))
rlut.at <- seq(-bound_rlut,bound_rlut,length.out=25)

# rlutcs
all_rlutcs <- rbind(rlutcs_CAM5ATRAS,rlutcs_CESM1,rlutcs_E3SM,rlutcs_GEOS,rlutcs_GFDL,rlutcs_GISS,rlutcs_NorESM2)
all_rlutcs <- drop_na(all_rlutcs)
bound_rlutcs <- max(abs(all_rlutcs$layer))
rlutcs.at <- seq(-bound_rlutcs,bound_rlutcs,length.out=25)

# rsut
#if (scenario == "shp-60p-red-1950") {
#  all_rsut <- rbind(rsut_CAM5ATRAS,rsut_E3SM,rsut_GEOS,rsut_GFDL,rsut_GISS,rsut_NorESM2)
#} else {
all_rsut <- rbind(rsut_CAM5ATRAS,rsut_CESM1,rsut_E3SM,rsut_GEOS,rsut_GFDL,rsut_GISS,rsut_NorESM2)
#}
all_rsut <- drop_na(all_rsut)
bound_rsut <- max(abs(all_rsut$layer))
rsut.at <- seq(-bound_rsut,bound_rsut,length.out=25)

# rsutcs - remove CESM1 for ind shift, outlier
#if (scenario == "shp-ind-shift" | scenario == "shp-ind-shift-1950") {
#  all_rsutcs <- rbind(rsutcs_CAM5ATRAS,rsutcs_E3SM,rsutcs_GEOS,rsutcs_GFDL,rsutcs_GISS,rsutcs_NorESM2)
#} else {
all_rsutcs <- rbind(rsutcs_CAM5ATRAS,rsutcs_CESM1,rsutcs_E3SM,rsutcs_GEOS,rsutcs_GFDL,rsutcs_GISS,rsutcs_NorESM2)
#}
all_rsutcs <- drop_na(all_rsutcs)
bound_rsutcs <- max(abs(all_rsutcs$layer))
rsutcs.at <- seq(-bound_rsutcs,bound_rsutcs,length.out=25)

# standardize flux scales
flux.at.list <- list(rlut.at,rlutcs.at,rsut.at,rsutcs.at)
if (scenario == "shp-30p-red") {
  flux.at.shp_30p_red <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-60p-red") {
  flux.at.shp_60p_red <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-60p-red-1950") {
  flux.at.shp_60p_red_1950 <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-atl-shift") {
  flux.at.shp_atl_shift <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-atl-shift-1950") {
  flux.at.shp_atl_shift_1950 <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-ind-shift") {
  flux.at.shp_ind_shift <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}
if (scenario == "shp-ind-shift-1950") {
  flux.at.shp_ind_shift_1950 <- flux.at.list[[which.max(sapply(flux.at.list,max))]]
}


## standarize flux scales between scenarios
all.flux.at.list <- list(flux.at.shp_30p_red,flux.at.shp_60p_red,flux.at.shp_60p_red_1950,
                         flux.at.shp_atl_shift,flux.at.shp_ind_shift,flux.at.shp_atl_shift_1950,
                         flux.at.shp_ind_shift_1950)
#flux.at <- all.flux.at.list[[which.max(sapply(flux.at.list,max))]]
flux.at <- flux.at.shp_atl_shift # above isn't working for some reason, but only for 60p 1950??

## Get plots
# CAM5ATRAS
clt_CAM5ATRAS_plot <- get_plot("clt","CAM5ATRAS",clt_CAM5ATRAS,clt.at,"Cloud Area Fraction (% Diff) - CAM-ATRAS")
loadbc_CAM5ATRAS_plot <- get_plot("loadbc","CAM5ATRAS",loadbc_CAM5ATRAS,loadbc.at,"Load of Black Carbon (% Diff) - CAM-ATRAS")
loadso2_CAM5ATRAS_plot <- get_plot("loadso2","CAM5ATRAS",loadso2_CAM5ATRAS,loadso2.at,"Load of SO2 (% Diff) - CAM-ATRAS")
loadso4_CAM5ATRAS_plot <- get_plot("loadso4","CAM5ATRAS",loadso4_CAM5ATRAS,loadso4.at,"Load of SO4 (% Diff) - CAM-ATRAS")
rlut_CAM5ATRAS_plot <- get_plot("rlut","CAM5ATRAS",rlut_CAM5ATRAS,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - CAM-ATRAS")
rlutcs_CAM5ATRAS_plot <- get_plot("rlutcs","CAM5ATRAS",rlutcs_CAM5ATRAS,flux.at,"Upwelling Clear-Sky Longwave \n Radiation at TOA (% Diff) - CAM-ATRAS")
rsut_CAM5ATRAS_plot <- get_plot("rsut","CAM5ATRAS",rsut_CAM5ATRAS,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - CAM-ATRAS")
rsutcs_CAM5ATRAS_plot <- get_plot("rsutcs","CAM5ATRAS",rsut_CAM5ATRAS,flux.at,"Upwelling Clear-Sky Longwave \n Radiation at TOA (% Diff) - CAM-ATRAS")

# CESM1
clt_CESM1_plot <- get_plot("clt","CESM1",clt_CESM1,clt.at,"Cloud Area Fraction (% Diff) - CESM1")
loadbc_CESM1_plot <- get_plot("loadbc","CESM1",loadbc_CESM1,loadbc.at,"Load of Black Carbon (% Diff) - CESM1")
loadso2_CESM1_plot <- get_plot("loadso2","CESM1",loadso2_CESM1,loadso2.at,"Load of SO2 (% Diff) - CESM1")
loadso4_CESM1_plot <- get_plot("loadso4","CESM1",loadso4_CESM1,loadso4.at,"Load of SO4 (% Diff) - CESM1")
rlut_CESM1_plot <- get_plot("rlut","CESM1",rlut_CESM1,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - CESM1")
rlutcs_CESM1_plot <- get_plot("rlutcs","CESM1",rlutcs_CESM1,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - CESM1")
rsut_CESM1_plot <- get_plot("rsut","CESM1",rsut_CESM1,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - CESM1")
rsutcs_CESM1_plot <- get_plot("rsutcs","CESM1",rsut_CESM1,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - CESM1")

# E3SM
clt_E3SM_plot <- get_plot("clt","E3SM",clt_E3SM,clt.at,"Cloud Area Fraction (% Diff) - E3SM")
if (scenario != "shp-atl-shift") {
  loadbc_E3SM_plot <- get_plot("loadbc","E3SM",loadbc_E3SM,loadbc.at,"Load of Black Carbon (% Diff) - E3SM")
  loadso2_E3SM_plot <- get_plot("loadso2","E3SM",loadso2_E3SM,loadso2.at,"Load of SO2 (% Diff) - E3SM")
  loadso4_E3SM_plot <- get_plot("loadso4","E3SM",loadso4_E3SM,loadso4.at,"Load of SO4 (% Diff) - E3SM")
}
rlut_E3SM_plot <- get_plot("rlut","E3SM",rlut_E3SM,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - E3SM")
rlutcs_E3SM_plot <- get_plot("rlutcs","E3SM",rlutcs_E3SM,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - E3SM")
rsut_E3SM_plot <- get_plot("rsut","E3SM",rsut_E3SM,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - E3SM")
rsutcs_E3SM_plot <- get_plot("rsutcs","E3SM",rsut_E3SM,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - E3SM")

# GEOS
clt_GEOS_plot <- get_plot("clt","GEOS",clt_GEOS,clt.at,"Cloud Area Fraction (% Diff) - GEOS")
loadbc_GEOS_plot <- get_plot("loadbc","GEOS",loadbc_GEOS,loadbc.at,"Load of Black Carbon (% Diff) - GEOS")
loadso2_GEOS_plot <- get_plot("loadso2","GEOS",loadso2_GEOS,loadso2.at,"Load of SO2 (% Diff) - GEOS")
loadso4_GEOS_plot <- get_plot("loadso4","GEOS",loadso4_GEOS,loadso4.at,"Load of SO4 (% Diff) - GEOS")
rlut_GEOS_plot <- get_plot("rlut","GEOS",rlut_GEOS,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - GEOS")
rlutcs_GEOS_plot <- get_plot("rlutcs","GEOS",rlutcs_GEOS,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - GEOS")
rsut_GEOS_plot <- get_plot("rsut","GEOS",rsut_GEOS,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - GEOS")
rsutcs_GEOS_plot <- get_plot("rsutcs","GEOS",rsut_GEOS,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - GEOS")

# GISS
clt_GISS_plot <- get_plot("clt","GISS",clt_GISS,clt.at,"Cloud Area Fraction (% Diff) - GISS-E2.1")
loadbc_GISS_plot <- get_plot("loadbc","GISS",loadbc_GISS,loadbc.at,"Load of Black Carbon (% Diff) - GISS-E2.1")
loadso2_GISS_plot <- get_plot("loadso2","GISS",loadso2_GISS,loadso2.at,"Load of SO2 (% Diff) - GISS-E2.1")
loadso4_GISS_plot <- get_plot("loadso4","GISS",loadso4_GISS,loadso4.at,"Load of SO4 (% Diff) - GISS-E2.1")
rlut_GISS_plot <- get_plot("rlut","GISS",rlut_GISS,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - GISS-E2.1")
rlutcs_GISS_plot <- get_plot("rlutcs","GISS",rlutcs_GISS,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - GISS-E2.1")
rsut_GISS_plot <- get_plot("rsut","GISS",rsut_GISS,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - GISS-E2.1")
rsutcs_GISS_plot <- get_plot("rsutcs","GISS",rsut_GISS,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - GISS-E2.1")

# GFDL
if (scenario == "shp-60p-red" || scenario == "shp-atl-shift" || scenario == "shp-ind-shift") {
  loadso2_GFDL_plot <- get_plot("loadso2","GFDL",loadso2_GFDL,loadso2.at,"Load of SO2 (% Diff) - GFDL-ESM4")
}
if (scenario != "shp-ind-shift-1950" & scenario != "shp-30p-red") {
  loadbc_GFDL_plot <- get_plot("loadbc","GFDL",loadbc_GFDL,loadbc.at,"Load of Black Carbon (% Diff) - GFDL-ESM4")
  loadso4_GFDL_plot <- get_plot("loadso4","GFDL",loadso4_GFDL,loadso4.at,"Load of SO4 (% Diff) - GFDL-ESM4")
}
rlut_GFDL_plot <- get_plot("rlut","GFDL",rlut_GFDL,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - GFDL-ESM4")
rlutcs_GFDL_plot <- get_plot("rlutcs","GFDL",rlutcs_GFDL,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - GFDL-ESM4")
rsut_GFDL_plot <- get_plot("rsut","GFDL",rsut_GFDL,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - GFDL-ESM4")
rsutcs_GFDL_plot <- get_plot("rsutcs","GFDL",rsut_GFDL,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - GFDL-ESM4")

# NorESM2
clt_NorESM2_plot <- get_plot("clt","NorESM2",clt_NorESM2,clt.at,"Cloud Area Fraction (% Diff) - NorESM2")
loadbc_NorESM2_plot <- get_plot("loadbc","NorESM2",loadbc_NorESM2,loadbc.at,"Load of Black Carbon (% Diff) - NorESM2")
loadso2_NorESM2_plot <- get_plot("loadso2","NorESM2",loadso2_NorESM2,loadso2.at,"Load of SO2 (% Diff) - NorESM2")
loadso4_NorESM2_plot <- get_plot("loadso4","NorESM2",loadso4_NorESM2,loadso4.at,"Load of SO4 (% Diff) - NorESM2")
rlut_NorESM2_plot <- get_plot("rlut","NorESM2",rlut_NorESM2,flux.at,"Upwelling Longwave Radiation \n at TOA (% Diff) - NorESM2")
rlutcs_NorESM2_plot <- get_plot("rlutcs","NorESM2",rlutcs_NorESM2,flux.at,"Clear-Sky Upwelling Longwave \n Radiation at TOA (% Diff) - NorESM2")
rsut_NorESM2_plot <- get_plot("rsut","NorESM2",rsut_NorESM2,flux.at,"Upwelling Shortwave Radiation \n at TOA (% Diff) - NorESM2")
rsutcs_NorESM2_plot <- get_plot("rsutcs","NorESM2",rsut_NorESM2,flux.at,"Clear-Sky Upwelling Shortwave \n Radiation at TOA (% Diff) - NorESM2")

# save plot function
save_plot <- function(plot) {
  setwd(paste0(file_dir,"/output/scenario/",scenario,"/",model))
  png(filename=paste0(deparse(substitute(plot)),"_per-diff.png"),
      width=844,height=620,units="px")
  print(plot)
  dev.off()
  
  # back to previous wd
  setwd(file_dir)
}

## save plots
# save CAM5ATRAS plots
model <- "CAM5ATRAS"
save_plot(clt_CAM5ATRAS_plot)
save_plot(loadbc_CAM5ATRAS_plot)
save_plot(loadso2_CAM5ATRAS_plot)
save_plot(loadso4_CAM5ATRAS_plot)
save_plot(rlut_CAM5ATRAS_plot)
save_plot(rlutcs_CAM5ATRAS_plot)
save_plot(rsut_CAM5ATRAS_plot)
save_plot(rsutcs_CAM5ATRAS_plot)

# save CESM1 plots
model <- "CESM1"
save_plot(clt_CESM1_plot)
save_plot(loadbc_CESM1_plot)
save_plot(loadso2_CESM1_plot)
save_plot(loadso4_CESM1_plot)
save_plot(rlut_CESM1_plot)
save_plot(rlutcs_CESM1_plot)
save_plot(rsut_CESM1_plot)
save_plot(rsutcs_CESM1_plot)

# save E3SM plots
model <- "E3SM"
save_plot(clt_E3SM_plot)
if (scenario != "shp-atl-shift") {
  save_plot(loadbc_E3SM_plot)
  save_plot(loadso2_E3SM_plot)
  save_plot(loadso4_E3SM_plot)
}
save_plot(rlut_E3SM_plot)
save_plot(rlutcs_E3SM_plot)
save_plot(rsut_E3SM_plot)
save_plot(rsutcs_E3SM_plot)

# save GEOS plots
model <- "GEOS"
save_plot(clt_GEOS_plot)
save_plot(loadbc_GEOS_plot)
save_plot(loadso2_GEOS_plot)
save_plot(loadso4_GEOS_plot)
save_plot(rlut_GEOS_plot)
save_plot(rlutcs_GEOS_plot)
save_plot(rsut_GEOS_plot)
save_plot(rsutcs_GEOS_plot)

# save GISS plots
model <- "GISS"
save_plot(clt_GISS_plot)
save_plot(loadbc_GISS_plot)
save_plot(loadso2_GISS_plot)
save_plot(loadso4_GISS_plot)
save_plot(rlut_GISS_plot)
save_plot(rlutcs_GISS_plot)
save_plot(rsut_GISS_plot)
save_plot(rsutcs_GISS_plot)

# save GFDL plots
model <- "GFDL"
if (scenario == "shp-60p-red" || scenario == "shp-atl-shift" || scenario == "shp-ind-shift") {
  save_plot(loadso2_GFDL_plot)
}
if (scenario != "shp-ind-shift-1950" & scenario != "shp-30p-red") {
  save_plot(loadbc_GFDL_plot)
  save_plot(loadso4_GFDL_plot)
}
save_plot(rlut_GFDL_plot)
save_plot(rlutcs_GFDL_plot)
save_plot(rsut_GFDL_plot)
save_plot(rsutcs_GFDL_plot)

# save NorESM2 plots
model <- "NorESM2"
save_plot(clt_NorESM2_plot)
save_plot(loadbc_NorESM2_plot)
save_plot(loadso2_NorESM2_plot)
save_plot(loadso4_NorESM2_plot)
save_plot(rlut_NorESM2_plot)
save_plot(rlutcs_NorESM2_plot)
save_plot(rsut_NorESM2_plot)
save_plot(rsutcs_NorESM2_plot)

## Arrange plots in grid
clt_plots <- grid.arrange(clt_CAM5ATRAS_plot,clt_CESM1_plot,clt_E3SM_plot,#clt_GEOS_plot,
                          clt_GISS_plot,clt_NorESM2_plot)

if (scenario == "shp-60p-red" || scenario == "shp-atl-shift" || scenario == "shp-ind-shift") {
  if (scenario == "shp-atl-shift") {
    loadso2_plots <- grid.arrange(loadso2_CAM5ATRAS_plot,loadso2_CESM1_plot,
                                  #loadso2_GEOS_plot,
                                  loadso2_GFDL_plot,loadso2_GISS_plot,
                                  loadso2_NorESM2_plot)
  } else {
    loadso2_plots <- grid.arrange(loadso2_CAM5ATRAS_plot,loadso2_CESM1_plot,
                                  loadso2_E3SM_plot,#loadso2_GEOS_plot,
                                  loadso2_GFDL_plot,loadso2_GISS_plot,
                                  loadso2_NorESM2_plot)
  }
} else {
  loadso2_plots <- grid.arrange(loadso2_CAM5ATRAS_plot,loadso2_CESM1_plot,
                                loadso2_E3SM_plot,#loadso2_GEOS_plot,
                                loadso2_GISS_plot,loadso2_NorESM2_plot)
}

if (scenario != "shp-ind-shift-1950" & scenario != "shp-30p-red") {
  if (scenario == "shp-atl-shift") {
    loadso4_plots <- grid.arrange(loadso4_CAM5ATRAS_plot,loadso4_CESM1_plot,
                                  #loadso4_GEOS_plot,
                                  loadso4_GFDL_plot,loadso4_GISS_plot,
                                  loadso4_NorESM2_plot)
  } else {
    loadso4_plots <- grid.arrange(loadso4_CAM5ATRAS_plot,loadso4_CESM1_plot,
                                  loadso4_E3SM_plot,#loadso4_GEOS_plot,
                                  loadso4_GFDL_plot,loadso4_GISS_plot,
                                  loadso4_NorESM2_plot)
  }
} else {
  loadso4_plots <- grid.arrange(loadso4_CAM5ATRAS_plot,loadso4_CESM1_plot,
                                loadso4_E3SM_plot,#loadso4_GEOS_plot,
                                loadso4_GISS_plot,loadso4_NorESM2_plot)
}

rlut_plots <- grid.arrange(rlut_CAM5ATRAS_plot,rlut_CESM1_plot,rlut_E3SM_plot,#rlut_GEOS_plot,
                           rlut_GFDL_plot,rlut_GISS_plot,
                           rlut_NorESM2_plot)

rlutcs_plots <- grid.arrange(rlutcs_CAM5ATRAS_plot,rlutcs_CESM1_plot,
                             rlutcs_E3SM_plot,rlutcs_GFDL_plot,#rlutcs_GEOS_plot,
                             rlutcs_GISS_plot,rlutcs_NorESM2_plot)

rsut_plots <- grid.arrange(rsut_CAM5ATRAS_plot,rsut_CESM1_plot,rsut_E3SM_plot,#rsut_GEOS_plot,
                           rsut_GFDL_plot,rsut_GISS_plot,
                           rsut_NorESM2_plot)

rsutcs_plots <- grid.arrange(rsutcs_CAM5ATRAS_plot,rsutcs_CESM1_plot,
                             rsutcs_E3SM_plot,rsutcs_GFDL_plot,#rsutcs_GEOS_plot,
                             rsutcs_GISS_plot,rsutcs_NorESM2_plot)

# Save plots as pdf
# setwd(paste0(file_dir,"/output/scenario/",scenario))
# pdf(paste0(scenario,'_maps_per-diff.pdf'), height = 11, width = 8.5, paper = "letter")
# 
# grid.draw(clt_plots)
# grid.newpage()
# grid.draw(loadso2_plots)
# grid.newpage()
# grid.draw(loadso4_plots)
# grid.newpage()
# grid.draw(rlut_plots)
# grid.newpage()
# grid.draw(rsut_plots)
# grid.newpage()
# grid.draw(rlutcs_plots)
# grid.newpage()
# grid.draw(rsutcs_plots)
# dev.off()



# get average rsut plot
rsut_all <- (rsut_CAM5ATRAS+rsut_CESM1+rsut_E3SM+rsut_GEOS+rsut_GISS+rsut_GFDL+rsut_NorESM2)/7
rsut_all.at <- drop_na(rsut_all)
bound_rsut_all <- max(abs(rsut_all.at$layer))
rsut_all.at <- seq(-bound_rsut_all,bound_rsut_all,length.out=25)

rsut_all_plot <- get_plot("rsut","all",rsut_all,rsut_all.at,"Upwelling Shortwave Radiation \n at TOA (%) - All Models")

# get average rlut plot
rlut_all <- (rlut_CAM5ATRAS+rlut_CESM1+rlut_E3SM+rlut_GEOS+rlut_GISS+rlut_GFDL+rlut_NorESM2)/7
rlut_all.at <- drop_na(rlut_all)
bound_rlut_all <- max(abs(rlut_all.at$layer))
rlut_all.at <- seq(-bound_rlut_all,bound_rlut_all,length.out=25)

rlut_all_plot <- get_plot("rlut","all",rlut_all,rsut_all.at,"Upwelling Longwave Radiation \n at TOA (%) - All Models")
