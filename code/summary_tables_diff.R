# ------------------------------------------------------------------------------
# Program Name: summary_tables_diff.R
# Authors: Ciara Donegan
# Date Last Modified: January 10, 2024
# Program Purpose: Produces summary tables of the difference between the
# perturbations and the reference case averaged over all years
# Input Files: ~Emissions-MIP/input/
# Output Files: ~Emissions-MIP/output/
# TODO:
# ------------------------------------------------------------------------------

# just run summary_plots_diff for now to get all needed info

# clear model data
rm(model_CAM_ATRAS,model_CESM,model_E3SM,model_GEOS,
   model_GFDL,model_GISS,model_NorESM)

# lm slopes
data <- loadso2_rad_combined
data_name <- deparse(substitute(loadso2_rad_combined))
model_CAM_ATRAS <- lm(value.y ~ value.x, data=filter(data,model=="CAM-ATRAS"))
model_CESM <- lm(value.y ~ value.x, data=filter(data,model=="CESM1"))
model_E3SM <- lm(value.y ~ value.x, data=filter(data,model=="E3SM"))
model_GEOS <- lm(value.y ~ value.x, data=filter(data,model=="GEOS"))
model_GFDL <- lm(value.y ~ value.x, data=filter(data,model=="GFDL-ESM4"))
model_GISS <- lm(value.y ~ value.x, data=filter(data,model=="GISS-E2.1"))
model_NorESM <- lm(value.y ~ value.x, data=filter(data,model=="NorESM2"))

get_slope <- function(model) {
  coef(model)[2]
}

get_r2 <- function(model) {
  summary(model)$r.squared
}

get_intercept <- function(model) {
  coef(model)[1]
}

# all models
df <- data.frame(models=c("CAM-ATRAS","CESM1","E3SM","GEOS","GFDL-ESM4",
                          "GISS-E2.1","NorESM2"),
                 resolution=c("129, 154, 180, 204m","124, 149, 173, 197m",
                              "25, 54, 72, 77, 82, 87m","58, 131, 65, 133m",
                              "35, 50, 75, 90, 120m","170, 190, 220, 240m",
                              "127, 152, 176, 201m"),
                 bottom.thickness=c(129,124,25,58,35,170,127),
                 layers=c(4,4,6,4,5,4,4),
                 slopes=c(get_slope(model_CAM_ATRAS),get_slope(model_CESM),
                          get_slope(model_E3SM),get_slope(model_GEOS),
                          get_slope(model_GFDL),get_slope(model_GISS),
                          get_slope(model_NorESM)),
                 intercept=c(get_intercept(model_CAM_ATRAS),get_intercept(model_CESM),
                             get_intercept(model_E3SM),get_intercept(model_GEOS),
                             get_intercept(model_GFDL),get_intercept(model_GISS),
                             get_intercept(model_NorESM)),
                 r.squared=c(get_r2(model_CAM_ATRAS),get_r2(model_CESM),
                             get_r2(model_E3SM),get_r2(model_GEOS),
                             get_r2(model_GFDL),get_r2(model_GISS),
                             get_r2(model_NorESM)))

# missing GFDL
df <- data.frame(models=c("CAM-ATRAS","CESM1","E3SM","GEOS",
                          "GISS modelE","NorESM2"),
                 resolution=c("129, 154, 180, 204m","124, 149, 173, 197m",
                              "25, 54, 72, 77, 82, 87m","58, 131, 65, 133m",
                              "170, 190, 220, 240m",
                              "127, 152, 176, 201m"),
                 bottom.thickness=c(129,124,25,58,170,127),
                 layers=c(4,4,6,4,4,4),
                 slopes=c(get_slope(model_CAM_ATRAS),get_slope(model_CESM),
                          get_slope(model_E3SM),get_slope(model_GEOS),
                          get_slope(model_GISS),
                          get_slope(model_NorESM)),
                 intercept=c(get_intercept(model_CAM_ATRAS),get_intercept(model_CESM),
                             get_intercept(model_E3SM),get_intercept(model_GEOS),
                             get_intercept(model_GISS),
                             get_intercept(model_NorESM)),
                 r.squared=c(get_r2(model_CAM_ATRAS),get_r2(model_CESM),
                             get_r2(model_E3SM),get_r2(model_GEOS),
                             get_r2(model_GISS),
                             get_r2(model_NorESM)))

write.csv(df,paste0("C:\\Users\\done231\\OneDrive - PNNL\\Desktop\\Phase1b_data\\output\\",region,"\\",region,"_table_",data_name,".csv"))
