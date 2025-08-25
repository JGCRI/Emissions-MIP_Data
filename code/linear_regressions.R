# Linear regression plots -- slopes, intercepts, R2 values

# Load in data from summary_plots_diff.R
#region <- "NH-pacific"
data <- loadso2_loadso4_combined
data_name <- "loadso2_loadso4_combined"

# Split data by model
data_CAMATRAS <- data[data$model=="CAM-ATRAS",]
data_CESM1 <- data[data$model=="CESM1",]
data_E3SM <- data[data$model=="E3SM",]
data_GFDL <- data[data$model=="GFDL-ESM4",]
data_GISS <- data[data$model=="GISS-E2.1",]
data_NorESM <- data[data$model=="NorESM2",]

# Linear regressions
model_CAMATRAS <- lm(value.y ~ value.x, data = data_CAMATRAS)
r_CAMATRAS <- summary(model_CAMATRAS)$r.squared
intercept_CAMATRAS <- model_CAMATRAS$coefficients[1]
slope_CAMATRAS <- model_CAMATRAS$coefficients[2]

model_CESM1 <- lm(value.y ~ value.x, data = data_CESM1)
r_CESM1 <- summary(model_CESM1)$r.squared
intercept_CESM1 <- model_CESM1$coefficients[1]
slope_CESM1 <- model_CESM1$coefficients[2]

model_E3SM <- lm(value.y ~ value.x, data = data_E3SM)
r_E3SM <- summary(model_E3SM)$r.squared
intercept_E3SM <- model_E3SM$coefficients[1]
slope_E3SM <- model_E3SM$coefficients[2]

model_GFDL <- lm(value.y ~ value.x, data = data_GFDL)
r_GFDL <- summary(model_GFDL)$r.squared
intercept_GFDL <- model_GFDL$coefficients[1]
slope_GFDL <- model_GFDL$coefficients[2]

model_GISS <- lm(value.y ~ value.x, data = data_GISS)
r_GISS <- summary(model_GISS)$r.squared
intercept_GISS <- model_GISS$coefficients[1]
slope_GISS <- model_GISS$coefficients[2]

model_NorESM <- lm(value.y ~ value.x, data = data_NorESM)
r_NorESM <- summary(model_NorESM)$r.squared
intercept_NorESM <- model_NorESM$coefficients[1]
slope_NorESM <- model_NorESM$coefficients[2]

# Put values in table
models <- c("CAM-ATRAS","CESM1","E3SM","GFDL-ESM4","GISS-E2.1","NorESM2")
slopes <- c(slope_CAMATRAS,slope_CESM1,slope_E3SM,slope_GFDL,slope_GISS,slope_NorESM) %>% as.numeric()
intercepts <- c(intercept_CAMATRAS,intercept_CESM1,intercept_E3SM,intercept_GFDL,intercept_GISS,intercept_NorESM) %>% as.numeric()
r_squared <- c(r_CAMATRAS,r_CESM1,r_E3SM,r_GFDL,r_GISS,r_NorESM)

table <- data.frame(models=models,
                    slopes=slopes,
                    intercepts=intercepts,
                    r_squared=r_squared)

write.csv(table,file=paste0(emi_dir,"/output/",region,"/",region,"_linregs_",data_name,".csv"))
