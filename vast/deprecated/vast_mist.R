# VAST with updated model version (streamlined)
# Interactions (MIST) version
# No density covariates
library(TMB)
library(VAST)
library(tidyverse)
map <- maps::map

# Load data
source(paste0(getwd(),"/vast/make_combined_data.R"))
fp = paste0(getwd(),'/vast/output/mist/')
dir.create(fp)

model_directory <- paste0(getwd(),'/vast/output/mist/model/')
dir.create(model_directory)

# Make settings
settings <- make_settings(n_x = 25,Region = "eastern_bering_sea",
                          purpose = "MICE",fine_scale = FALSE,
                          strata.limits = data.frame(STRATA = "All_areas"),
                          n_categories = length(unique(dat_combined$spp)))

# Calculate density covariates

# X_gtp <- make_density_covariates(dat = dat_combined,Region = settings$Region,strata.limits = settings$strata.limits,
#                                  grid_size_km = settings$grid_size_km,n_x = settings$n_x,Method = settings$Method,Lon_i = dat_combined[,'lon'],
#                                  Lat_i=dat_combined[,'lat'],fine_scale=settings$fine_scale)


# Customize some settings
settings
settings$VamConfig
settings$VamConfig[3] <- 0
settings$RhoConfig <- c(3,3,4,0)
settings$FieldConfig
settings$Options[c(5,6)] <- FALSE
settings$bias.correct=FALSE # no bias correction for now

# Run model
fit <- fit_model(settings = settings,
                 Lat_i = dat_combined[,'lat'],Lon_i = dat_combined[,'lon'],
                 t_iz = dat_combined[,'year'],c_iz = as.numeric(dat_combined[,'spp'])-1,
                 b_i = dat_combined[,'abun'],a_i = dat_combined[,'area_km2'],
                 # X_itp = array(0,dim=c(length(dat_combined[,'abun']),dim(X_gtp)[2:3])),
                 working_dir = model_directory,silent = FALSE)
results = plot_results(fit = fit,settings = settings,working_dir = fp)