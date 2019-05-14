# VAST with updated model version (streamlined)
# Constant intercepts, ordination version
# No density covariates
library(TMB)
library(VAST)
library(tidyverse)
map <- maps::map

# Load data
source(paste0(getwd(),"/vast/make_combined_data.R"))
fp = paste0(getwd(),'/vast/output/ordination/')
dir.create(fp)

model_directory <- paste0(getwd(),'/vast/output/ordination/model/')
dir.create(model_directory)

example <- load_example(data_set="five_species_ordination")

# Make settings
settings <- make_settings(n_x = 100,Region = "eastern_bering_sea",
                          purpose = "ordination",fine_scale = TRUE,
                          strata.limits = data.frame(STRATA = "All_areas"),
                          zone=NA,n_categories = 2,treat_nonencounter_as_zero = TRUE)

# Calculate density covariates

# X_gtp <- make_density_covariates(dat = dat_combined,Region = settings$Region,strata.limits = settings$strata.limits,
#                                  grid_size_km = settings$grid_size_km,n_x = settings$n_x,Method = settings$Method,Lon_i = dat_combined[,'lon'],
#                                  Lat_i=dat_combined[,'lat'],fine_scale=settings$fine_scale)


# Customize some settings
settings
settings$RhoConfig[c('Beta1','Beta2')]=3 # make Betas constant among years; Epsilons follow an AR1 structure
settings$ObsModel=c(2,1)
settings$bias.correct # no bias correction for now

# Run model
fit <- fit_model(settings = settings,
                 Lat_i = dat_combined[,'lat'],Lon_i = dat_combined[,'lon'],
                 t_iz = dat_combined[,'year'],c_iz = as.numeric(dat_combined[,'spp'])-1,
                 b_i = dat_combined[,'abun'],a_i = dat_combined[,'area_km2'],X_gtp = NULL,
                 # X_itp = array(0,dim=c(length(dat_combined[,'abun']),dim(X_gtp)[2:3])),
                 working_dir = model_directory,silent = FALSE)
results = plot_results(fit = fit,settings = settings,working_dir = fp)