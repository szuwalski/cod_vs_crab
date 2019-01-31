## Running VAST multi-species example
# install_github("nwfsc-assess/geostatistical_delta-GLMM") 

## from VAST package vignette
library(TMB)
library(VAST)

#### Settings ####
Version = "VAST_v4_4_0"

#spatial settings
Method <- "Mesh"
grid_size_km <- 50
n_x <- 100 # number of stations
Kmeans_Config <- list("randomseed"=1, "nstart"=100,"iter.max"=1e3)

# whether to include spatial and s-t variation, whether it is autocorrelated, and whether there is overdispersion
FieldConfig = c(Omega1 = 3, Epsilon1 = 3, Omega2 = 3,
                Epsilon2 = 3)
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Vessel = 0, VesselYear = 0)
ObsModel = c(2, 0)

#outputs we want
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

strata.limits <- data.frame(STRATA = "All_areas")

# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Atheresthes stomias","Gadus chalcogrammus","Hippoglossoides elassodon")

# save settings and data in a list
DateFile = paste0(getwd(),'/vast/output/')
dir.create(DateFile)
Record = list(Version = Version, Method = Method, grid_size_km = grid_size_km,
              n_x = n_x, FieldConfig = FieldConfig, RhoConfig = RhoConfig,
              OverdispersionConfig = OverdispersionConfig, ObsModel = ObsModel,
              Kmeans_Config = Kmeans_Config, Region = Region,
              Species_set = Species_set, strata.limits = strata.limits)
save(Record, file = file.path(DateFile, "Record.RData"))
capture.output(Record, file = paste0(DateFile, "Record.txt"))

#### Prepare Data ####
# Depending upon the Data_Set chosen, we load archived data sets that are distributed with the package.
# Each archived data set is then reformatted to create a data-frame Data_Geostat with a standardized set of
# columns. For a new data set, the user is responsible for formatting Data_Geostat appropriately to match
# this format. We show the first six rows of Data_Geostat given that Data_Set = Data_Set.

DF = FishData::download_catch_rates(survey = "Eastern_Bering_Sea",species_set = Species_set)
Data_Geostat = data.frame(spp = DF[, "Sci"], Year = DF[,"Year"], 
                          Catch_KG = DF[, "Wt"], 
                          AreaSwept_km2 = 0.01,
                          Vessel = 0, 
                          Lat = DF[, "Lat"], 
                          Lon = DF[, "Long"])

# make the extrapolatoin grid, builiding an object used to determine areas to extrapoolate densities to when calaculating indices
Extrapolation_List = make_extrapolation_info(Region = Region,strata.limits = strata.limits)

# generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`
Spatial_List = make_spatial_info( grid_size_km=grid_size_km, 
                                  n_x=n_x, 
                                  Method=Method, 
                                  Lon=Data_Geostat[,'Lon'], 
                                  Lat=Data_Geostat[,'Lat'], 
                                  Extrapolation_List=Extrapolation_List, 
                                  DirPath=DateFile, Save_Results=FALSE )

# Add the knots to the the data
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i)

#### Build and Run Model ####

# in order to estimate params, build a list of data-inputs used for param estimation
TmbData = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
                  "s_i"=Data_Geostat[,'knot_i']-1, 
                  "t_i"=Data_Geostat[,'Year'], 
                  "a_xl"=Spatial_List$a_xl, 
                  "MeshList"=Spatial_List$MeshList, 
                  "GridList"=Spatial_List$GridList, 
                  "Method"=Spatial_List$Method, 
                  "Options"=Options )

## Build TMB model object
TmbList = Build_TMB_Fn("TmbData"=TmbData, 
                       "RunDir"=DateFile, 
                       "Version"=Version, 
                       "RhoConfig"=RhoConfig, 
                       "loc_x"=Spatial_List$loc_x, 
                       "Method"=Method)
Obj = TmbList[["Obj"]]

# Use gradient-based nonlinear minimizer to identify maximum likelihood esimates for fixed effects
Opt = TMBhelper::Optimize( obj=Obj, 
                           lower=TmbList[["Lower"]], 
                           upper=TmbList[["Upper"]], 
                           getsd=TRUE, 
                           savedir=DateFile, 
                           bias.correct=TRUE, 
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl"), newtonsteps=1 )

## bundle and save output!
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

#### Visualize and Diagnose Model Outputs####
load(paste0(DateFile,"Save.RData"))
for(i in 1:length(Save)) assign(names(Save)[i], Save[[i]])
# It is always good practice to conduct exploratory analysis of data.  Here, I visualize the spatial distribution of data.  
# Spatio-temporal models involve the assumption that the probability of sampling a given location is statistically independent 
# of the probability distribution for the response at that location.  
# So if sampling "follows" changes in density, then the model is probably not appropriate!
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )

# Here I print the diagnostics generated during parameter estimation, and I confirm that 
# (1) no parameter is hitting an upper or lower bound and 
# (2) the final gradient for each fixed-effect is close to zero. 
# For explanation of parameters, please see `?Data_Fn`.

Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## Encounter probability
# we check whether observed encounter frequencies for either low or high probability samples 
# are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)

## Fit to residuals of catch-rates (Q-Q plot)
Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP=NULL,
                              FileName_Phist=NULL, 
                              FileName_QQ="Q-Q_plot", FileName_Qhist=NULL, DateFile=DateFile) 

## finally we can visualize residuals on a map
# Define years to plot and labels for each plotted year
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

plot_residuals(Lat_i=Data_Geostat[,'Lat'], 
               Lon_i=Data_Geostat[,'Lon'], 
               TmbData=TmbData, 
               Report=Report, 
               Q=Q, 
               savedir=DateFile, 
               MappingDetails=MapDetails_List[["MappingDetails"]], 
               PlotDF=MapDetails_List[["PlotDF"]],
               MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
               Xlim=MapDetails_List[["Xlim"]], 
               Ylim=MapDetails_List[["Ylim"]], 
               FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], 
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)
