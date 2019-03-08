### VAST relating immature snow crab and Pacific cod###

#### Settings ####
library(TMB)
library(VAST)
library(tidyverse)
map <- maps::map

Version = get_latest_version( package="VAST" )
# VAST v7_0_0
# Version = "VAST_v4_4_0"

#### Import Data ####
load("data/longform_opilio2.Rdata")
load("data/cod_dat_clean.Rdata")

fp = paste0(getwd(),'/vast/output/opilio_mf/')
dir.create(fp)

# filter to just immature crabs, summarise total by station/year, and row-bind cod data
dat_opilio <- opi_dat_long %>% 
  filter(maturity=="Mature",sex=="Female",units=="kilos") %>% 
  group_by(towid,year,lon,lat,area_km2) %>% 
  summarise(biomass=sum(value,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(spp="Chionoecetes opilio",vessel=0)

# combine and add missing zeros for missing tows
dat <- dat_opilio %>% mutate(spp=factor(spp)) %>% 
  filter(!is.na(lat),!is.na(lon),year<2018) %>% 
  # do we fill in zeroes? for all tow/spp combinations that were previously missing
  complete(spp,nesting(year,lat,lon),fill=list(biomass=0,area_km2=0.01,vessel=0)) %>%
  ungroup() %>%
  # complete(spp,nesting(year,lat,lon),fill=list(area_km2=0.01,vessel=0))
  select(spp,year,lon,lat,area_km2,biomass,vessel)

# data check- number of zeroes and NAs
zeroes_check <- dat %>% 
  group_by(year,spp) %>% 
  summarise(num_zeroes=sum(biomass==0),num_NA=sum(is.na(biomass))) %>% 
  arrange(spp,year)


# make the extrapolatoin grid, building an object used to determine areas to extrapolate densities to when calculating indices
# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Chionoecetes opilio")

strata.limits <- data.frame(STRATA = "All_areas")

Extrapolation_List = make_extrapolation_info(Region = Region,strata.limits = strata.limits)

# generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`

# Stochastic partial differential equation (SPDE) with geometric anisotropy
Method <- "Mesh"
Aniso <- 1
grid_size_km=25
# number of 'knots'
n_x <- 100

Spatial_List = make_spatial_info(grid_size_km=grid_size_km,
                                 n_x=n_x, 
                                 Method=Method, 
                                 Lon_i=dat$lon, 
                                 Lat_i=dat$lat, 
                                 Extrapolation_List=Extrapolation_List, 
                                 randomseed=1,nstart=100,iter.max=1e3,
                                 DirPath=fp, Save_Results=FALSE)
save(Spatial_List,file= paste0(fp,"Spatial_List.Rdata"))

# Add the knots to the the data
dat <- dat %>% mutate(knot_i=Spatial_List$knot_i)

# Spatial settings

# whether to include spatial (omega) and s-t (epsilon) variation
# for two linear predictors- 1. encounter probability, and 2. positive catch rate
# for a univariate model, these are 0 ("turned off") or 1
# for a multivariate model, can be any whole number 0:C, where C is the number of categories
# indicating the number of factors to estimate
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1,
                Epsilon2 = 1)
# is there temporal correlation in the intercepts (Beta) or s-t variation (Epsilon)?
# 0: each year as fixed effect; 
# 1: each year as random following IID distribution; 
# 2: each year as random following a random walk; 
# 3: constant among years as fixed effect; 
# 4: each year as random following AR1 process
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)

# is there overdispersion? often attributable to 'vessel effects'
# 0 means no overdispersion
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

# Options to estimate interactions

# Options to estimate interactions, where 
# first slot selects method for forming interaction matrix
# second indicates rank
# third indicates whether to incorporate effect of F_ct, 
# fourth indicates whether to add a new "t=0" year (while incrementing all t_i inputs) which represents B0
# Method = 2 means Real-eigenvalues
VamConfig	<- c("Method"=0,"Rank"=0, "Timing" =0)

# what distributions and link functions to use?
# first is the distribution for positive catch rates, 
# second is the functional form for encounter probabilities,
# we choose a lognormal distribution for positive catch rates and a Poisson delta-model using log-link for enc prob and log-link for catch rates
ObsModel = c(1, 1)

#outputs we want
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

# save settings and data in a list

Record = list(Version = Version, 
              Method = Method,
              n_x = n_x, 
              FieldConfig = FieldConfig, 
              RhoConfig = RhoConfig,
              OverdispersionConfig = OverdispersionConfig, 
              VamConfig=VamConfig, 
              ObsModel = ObsModel,
              Region = Region,
              Species_set = Species_set, 
              strata.limits = strata.limits)
save(Record, file = file.path(fp, "Record.RData"))
capture.output(Record, file = paste0(fp, "Record.txt"))


#### Build and Run Model ####

# in order to estimate params, build a list of data-inputs used for param estimation
TmbData = make_data("Version"=Version, 
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "VamConfig" = VamConfig,
                  "ObsModel"=ObsModel, 
                  "c_iz"=as.numeric(dat$spp)-1,
                  "b_i"=dat$biomass, 
                  "a_i"=dat$area_km2, 
                  "v_i"=as.numeric(dat$vessel)-1, 
                  "s_i"=dat$knot_i-1, 
                  "t_iz"=dat$year,
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Aniso"=Aniso,
                  "Options"=Options)

## Build TMB model object
TmbList = make_model("build_model"=TRUE,
                       "TmbData"=TmbData, 
                       "RunDir"=fp, 
                       "Version"=Version, 
                       "RhoConfig"=RhoConfig, 
                       "loc_x"=Spatial_List$loc_x, 
                       "Method"=Method)
Map = TmbList$Map
Parameters = TmbList$Parameters

# Starting values for eigenvalues of B_cc
if( VamConfig['Method']==2 ){
  if( RhoConfig[3]==4 ){
    Parameters$Epsilon_rho1_f[] = 0.5
    if( VamConfig['Rank']>=1 ) Parameters[["Psi_fr"]][cbind(1:VamConfig['Rank'],1:VamConfig['Rank'])] = 0.1*(1:VamConfig['Rank']) - mean(Parameters$Epsilon_rho1_f[])
  }
}

# Optimize
Obj = TmbList[["Obj"]]
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# # Bounds
# if( VamConfig['Method']==2 ){
#   if( VamConfig['Rank']>0 ){
#     TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Psi_fr")[length(which(names(TmbList[["Lower"]])=="Psi_fr"))+1-1:VamConfig['Rank']]] = 0
#     TmbList[["Upper"]][which(names(TmbList[["Lower"]])=="Psi_fr")[length(which(names(TmbList[["Lower"]])=="Psi_fr"))+1-1:VamConfig['Rank']]] = 1
#   }
#   if( RhoConfig['Epsilon1'] %in% c(4,5) ){
#     TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Epsilon_rho1")] = 0
#     TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Epsilon_rho1_f")] = 0
#   }
# }

# Optimize                                         #
# Bias Correct?
BiasCorr= TRUE
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], savedir=fp, getsd=TRUE, 
                           bias.correct=BiasCorr, newtonsteps=1, 
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=c("Index_cyl")) )
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=dat, "Map"=Map)
save(Save, file=paste0(fp,"Save.RData"))
if("opt" %in% names(Opt)) capture.output( Opt$opt, file=paste0(fp,"parameter_estimates.txt"))
# Use gradient-based nonlinear minimizer to identify maximum likelihood esimates for fixed effects

#### Model Output Diagnostics ####

names(dat) <- c("spp","Year","Lon","Lat","area_km2", "Catch_KG","vessel","knot_i")
# diagnostic plots
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=dat, PlotDir=fp )

# convergence
pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] ) 

## Check encounter probabilities
# Check whether observed encounter frequencies for either low or high probability samples 
# are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=dat, DirName=fp)

# Check positive catch-rate
#nWe can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
# A good Q-Q plot will have residuals along the one-to-one line.  

Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=fp )

## Plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(dat$Year),max(dat$Year))
Years2Include = which( Year_Set %in% sort(unique(dat$Year)))

# Plot Pearson residuals.  
# If there are visible patterns (areas with consistently positive or negative residuals accross or within years) 
# then this is an indication of the model "overshrinking" results towards the intercept, and model results should then be treated with caution.  

plot_residuals(Lat_i=dat[,'Lat'], Lon_i=dat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, 
               savedir=fp, MappingDetails=MapDetails_List[["MappingDetails"]], 
               PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
               Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=fp, Year_Set=Year_Set, 
               Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], 
               Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

## Direction of "geometric anisotropy"

# We can visualize which direction has faster or slower decorrelation (termed "geometric anisotropy")
plot_anisotropy( FileName=paste0(fp,"Aniso.png"), Report=Report, TmbData=TmbData )


#### View Model Outputs ####

## Density surface for each year
# predicted density, but other options are obtained via other integers passed to `plot_set` as described in `?plot_maps`
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, 
                    Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
                    Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=fp, Year_Set=Year_Set, 
                    Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], 
                    Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0),
                    cex=1.8, plot_legend_fig=FALSE)

# We can also extract density predictions at different locations, for use or plotting in other software. 
# This is output in UTM using zone `r Extrapolation_List$zone-ifelse(Extrapolation_List$flip_around_dateline,30,0)`
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], 
                 "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], 
                 "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## Index of abundance
# The index of abundance is generally most useful for stock assessment models.
Index = plot_biomass_index( DirName=fp, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] ) 

# Range expansion/contraction
plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=fp, Year_Set=Year_Set)
