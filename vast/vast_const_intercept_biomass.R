### VAST relating immature and spawner crabs and size-structured Pacific cod###
### We do this based on numbers (not biomass) ###
## with temperature and depth as covariates ##

#### Settings ####
library(TMB)
library(VAST)
library(tidyverse)
map <- maps::map

fp = paste0(getwd(),'/vast/output/const_intercept_bio/')
dir.create(fp)

model_directory <- paste0(getwd(),'/vast/output/const_intercept_bio/model/')
dir.create(model_directory)

####
Version = get_latest_version( package="VAST" )

#### Import Data ####
source(here('cod_crab_biomass_combined.R'))

# make the extrapolatoin grid, building an object used to determine areas to extrapolate densities to when calculating indices
# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Opilio immature","Opilio spawner","Pacific cod")

dat <- dat_combined %>% filter(spp %in% Species_set)

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
                                 DirPath=fp,Save_Results=FALSE)
save(Spatial_List,file= paste0(fp,"Spatial_List.Rdata"))

# Add the knots to the the data
dat <- dat %>% mutate(knot_i=Spatial_List$knot_i)

# Spatial settings

# whether to include spatial (omega) and s-t (epsilon) variation
# for two linear predictors- 1. encounter probability, and 2. positive catch rate
# for a univariate model, these are 0 ("turned off") or 1
# for a multivariate model, can be any whole number 0:C, where C is the number of categories
# indicating the number of factors to estimate
FieldConfig = c(Omega1 = 3, Epsilon1 = 3, Omega2 = 3,
                Epsilon2 = 3)
# is there temporal correlation in the intercepts (Beta) or s-t variation (Epsilon)?
# 0: each year as fixed effect; 
# 1: each year as random following IID distribution; 
# 2: each year as random following a random walk; 
# 3: constant among years as fixed effect; 
# 4: each year as random following AR1 process
RhoConfig = c(Beta1 = 3, Beta2 = 3, Epsilon1 = 0, Epsilon2 = 0)

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
# we choose a gamma distribution for positive catch rates and a Poisson-link delta-model for enc prob
# 
ObsModel = c(2, 1)

#outputs we want
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

# save settings and data in a list

Record = list(Version = Version, 
              Method = Method,
              n_x = n_x, 
              Aniso = Aniso,
              FieldConfig = FieldConfig, 
              RhoConfig = RhoConfig,
              OverdispersionConfig = OverdispersionConfig, 
              VamConfig=VamConfig, 
              ObsModel = ObsModel,
              Region = Region,
              Species_set = Species_set, 
              strata.limits = strata.limits,
              Options=Options)
save(Record, file = file.path(fp, "Record.RData"))
capture.output(Record, file = paste0(fp, "Record.txt"))

# Finally, we need the density covariates (temperature and depth)
# this will go into VAST::make_data as a 3 dimensional array
# with dimensions nknots x nyears x nvars
# We use temperature and depth data interpolated in another script
library(RANN)

X_gtp <- make_density_covariates(type="loc_x")


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
                    "X_xtp"=X_gtp,
                    "MeshList"=Spatial_List$MeshList,
                    "GridList"=Spatial_List$GridList,
                    "Method"=Spatial_List$Method,
                    "Aniso"=Aniso,
                    "Options"=Options)

## Build TMB model object
TmbList = make_model("build_model"=TRUE,
                     "TmbData"=TmbData, 
                     "RunDir"=model_directory, 
                     "Version"=Version, 
                     "RhoConfig"=RhoConfig, 
                     "loc_x"=Spatial_List$loc_x, 
                     "Method"=Method)
Map = TmbList$Map
Parameters = TmbList$Parameters

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
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=c("Index_cyl")))
Opt$time_for_run
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=dat, "Map"=Map)
save(Save, file=paste0(fp,"Save.RData"))
if("opt" %in% names(Opt)) capture.output( Opt$opt, file=paste0(fp,"parameter_estimates.txt"))
# Use gradient-based nonlinear minimizer to identify maximum likelihood esimates for fixed effects

#### Model Output Diagnostics ####
# if not re-running, load
load(paste0(fp,"Record.Rdata"))
load(paste0(fp,"Save.RData"))
load(paste0(fp,"Spatial_List.Rdata"))
load(paste0(fp,"parameter_estimates.Rdata"))
dat <- Save$Data
list2env(Record,envir = environment())

names(dat) <- c("spp","Year","Lon","Lat","area_km2", "Catch_KG","vessel","knot_i")
source("C:/Users/owenr/Documents/github/cod_vs_crab/vast/output__helper_fxns.R")
# diagnostic plots
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=dat, PlotDir=fp )

# convergence
pander::pandoc.table(Save$Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] ) 

## Check encounter probabilities
# Check whether observed encounter frequencies for either low or high probability samples 
# are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic( Report=Save$Report, Data_Geostat=dat, DirName=fp)

# Check positive catch-rate
#nWe can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
# A good Q-Q plot will have residuals along the one-to-one line.  

Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=fp )

save(Q,file=paste0(fp,"quantile_diag.Rdata"))

## Plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(dat$Year),max(dat$Year))
Years2Include = which( Year_Set %in% sort(unique(dat$Year)))

# Plot Pearson residuals.  
# If there are visible patterns (areas with consistently positive or negative residuals accross or within years) 
# then this is an indication of the model "overshrinking" results towards the intercept, and model results should then be treated with caution.  
TmbData$n_x <- TmbData$n_g
res <- plot_residuals(Lat_i=dat[,'Lat'], Lon_i=dat[,'Lon'], TmbData=TmbData, Report=Save$Report, Q=Q, 
                      savedir=fp, MappingDetails=MapDetails_List[["MappingDetails"]], 
                      PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
                      Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=fp, Year_Set=Year_Set, 
                      Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], 
                      Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)
save(res,file=paste0(fp,"residuals.Rdata"))

## Direction of "geometric anisotropy"

# We can visualize which direction has faster or slower decorrelation (termed "geometric anisotropy")
plot_anisotropy( FileName=paste0(fp,"Aniso.png"), Report=Report, TmbData=TmbData )


#### View Model Outputs ####

## Plot spatial and spatio-temporal covariance

# Spatial and spatio-temporal covariance among species in encounter probability and positive catch rates 
# (depending upon what is turned on via `FieldConfig`)
Cov_List = Summarize_Covariance( Report=Report, ParHat=Obj$env$parList(), Data=TmbData, 
                                 SD=Save$Opt$SD, plot_cor=FALSE, category_names=levels(dat$spp), 
                                 plotdir=fp, mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

## Density surface for each year
# predicted density, but other options are obtained via other integers passed to `plot_set` as described in `?plot_maps`
plot_density(Region,Spatial_List,Report,Extrapolation_List,dat,fp,saveplots = T)
# Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, 
#                     Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
#                     Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=fp, Year_Set=Year_Set, 
#                     Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], category_names = levels(dat$spp),
#                     Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0),
#                     cex=1.8, plot_legend_fig=FALSE)

## Index of abundance
# The index of abundance is generally most useful for stock assessment models.
Index = plot_biomass_index( DirName=fp, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, 
                            category_names=levels(dat$spp),Years2Include=Years2Include, use_biascorr=TRUE )
pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] ) 

# Range expansion/contraction
plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=fp, 
                 category_names = levels(dat$spp),Year_Set=Year_Set)


## Plot factors
# Finally, we can inspect the factor-decomposition for community-level patterns.  
# This generates many plots
factors <- Plot_factors( Report=Save$Report, ParHat=Save$ParHat, Data=TmbData, SD=Opt$SD, mapdetails_list=MapDetails_List, 
                         Year_Set=Year_Set, category_names=levels(dat$spp), plotdir=fp )
fct_loadings <- plot_fct_loadings(factors,dat,rotated = TRUE,fp)

# plot 2-dimensional loadings
library(ggsci)
pca <- fct_loadings %>% 
  filter(fct_num < 3) %>% 
  select(spp,fct,fct_num,load) %>%
  unite("fct",fct,fct_num) %>% 
  group_by(spp) %>% 
  spread(fct,load) %>% 
  mutate(spp_sym= ifelse(spp %in% c('Opilio Immature','Opilio Spawner'),'crab','cod'))
pca_o1_plot <- pca %>% 
  ggplot(aes(Omega1_1,Omega1_2,col=spp,shape=spp_sym))+
  geom_point(size=5)+
  coord_equal()+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  xlim(-2,2)+ylim(-2,2)+
  scale_color_npg()+
  guides(shape='none')+
  theme(legend.position = c(0.8,0.2))+
  labs(title="Spatial Variation, Encounter Probability",x="Factor 1",y="Factor 2",col='')
pca_o1_plot
pca_o2_plot <- pca %>% 
  ggplot(aes(Omega2_1,Omega2_2,col=spp,shape=spp_sym))+
  geom_point(size=5)+
  coord_equal()+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  xlim(-2,2)+ylim(-2,2)+
  scale_color_npg()+
  guides(col='none',shape='none')+
  labs(title="Spatial Variation, Positive Abundance",x="Factor 1",y="Factor 2",col='')
pca_o2_plot
pca_e1_plot <- pca %>% 
  ggplot(aes(Epsilon1_1,Epsilon1_2,col=spp,shape=spp_sym))+
  geom_point(size=5)+
  coord_equal()+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  scale_color_npg()+
  guides(col='none',shape='none')+
  labs(title="Spatiotemporal Variation, Encounter Probability",x="Factor 1",y="Factor 2",col='')
pca_e1_plot
pca_e2_plot <- pca %>% 
  ggplot(aes(Epsilon2_1,Epsilon2_2,col=spp,shape=spp_sym))+
  geom_point(size=5)+
  coord_equal()+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  scale_color_npg()+
  guides(col='none',shape='none')+
  labs(title="Spatiotemporal Variation, Positive Abundance",x="Factor 1",y="Factor 2",col='')
pca_e2_plot

library(gridExtra)
pca_combined <- grid.arrange(pca_o1_plot,pca_o2_plot,pca_e1_plot,pca_e2_plot,nrow=2)

ggsave(plot=pca_combined,h=10,w=10,filename = paste0(fp,"pca_combined.png"))

#correlations
corr_enc <-  plot_category_correlations(Report=Save$Report,Spatial_List,Extrapolation_List,dat,type='enc')
corr_enc
ggsave(plot=corr_enc,filename=paste0(fp,'correlation_encounter_rate.png'),w=6,h=6)
corr_pos <-  plot_category_correlations(Report=Save$Report,Spatial_List,Extrapolation_List,dat,type='pos')
corr_pos
ggsave(plot=corr_pos,filename=paste0(fp,'correlation_pos_abundance.png'),w=6,h=6)
corr_dens <-  plot_category_correlations(Report=Save$Report,Spatial_List,Extrapolation_List,dat,type='density')
corr_dens
ggsave(plot=corr_dens,filename=paste0(fp,'correlation_predicted_density.png'),w=6,h=6)

# plot co-encounter probability and station/density ##
# pull out encounter rate and log estimated density
enc_est <- Save$Report$R1_xcy
dims <- dim(enc_est)
nstations <- dims[1]
nyrs <- dims[3]
ncats <- dims[2]
dim(enc_est) <- c(nstations*nyrs,ncats)

# rearrange into dataframe
dens_est <- log(Save$Report$D_xcy)
dim(dens_est) <- c(nstations*nyrs,ncats)

cats <- tools::toTitleCase(levels(dat$spp))
cats <- gsub(" ","_",cats)
colnames(dens_est) <- cats
colnames(enc_est) <- cats
dens_est <- as.data.frame(dens_est) %>% mutate(type='log_dens')
enc_est <- as.data.frame(enc_est) %>% mutate(type='encounter')
dens_enc <- bind_rows(dens_est,enc_est)

enc_dens_df <- tibble(station=rep(1:nstations,nyrs*2),year=rep(rep(Year_Set,each=nstations),2)) %>% 
  bind_cols(dens_enc)

# add depth and temperature information
depth_temp <- X_gtp
dim(depth_temp) <- c(dims[1]*dims[3],2)
enc_dens_df <- enc_dens_df %>% 
  mutate(temp=rep(depth_temp[,1],2),depth=rep(depth_temp[,2],2)) %>% 
  # add depth domain category
  mutate(domain=case_when(
    depth<50 ~ "Inner",
    depth>=50 & depth<=100 ~ "Middle",
    depth>100 ~ "Outer"
  ))

enc_temp_depth <- enc_dens_df %>% 
  filter(type=='encounter') %>% 
  # joint encounter probability for immature opilio and medium cod
  mutate(opi_medcod=Opilio_Immature*Pacific_Cod)

quants <- c(0.05,0.5,0.95)
enc_temp_plot <- enc_temp_depth %>% 
  filter(domain=='Middle') %>% 
  ggplot(aes(temp,opi_medcod))+
  geom_point()+geom_quantile(quantiles=quants,size=2)+
  labs(x='Temperature',y='Joint Encounter Probability',title="Immature Crab and Medium Cod\nJoint Encounter Probability")
# facet_wrap(~domain,nrow=2)
enc_temp_plot


# correlations
e
