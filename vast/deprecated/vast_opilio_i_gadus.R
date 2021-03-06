### VAST relating immature snow crab and Pacific cod###

#### Settings ####
library(TMB)
library(VAST)
library(tidyverse)

# Version = get_latest_version( package="VAST" )
# VAST v5_4_0
Version = "VAST_v4_4_0"

#### Import Data ####
load("data/longform_opilio2.Rdata")
load("data/cod_dat_clean.Rdata")

fp = paste0(getwd(),'/vast/output/opilio_i_gadus/')
dir.create(fp)

# filter to just immature crabs, summarise total by station/year, and row-bind cod data
dat_opilio <- opi_dat_long %>% 
  filter(maturity=="Immature",units=="kilos") %>% 
  group_by(towid,year,lon,lat,area_km2) %>% 
  summarise(biomass=sum(value,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(spp="Chionoecetes opilio",vessel=0)
dat_gadus <- cod_dat_clean %>% 
  rename(lon=midlon,lat=midlat) %>% 
  group_by(towid,year,lon,lat,area_km2) %>% 
  summarise(biomass=sum(weight,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(spp="Gadus macrocephalus",vessel=0)

# combine and add missing zeros for missing tows
dat_combined <- bind_rows(dat_opilio,dat_gadus) %>% mutate(spp=factor(spp)) %>% 
  filter(!is.na(lat),!is.na(lon),year<2018) %>% 
  # do we fill in zeroes? for all tow/spp combinations that were previously missing
  complete(spp,nesting(year,lat,lon),fill=list(biomass=0,area_km2=0.01,vessel=0)) %>%
  ungroup() %>%
  # complete(spp,nesting(year,lat,lon),fill=list(area_km2=0.01,vessel=0))
  select(spp,year,lon,lat,area_km2,biomass,vessel)

# data check- number of zeroes and NAs
zeroes_check <- dat_combined %>% 
  group_by(year,spp) %>% 
  summarise(num_zeroes=sum(biomass==0),num_NA=sum(is.na(biomass))) %>% 
  arrange(spp,year)

# make the extrapolatoin grid, building an object used to determine areas to extrapolate densities to when calculating indices
# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Chionoecetes opilio","Gadus macrocephalus")

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
                                 Lon_i=dat_combined$lon, 
                                 Lat_i=dat_combined$lat, 
                                 Extrapolation_List=Extrapolation_List, 
                                 randomseed=1,nstart=100,iter.max=1e3,
                                 DirPath=fp, Save_Results=FALSE)

# Add the knots to the the data
dat <- dat_combined %>% mutate(knot_i=Spatial_List$knot_i)

# Spatial settings

# whether to include spatial (omega) and s-t (epsilon) variation
# for two linear predictors- 1. encounter probability, and 2. positive catch rate
# for a univariate model, these are 0 ("turned off") or 1
# for a multivariate model, can be any whole number 0:C, where C is the number of categories
# indicating the number of factors to estimate
FieldConfig = c(Omega1 = 2, Epsilon1 = 2, Omega2 = 2,
                Epsilon2 = 2)
# is there temporal correlation in the intercepts (Beta) or s-t variation (Epsilon)?
# 0: each year as fixed effect;
# 1: each year as random following IID distribution; 
# 2: each year as random following a random walk; 
# 3: constant among years as fixed effect; 
# 4: each year as random following AR1 process

RhoConfig = c(Beta1 = 2, Beta2 = 2, Epsilon1 = 4, Epsilon2 = 4)

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

VamConfig	<- c("Method"=2,"Rank"=1, "Timing" =0)

# what distributions and link functions to use?
# first is the distribution for positive catch rates, 
# second is the functional form for encounter probabilities,
# we choose a lognormal distribution for positive catch rates and a Poisson delta-model using log-link for enc prob and log-link for catch rates
BiasCorr= TRUE
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
TmbData = Data_Fn("Version"=Version, 
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
TmbList = Build_TMB_Fn("build_model"=TRUE,
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
Obj_orig = TmbList[["Obj"]]
Obj_orig$fn( Obj_orig$par )
Obj_orig$gr( Obj_orig$par )

# Bounds
if( VamConfig['Method']==2 ){
  if( VamConfig['Rank']>0 ){
    TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Psi_fr")[length(which(names(TmbList[["Lower"]])=="Psi_fr"))+1-1:VamConfig['Rank']]] = 0
    TmbList[["Upper"]][which(names(TmbList[["Lower"]])=="Psi_fr")[length(which(names(TmbList[["Lower"]])=="Psi_fr"))+1-1:VamConfig['Rank']]] = 1
  }
  if( RhoConfig['Epsilon1'] %in% c(4,5) ){
    TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Epsilon_rho1")] = 0
    TmbList[["Lower"]][which(names(TmbList[["Lower"]])=="Epsilon_rho1_f")] = 0
  }
}

# Optimize                                         #
Opt = TMBhelper::Optimize( obj=Obj_orig, 
                           lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], 
                           savedir=fp, 
                           getsd=TRUE, 
                           bias.correct=BiasCorr, 
                           newtonsteps=1, 
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=c("Index_cyl")))
Report = Obj_orig$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj_orig$env$parList(Opt$par), "Data"=dat, "Map"=Map)
save(Save, file=paste0(fp,"Save.RData"))
if("opt" %in% names(Opt)) capture.output( Opt$opt, file=paste0(fp,"parameter_estimates.txt"))




# Use gradient-based nonlinear minimizer to identify maximum likelihood esimates for fixed effects

#### Visualize and Diagnose Model Outputs####

# if loading of past model run needed...
# load(paste0(fp,"Save.RData"))
# for(i in 1:length(Save)) assign(names(Save)[i], Save[[i]])

# Spatial distribution of data
names(dat) <- c("Species","Year","Lon","Lat","area_km2","Catch_KG","vessel","knot_i")
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=dat, PlotDir=fp )

# Here I print the diagnostics generated during parameter estimation, and I confirm that 
# (1) no parameter is hitting an upper or lower bound and 
# (2) the final gradient for each fixed-effect is close to zero. 
# For explanation of parameters, please see `?Data_Fn`.

pander::pandoc.table(Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')])

## Encounter probability
# we check whether observed encounter frequencies for either low or high probability samples 
# are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report=Report, Data_Geostat=dat, DirName=fp)

## Fit to residuals of catch-rates (Q-Q plot)
Q = plot_quantile_diagnostic( TmbData=TmbData, 
                              Report=Report, 
                              DateFile=fp) 

## finally we can visualize residuals on a map
# Define years to plot and labels for each plotted year
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(dat[,'Year']),max(dat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(dat$Year)))

plot_residuals(Lat_i=dat[,'Lat'], 
               Lon_i=dat[,'Lon'], 
               TmbData=TmbData, 
               Report=Report, 
               Q=Q, 
               savedir=fp, 
               MappingDetails=MapDetails_List[["MappingDetails"]], 
               PlotDF=MapDetails_List[["PlotDF"]],
               MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
               Xlim=MapDetails_List[["Xlim"]], 
               Ylim=MapDetails_List[["Ylim"]], 
               FileName=fp, 
               Year_Set=Year_Set, 
               Years2Include=Years2Include, 
               Rotate=MapDetails_List[["Rotate"]], 
               Cex=MapDetails_List[["Cex"]], 
               Legend=MapDetails_List[["Legend"]], 
               zone=MapDetails_List[["Zone"]], 
               mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

### range indices- center of gravity, kernel-area occupied, and effective-area occupied ###
Sdreport <-TMB::sdreport(Obj)
range_indices <-plot_range_index(Sdreport=Sdreport,
                                 Report=Report,
                                 TmbData=TmbData,
                                 Year_Set=Year_Set,
                                 PlotDir=fp)

# figure theme
5+15+9

effarea_plot <- effarea %>% 
  ggplot(aes(Year,EffectiveArea,ymin=EffectiveArea-SE,ymax=EffectiveArea+SE))+
  geom_ribbon(alpha=0.5,fill='darkgreen')+
  geom_line()+
  scale_x_continuous(breaks=seq(1985,2015,by=5))+
  labs(y="Area Occupied (1000s km2)")
effarea_plot

cog1d <-cog %>% as_tibble() %>% 
  mutate(direction=ifelse(m==1,"eastings","northings")) %>% 
  rename(cog=COG_hat,se=SE)

cog_plot1d <- cog1d %>% 
  ggplot(aes(Year,cog))+
  geom_ribbon(aes(ymin=cog-se,ymax=cog+se),alpha=0.5,fill='red')+
  geom_line()+
  facet_wrap(~direction,scales='free_y',nrow=2)+
  scale_x_continuous(breaks=seq(1985,2015,by=5))+
  labs(y="Center of Gravity")+
  theme(axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))
cog_plot1d

# add temperature to this one, for the "middle" domain (the cold pool)
temperature_dat <- opi_dat_long %>% 
  select(year,lat,lon,temp,depth) %>% 
  distinct() %>% 
  filter(!is.na(temp),!is.na(depth)) %>% 
  mutate(domain=case_when(
    depth<50 ~ "inner",
    depth>=50 & depth<100 ~ "middle",
    depth>=100 ~ "outer"
  )) %>% 
  mutate(domain=factor(domain, levels=c("inner","middle","outer"))) %>% 
  group_by(year,domain) %>% 
  summarise(mean_temp=mean(temp,na.rm=T),sd_temp=sd(temp,na.rm=T)) %>% 
  ungroup()

# plot
mid_domain_temp_plot <- temperature_dat %>% 
  filter(domain=='middle') %>% 
  ggplot(aes(year,mean_temp,ymin=mean_temp-sd_temp,ymax=mean_temp+sd_temp))+
  geom_ribbon(alpha=0.5,fill='blue')+
  geom_line()+
  geom_hline(yintercept=0,linetype=2)+
  scale_x_continuous(breaks=seq(1985,2015,by=5))+
  labs(x='Year',y='Temp at 50-100m')+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))
mid_domain_temp_plot

# join center of gravity and temperature plots
library(gridExtra)
lay <- matrix(c(1,1,2,3),nrow=4)
cog_and_temp_plot <- grid.arrange(cog_plot1d,mid_domain_temp_plot,effarea_plot,ncol=1,layout_matrix=lay)

ggsave(plot=cog_and_temp_plot, file="vast/output/opilio_i_gadus/cog_temp_effarea.png",height=8,width=6)

cog2d<-cog %>% as_tibble() %>% 
  mutate(m=factor(m)) %>% 
  split(.$m) %>% bind_cols() %>% 
  rename(east=COG_hat,north=COG_hat1,east_se=SE,north_se=SE1) %>% 
  select(Year,east,east_se,north,north_se) %>% 
  mutate(east_upper=east+east_se,east_lower=east-east_se,
         north_upper=north+north_se,north_lower=north-north_se)

cog_plot2d <- cog2d %>% 
  ggplot(aes(east,north,ymin=north_lower,ymax=north_upper,xmin=east_lower,xmax=east_upper,col=Year))+
  geom_pointrange()+
  geom_errorbarh()+
  geom_segment(aes(xend=lead(east),yend=lead(north)),arrow=arrow(angle=20,length = unit(0.5, "cm"),type='closed'))+
  scale_color_viridis()+
  labs(x="eastings",y="northings")
cog_plot2d
ggsave(plot=cog_plot2d, file="vast/output/opilio_i_gadus/cog2d.png")

# correlations
cog_effarea_temp <- cog2d %>% 
  left_join(temperature_dat,by=c('Year'='year')) %>% 
  left_join(effarea,by='Year')

# corrs between middle domain temperature and effective area/center of gravity
library(GGally)
cog_effarea_temp_corr_plot<-ggpairs(cog_effarea_temp,columns = c(2,4,11,13),
                                    lower=list(continuous=wrap("smooth", colour="blue",alpha=0.8)),
                                    diag=list(continuous=wrap("densityDiag",fill='red',alpha=0.5)),
                                    columnLabels = c("East COG","North COG","Temperature","Eff Area"))
ggsave(plot=cog_effarea_temp_corr_plot, file="vast/output/opilio_i_gadus/cog_effarea_temp_corr.png")
