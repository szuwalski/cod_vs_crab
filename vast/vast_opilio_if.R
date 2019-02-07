### VAST for immature female snow crab ###

#### Settings ####
library(TMB)
library(VAST)
library(tidyverse)

Version = "VAST_v4_4_0"

# Spatial settings

# Stochastic partial differential equation (SPDE) with geometric anisotropy
Method <- "Mesh"
Aniso <- 1
grid_size_km <- 50
# number of 'knots'
n_x <- 300

# whether to include spatial (omega) and s-t (epsilon) variation
# for a univariate model, these are 0 ("turned off") or 1
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1,
                Epsilon2 = 1)
# is there temporal correlation in the intercepts or s-t variation?
# 0 means each year is a fixed effect, which is what we choose for now
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)

# is there overdispersion? often attributable to 'vessel effects'
# 0 means no overdispersion
OverdispersionConfig = c(Vessel = 0, VesselYear = 0)

# what distributions and link functions to use?
# first is the distribution for positive catch rates, second is the functional form for encounter probabilities
# we choose a gamma distribution and a Poisson log-link function
ObsModel = c(2, 1)

#outputs we want
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

strata.limits <- data.frame(STRATA = "All_areas")

# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Chionoecetes opilio")

# save settings and data in a list
fp = paste0(getwd(),'/vast/output/opilio_if/')
dir.create(fp)
Record = list(Version = Version, Method = Method, grid_size_km = grid_size_km,
              n_x = n_x, FieldConfig = FieldConfig, RhoConfig = RhoConfig,
              OverdispersionConfig = OverdispersionConfig, ObsModel = ObsModel,
              Region = Region,
              Species_set = Species_set, strata.limits = strata.limits)
save(Record, file = file.path(fp, "Record.RData"))
capture.output(Record, file = paste0(fp, "Record.txt"))

#### Import Data ####
load("data/longform_opilio2.Rdata")

# filter to just immature females and summarise total by station/year
dat <- opi_dat_long %>% 
  filter(maturity=="Immature",sex=="Female",units=="kilos") %>% 
  group_by(station,year,lon,lat,area_km2,temp) %>% 
  summarise(biomass=sum(value,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(spp=factor("Chionoecetes opilio"),vessel=0)

# make the extrapolatoin grid, building an object used to determine areas to extrapolate densities to when calculating indices
Extrapolation_List = make_extrapolation_info(Region = Region,strata.limits = strata.limits)

# generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`
Spatial_List = make_spatial_info(grid_size_km=grid_size_km,
                                  n_x=n_x, 
                                  Method=Method, 
                                  Lon_i=dat$lon, 
                                  Lat_i=dat$lat, 
                                  Extrapolation_List=Extrapolation_List, 
                                  randomseed=1,nstart=100,iter.max=1e3,
                                  DirPath=fp, Save_Results=FALSE)

# Add the knots to the the data
dat <- cbind(dat, "knot_i"=Spatial_List$knot_i)

#### Build and Run Model ####

# in order to estimate params, build a list of data-inputs used for param estimation
TmbData = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel_ez"=matrix(ObsModel,ncol=2), 
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
TmbList = Build_TMB_Fn("TmbData"=TmbData, 
                       "RunDir"=fp, 
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
                           savedir=fp, 
                           bias.correct=TRUE, 
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl"), newtonsteps=1 )

## bundle and save output!
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(fp,"Save.RData"))

#### Visualize and Diagnose Model Outputs####

# if loading of past model run needed...
# load(paste0(fp,"Save.RData"))
# for(i in 1:length(Save)) assign(names(Save)[i], Save[[i]])

# Spatial distribution of data
names(dat) <- c("station","Year","Lon","Lat","area_km2", "temp","Catch_KG","spp","vessel","knot_i")
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
Years2Include = which( Year_Set %in% sort(unique(dat[,'Year'])))

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
library(viridis)
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=12,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=10,color="black"),
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)

# center of gravity, effective area occupied, and temperature
cog <- range_indices$COG_Table
effarea <- range_indices$EffectiveArea_Table %>% 
  as_tibble()
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

ggsave(plot=cog_and_temp_plot, file="vast/output/opilio_if/cog_temp_effarea.png",height=8,width=6)

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
ggsave(plot=cog_plot2d, file="vast/output/opilio_if/cog2d.png")

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
ggsave(plot=cog_effarea_temp_corr_plot, file="vast/output/opilio_if/cog_effarea_temp_corr.png")
