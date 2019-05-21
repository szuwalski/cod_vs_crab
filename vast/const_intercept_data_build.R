# Build TMB for 2.1.0 version of VAST
### VAST relating immature and spawner crabs and size-structured Pacific cod###
### We do this based on numbers (not biomass) ###
## with temperature and depth as covariates ##

#### Settings ####
library(TMB)
library(VAST)
library(tidyverse)
map <- maps::map

fp = paste0(getwd(),'/vast/output/const_intercept/')
dir.create(fp)

model_directory <- paste0(getwd(),'/vast/output/const_intercept/model/')
dir.create(model_directory)

####
Version = get_latest_version( package="VAST" )
# VAST v7_0_0
# Version = "VAST_v4_4_0"

#### Import Data ####
load("data/longform_opilio2.Rdata")
load("data/cod_dat_size_number.Rdata")

# size frequency of immature crabs
opi_imm_n <- opi_dat_long %>% filter(maturity=="Immature",units=="numbers") %>% 
  group_by(size) %>% 
  summarise(totn=sum(value)) %>% 
  ungroup()
opi_imm_n %>% ggplot(aes(size,totn))+geom_bar(stat='identity')+labs(x="Carapace Width",y="Total Number",title="Immature Snow Crab Size Frequency")

# filter for immature or mature female crabs, summarise total numbers by station/year, and row-bind cod data
# In cod stomachs, 95% of crabs range between 8 and 57mm carapace width. So we will also use this as a cutoff
dat_opilio <- opi_dat_long %>% 
  mutate(spp=case_when(
    maturity=="Immature"&size<58 ~ "Opilio immature",
    maturity=="Mature" & sex=="Female" ~ "Opilio spawner",
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(spp),units=="numbers") %>% 
  group_by(spp,year,lon,lat,area_km2) %>% 
  summarise(abun=round(sum(value,na.rm=T))) %>% 
  ungroup() %>% 
  mutate(vessel=0)

# For cod- split by 3 size classes (0-200mm FL, 200-800mm, >800mm) following Burgos et al. (2013)
# for cod sizes expected to not eat crab, eat crab often, and eat crab rarely, respectively
dat_gadus <- cod_dat %>% 
  rename(lon=midlon,lat=midlat) %>% 
  mutate(size_class=case_when(
    length<=200 ~ "small cod",
    length>200 & length<=800 ~ "medium cod",
    length>800 ~ "large cod"
  )) %>% 
  group_by(size_class,year,lon,lat) %>% 
  summarise(abun=sum(frequency,na.rm=T)) %>% 
  ungroup() %>% 
  # add area swept for reference
  mutate(spp=size_class,vessel=0,area_km2=0.01) %>% 
  select(spp,year,lon,lat,area_km2,abun,vessel)

# combine and add missing zeros for missing tows
dat_combined <- bind_rows(dat_opilio,dat_gadus) %>%
  mutate(spp=factor(spp,levels=c("Opilio immature","Opilio spawner","small cod","medium cod","large cod"))) %>% 
  filter(!is.na(lat),!is.na(lon),year<2018) %>% 
  # do we fill in zeroes? for all tow/spp combinations that were previously missing
  complete(spp,nesting(year,lat,lon),fill=list(abun=0,area_km2=0.01,vessel=0)) %>%
  ungroup() %>%
  # complete(spp,nesting(year,lat,lon),fill=list(area_km2=0.01,vessel=0))
  select(spp,year,lon,lat,area_km2,abun,vessel)

# data check- number of zeroes and NAs
zeroes_check <- dat_combined %>% 
  group_by(year,spp) %>% 
  summarise(num_zeroes=sum(abun==0),num_NA=sum(is.na(abun)),total_n=n(),perc=num_zeroes/total_n) %>% 
  arrange(spp,year)

# abundance (number) across all years, cod
# cod_hist <- dat_combined %>% filter(spp %in% c("small cod","medium cod","large cod")) %>% 
#   ggplot(aes(abun,fill=spp))+
#   geom_histogram(position='dodge',binwidth=1)+
#   xlim(0,50)
# cod_hist

# make the extrapolatoin grid, building an object used to determine areas to extrapolate densities to when calculating indices
# species 
Region = "Eastern_Bering_Sea"
Species_set = c("Opilio immature","Opilio spawner","small cod","medium cod","large cod")

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
# we choose a normal distribution for positive catch rates and a conventional delta-model for enc prob
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
load('data/depth_interpolated.RData')
load('data/nbt_interpolated.RData')

make_density_covariates <- function(type="loc_x") {
  # locations of knots from Spatial_List object
  loc_query <- Spatial_List[[type]]
  # years
  yrs <- seq(min(dat$year),max(dat$year))
  # locations of interpolated temperature and depth info
  loc_covars <- nbt.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_query,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_xtp <- array(data=NA,dim=c(nrow(loc_query),length(yrs),2))
  # fill in covariates for each knot in each year
  # depth is constant across years
  for(j in 1:length(yrs)) {
    yr_temp <- nbt.interpolated %>% filter(year==yrs[j])
    X_xtp[,j,1] <- yr_temp$temp[which_nn]
    X_xtp[,j,2] <- depth.interpolated$depth[which_nn]
  }
  return(X_xtp)
}

X_gtp <- make_density_covariates(type="loc_x")
X_itp <- make_density_covariates(type="loc_i")

#### Build and Run Model ####

# in order to estimate params, build a list of data-inputs used for param estimation
TmbData = make_data("Version"=Version, 
                    "FieldConfig"=FieldConfig,
                    "OverdispersionConfig"=OverdispersionConfig, 
                    "RhoConfig"=RhoConfig, 
                    "VamConfig" = VamConfig,
                    "ObsModel"=ObsModel,
                    "spatial_list"=Spatial_List,
                    "c_iz"=as.numeric(dat$spp)-1,
                    "b_i"=dat$abun, 
                    "a_i"=dat$area_km2, 
                    "v_i"=as.numeric(dat$vessel)-1, 
                    "s_i"=dat$knot_i-1, 
                    "t_iz"=dat$year,
                    "a_xl"=Spatial_List$a_gl,
                    "X_gtp"=X_gtp,
                    "X_itp"=X_itp,
                    "MeshList"=Spatial_List$MeshList,
                    "GridList"=Spatial_List$GridList,
                    "Method"=Spatial_List$Method,
                    "Aniso"=Aniso,
                    "Options"=Options)