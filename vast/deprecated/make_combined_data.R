# Quick script to create clean cod and crab data ready for use in VAST models
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

dat_combined= as.data.frame(dat_combined)
## We also make the density covariates
# Temperature and depth
# this will go into VAST::make_data as a 3 dimensional array
# with dimensions nknots x nyears x nvars
# We use temperature and depth data interpolated in another script
library(RANN)
library(VAST)
load('data/depth_interpolated.RData')
load('data/nbt_interpolated.RData')


make_density_covariates <- function(dat,Region,strata.limits,grid_size_km,n_x,Method,Lon_i,Lat_i,fine_scale) {
  # make extrapolation info and spatial list (from VAST)
  Extrapolation_List = make_extrapolation_info(Region = 'eastern_bering_sea',strata.limits = data.frame(STRATA = "All_areas"))
  Spatial_List = make_spatial_info(grid_size_km=grid_size_km,n_x=n_x, Method=Method, Lon_i=dat$lon, Lat_i=dat$lat, 
                                   Extrapolation_List=Extrapolation_List, fine_scale=fine_scale,
                                   randomseed=1,nstart=100,iter.max=1e3,
                                   Save_Results=FALSE)
  
  # locations of knots from Spatial_List object
  loc_knots <- Spatial_List$loc_g
  # years
  yrs <- seq(min(dat$year),max(dat$year))
  # locations of interpolated temperature and depth info
  loc_covars <- nbt.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_knots,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_gtp <- array(data=NA,dim=c(nrow(loc_knots),length(yrs),2))
  # fill in covariates for each knot in each year
  # depth is constant across years
  for(j in 1:length(yrs)) {
    yr_temp <- nbt.interpolated %>% filter(year==yrs[j])
    X_gtp[,j,1] <- yr_temp$temp[which_nn]
    X_gtp[,j,2] <- depth.interpolated$depth[which_nn]
  }
  return(X_gtp)
}
