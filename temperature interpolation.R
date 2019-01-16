# Temperature interpolation
# Based on NBT in the NMFS Summer Trawl survey
# Interpolation, snow crab data
library(raster)
library(rasterVis)
library(rgdal)
library(tidyverse)
library(gstat)
library(sf)
library(viridis)

select <- dplyr::select

# figure theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)

# import data
load('data/longform_opilio2.Rdata')

nbt_dat <- opi_dat_long %>% 
  select(year,lat,lon,temp) %>% 
  distinct() %>% 
  filter(!is.na(temp))

# convert to simple features (spatial points) object
# transform to a good Bering Sea projection (EPSG code 3571, Lambert azimuthal equal-area)
dat.sf <- st_as_sf(nbt_dat,coords=c('lon','lat'),crs=4326) %>% st_transform(3571) 
dat.sp <- dat.sf %>% as_Spatial()

# convex hull of entire survey area
dat.ch <- st_convex_hull(dat.sf %>% st_union())


#raster
r <- raster(dat.ch %>% as_Spatial())
res(r) <- 10000 #1 km

# state outline
ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(3571)

bbox <- st_bbox(dat.ch)
bbox <- bbox*c(0.9,1.1,1.1,0.9) #expand by 10%

#interpolate using inverse distance weighting
idm<-gstat(formula=temp~1,locations=dat.sf %>% filter(year==1995) %>% as_Spatial())
dat.idw <- interpolate(r,idm) %>% 
  mask(as_Spatial(dat.ch)) %>% 
  as.data.frame(xy=TRUE) %>% 
  mutate(year=1995)
# # Test plot
ggplot()+
  geom_raster(data=dat.idw,aes(x,y,fill=var1.pred),na.rm=T,alpha=0.8,interpolate=TRUE)+
  geom_sf(data=ak,fill='gray80')+
  # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>%
  # gplot(dat.idw)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Temperature',title="Test:NBT")+
  scale_fill_viridis(na.value=NA,option="C")

## function to create a map like the above
interpolate_data <- function(df,yr) {
  #subset data
  dat <- df %>% 
    filter(year==yr)
  
  # convert to simple features (spatial points) object
  dat.sf <- st_as_sf(dat,coords=c('lon','lat'),crs=4326) %>% st_transform(3571)
  dat.sp <- dat.sf %>% as_Spatial()
  
  #interpolate using inverse distance weighting
  idm<-gstat(formula=temp~1,locations=dat.sp)
  dat.idw <- interpolate(r,idm) %>% 
    mask(as_Spatial(dat.ch)) %>% 
    as.data.frame(xy=TRUE) %>% 
    mutate(year=yr)
  
  dat.idw
}

### Calculate and save plots (takes a long time)
# For mature males
purrr::map(1982:2017, ~{
  interpolate_data(df=nbt_dat,yr=.x)
}) %>% bind_rows() -> nbt.interpolated


## With gganimate
library(gganimate)
nbt.gif<-nbt.interpolated %>% 
  ggplot()+
  geom_raster(aes(x,y,fill=var1.pred,frame=year),na.rm=T,alpha=0.8,interpolate = TRUE)+
  geom_sf(data=ak,fill='gray80')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Near Bottom\nTemperature (deg C)',title='Temperature, {frame_time}')+
  scale_fill_viridis(na.value=NA,option="C",limits=c(-2,10))+
  theme(plot.title = element_text(size=16))+
  
  transition_time(year)

animate(nbt.gif,fps=1,nframes=length(unique(nbt_dat$year)),width=800,height=600)
anim_save(filename="plots/gifs/nbt.gif")

# With a time series by depth band ()
nbt_ts <- opi_dat_long %>% 
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
  summarise(mean_temp=mean(temp,na.rm=T)) %>% 
  ungroup()
  