# Interpolation, snow crab data
library(tidyverse)
library(gstat)
library(sf)
library(raster)
library(rasterVis)
library(rgdal)
library(viridis)

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

# for now, let's start with one year of mature male crabs, collapsing all size classes
dat <- opi_dat_long %>% filter(year==1996,sex=='Female',maturity=='Immature') %>% 
  group_by(year,lat,lon) %>% 
  summarise(density_weight=sum(value/area_km2/1e6,na.rm=T)) %>% 
  ungroup()
yr <- "1996"

# convert to simple features (spatial points) object
dat.sf <- st_as_sf(dat,coords=c('lon','lat'),crs=4326) %>% st_transform(3571)
dat.sp <- dat.sf %>% as_Spatial()

# convex hull of entire survey area

dat.ch <- opi_dat_long %>% distinct(station,lat,lon) %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
  st_transform(3571) %>% 
  st_union() %>% 
  st_convex_hull()
# 
# # so we know what we're dealing with
# dat.sf %>%
#   ggplot()+
#   geom_sf(aes(size=log(density_weight)),alpha=0.3)

#raster
r <- raster(dat.ch %>% as_Spatial())
res(r) <- 10000 #1 km

#interpolate using inverse distance weighting
idm<-gstat(formula=density_weight~1,locations=dat.sp)
dat.idw <- interpolate(r,idm) %>% mask(as_Spatial(dat.ch))%>% 
  as.data.frame(xy=TRUE)

# state outline
ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(3571)

bbox <- st_bbox(dat.ch)
bbox <- bbox*c(0.9,1.1,1.1,0.9) #expand by 10%
ggplot()+
  geom_raster(data=dat.idw,aes(x,y,fill=var1.pred),na.rm=T,alpha=0.8,interpolate=TRUE)+
  geom_sf(data=ak,fill='gray80')+
  # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>% 
  # gplot(dat.idw)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(MT/km2)',title=yr)+
  scale_fill_viridis(na.value=NA,option="C")

## Loop and create gif

## function to create a map like the above
interpolate_data <- function(df,yr,sex,maturity) {
  #subset data
  dat <- df %>% 
    filter(year==yr,sex==sex,maturity==maturity) %>% 
    group_by(year,lat,lon) %>% 
    summarise(density_weight=sum(value/area_km2/1e6,na.rm=T)) %>% 
    ungroup()
  
  # convert to simple features (spatial points) object
  dat.sf <- st_as_sf(dat,coords=c('lon','lat'),crs=4326) %>% st_transform(3571)
  dat.sp <- dat.sf %>% as_Spatial()
  
  #interpolate using inverse distance weighting
  idm<-gstat(formula=density_weight~1,locations=dat.sp)
  dat.idw <- interpolate(r,idm) %>% 
    
    # mask with convex hull from above
    mask(as_Spatial(dat.ch)) %>% 
    as.data.frame(xy=TRUE) %>% 
    mutate(year=yr)
  
  dat.idw
}

# appropriate scale limits for different categories

#mature males c(0,10)
#immature males c(0,10)
#mature females c(0,25)
#immature females c(0,10)

### Calculate and save plots (takes a long time)
# For mature males
purrr::map(1982:2017, ~{
  interpolate_data(df=opi_dat_long,yr=.x,sex="Male",maturity="Mature")
}) %>% bind_rows() -> mm.interpolated
# For mature females
purrr::map(1982:2017, ~{
  interpolate_data(df=opi_dat_long,yr=.x,sex="Female",maturity="Mature")
}) %>% bind_rows() -> mf.interpolated
# For immature males
purrr::map(1982:2017, ~{
  interpolate_data(df=opi_dat_long,yr=.x,sex="Male",maturity="Immature")
}) %>% bind_rows() -> im.interpolated
# For immature females
purrr::map(1982:2017, ~{
  interpolate_data(df=opi_dat_long,yr=.x,sex="Female",maturity="Immature")
}) %>% bind_rows() -> if.interpolated

## With gganimate
library(gganimate)
mm.gif<-mm.interpolated %>% 
  ggplot()+
    geom_raster(aes(x,y,fill=var1.pred),na.rm=T,alpha=0.8,interpolate = TRUE)+
    geom_sf(data=ak,fill='gray80')+
    xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
    labs(x='',y='',fill='Density\n(1000s MT/km2)',title='Mature Males, {frame_time}')+
    scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))+
    theme(plot.title = element_text(size=16))+
  
    transition_time(year)

# mm.interpolated %>% filter(year==1991) %>% 
# ggplot()+
#   geom_raster(aes(x,y,fill=var1.pred),na.rm=T,alpha=0.8,interpolate=TRUE)+
#   geom_sf(data=ak,fill='gray80')+
#   # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>%
#   # gplot(dat.idw)+
#   xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
#   labs(x='',y='',fill='Density\n(1000s MT/km2)',title=yr)+
#   scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))

animate(mm.gif,fps=1,nframes=length(unique(mm.interpolated$year)),width=800,height=600)
anim_save(filename="plots/gifs/mature_males.gif")

mf.gif<-mf.interpolated %>% 
  ggplot()+
  geom_raster(aes(x,y,fill=var1.pred,frame=year),na.rm=T,alpha=0.8,interpolate = TRUE)+
  geom_sf(data=ak,fill='gray80')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(1000s MT/km2)',title='Mature Females, {frame_time}')+
  scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))+
  theme(plot.title = element_text(size=16))+
  
  transition_time(year)

animate(mf.gif,fps=1,nframes=length(unique(mm.interpolated$year)),width=800,height=600)
anim_save(filename="plots/gifs/mature_females.gif")

im.gif<-im.interpolated %>% 
  ggplot()+
  geom_raster(aes(x,y,fill=var1.pred,frame=year),na.rm=T,alpha=0.8,interpolate = TRUE)+
  geom_sf(data=ak,fill='gray80')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(1000s MT/km2)',title='Immature Males, {frame_time}')+
  scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))+
  theme(plot.title = element_text(size=16))+
  
  transition_time(year)

animate(im.gif,fps=1,nframes=length(unique(mm.interpolated$year)),width=800,height=600)
anim_save(filename="plots/gifs/immature_males.gif")

if.gif<-if.interpolated %>% 
  ggplot()+
  geom_raster(aes(x,y,fill=var1.pred,frame=year),na.rm=T,alpha=0.8,interpolate = TRUE)+
  geom_sf(data=ak,fill='gray80')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(1000s MT/km2)',title='Immature Females, {frame_time}')+
  scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))+
  theme(plot.title = element_text(size=16))+
  
  transition_time(year)

animate(if.gif,fps=1,nframes=length(unique(mm.interpolated$year)),width=800,height=600)
anim_save(filename="plots/gifs/immature_females.gif")

