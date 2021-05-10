# Interpolation, cod data
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
load('data/cod_dat_clean.Rdata')

# test data
testdat <- cod_dat_clean %>% filter(Year==1982)

# convert to simple features (spatial points) object
# transform to a good Bering Sea projection (EPSG code 3571, Lambert azimuthal equal-area)
dat.sf <- st_as_sf(cod_dat_clean,coords=c('midlon','midlat'),crs=4326) %>% st_transform(3571) 
dat.sp <- dat.sf %>% as_Spatial()

# convex hull of entire survey area

dat.ch <- cod_dat_clean %>% distinct(station,midlat,midlon) %>% 
  st_as_sf(coords=c('midlon','midlat'),crs=4326) %>% 
  st_transform(3571) %>% 
  st_union() %>% 
  st_convex_hull()


#raster
r <- raster(dat.ch %>% as_Spatial())
res(r) <- 10000 #1 km

#interpolate using inverse distance weighting
idm<-gstat(formula=dens_weight~1,locations=dat.sp)
dat.idw <- interpolate(r,idm) %>% mask(as_Spatial(dat.ch))%>% 
  as.data.frame(xy=TRUE)

# state outline
ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(3571)

bbox <- st_bbox(dat.ch)
bbox <- bbox*c(0.9,1.1,1.1,0.9) #expand by 10%

# Test plot
ggplot()+
  geom_raster(data=dat.idw,aes(x,y,fill=var1.pred/1000),na.rm=T,alpha=0.8,interpolate=TRUE)+
  geom_sf(data=ak,fill='gray80')+
  # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>% 
  # gplot(dat.idw)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(MT/km2)',title="Test: Cod distribution 1982 interpolated")+
  scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))

## function to create a map like the above
interpolate_data <- function(df,yr) {
  #subset data
  dat <- df %>% 
    filter(Year==yr)
  
  # convert to simple features (spatial points) object
  dat.sf <- st_as_sf(dat,coords=c('midlon','midlat'),crs=4326) %>% st_transform(3571)
  dat.sp <- dat.sf %>% as_Spatial()
  
  #interpolate using inverse distance weighting
  idm<-gstat(formula=dens_weight~1,locations=dat.sp)
  dat.idw <- interpolate(r,idm) %>% 
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
  interpolate_data(df=cod_dat_clean,yr=.x)
}) %>% bind_rows() -> cod.interpolated


## With gganimate
library(gganimate)
cod.gif<-cod.interpolated %>% 
  ggplot()+
  geom_raster(aes(x,y,fill=var1.pred/1000,frame=year),na.rm=T,alpha=0.8,interpolate = TRUE)+
  geom_sf(data=ak,fill='gray80')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(MT/km2)',title='Cod Density, {frame_time}')+
  scale_fill_viridis(na.value=NA,option="C",limits=c(0,25))+
  theme(plot.title = element_text(size=16))+
  
  transition_time(year)

animate(cod.gif,fps=1,nframes=length(unique(cod.interpolated$year)),width=800,height=600)
anim_save(filename="plots/gifs/cod.gif")
