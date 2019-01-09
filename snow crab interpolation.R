# Kriging, snow crab data
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
load('data/longform_opilio.Rdata')

# for now, let's start with one year of mature male crabs, collapsing all size classes
dat <- opi_dat %>% filter(Year==1982,Sex=='Male',Maturity=='Mature') %>% 
  group_by(Year,Lat,Lon) %>% 
  summarise(density_weight=sum(value/AreaSwept_km2,na.rm=T)) %>% 
  ungroup()
yr <- "1982"

# convert to simple features (spatial points) object
dat.sf <- st_as_sf(dat,coords=c('Lon','Lat'),crs=4326) %>% st_transform(3571)
dat.sp <- dat.sf %>% as_Spatial()
dat.ch <- st_convex_hull(st_union(dat.sf))
# 
# # so we know what we're dealing with
# dat.sf %>%
#   ggplot()+
#   geom_sf(aes(size=log(density_weight)),alpha=0.3)

#raster
r <- raster(dat.sp)
res(r) <- 1000 #1 km

#interpolate using inverse distance weighting
idm<-gstat(formula=density_weight~1,locations=dat.sp)
dat.idw <- interpolate(r,idm) %>% mask(as_Spatial(dat.ch))%>% 
  as.data.frame(xy=TRUE)

# state outline
ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(3571)

bbox <- st_bbox(dat.sf)
ggplot()+
  geom_raster(data=dat.idw,aes(x,y,fill=var1.pred/1000),na.rm=T,alpha=0.8)+
  geom_sf(data=ak,fill='gray80')+
  # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>% 
  # gplot(dat.idw)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  labs(x='',y='',fill='Density\n(MT/km2)',title=yr)+
  scale_fill_viridis(na.value=NA,option="C")

## Loop and create gif

## function to create a map like the above
interpolate_data <- function(df,yr,sex,maturity,scalelimits,save.plot=TRUE) {
  #subset data
  dat <- df %>% 
    filter(Year==yr,Sex==sex,Maturity==maturity) %>% 
    group_by(Year,Lat,Lon) %>% 
    summarise(density_weight=sum(value/AreaSwept_km2/1e6,na.rm=T)) %>% 
    ungroup()
  
  # convert to simple features (spatial points) object
  dat.sf <- st_as_sf(dat,coords=c('Lon','Lat'),crs=4326) %>% st_transform(3571)
  dat.sp <- dat.sf %>% as_Spatial()
  dat.ch <- st_convex_hull(st_union(dat.sf))
  
  #interpolate using inverse distance weighting
  idm<-gstat(formula=density_weight~1,locations=dat.sp)
  dat.idw <- interpolate(r,idm) %>% 
    mask(as_Spatial(dat.ch)) %>% 
    as.data.frame(xy=TRUE)
  
  plot.title <- paste0(maturity," ",sex,"s, ",yr)
  
  out<-ggplot()+
    geom_raster(data=dat.idw,aes(x,y,fill=var1.pred),na.rm=T,alpha=0.8)+
    geom_sf(data=ak,fill='gray80')+
    # projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs") %>% 
    # gplot(dat.idw)+
    xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
    labs(x='',y='',fill='Density\n(1000s MT/km2)',title=plot.title)+
    scale_fill_viridis(na.value=NA,option="C",limits=scalelimits)
  
  if(save.plot) {
    suppressMessages(ggsave(plot=out,filename=paste0("plots/",plot.title,".png")))
    print(paste(plot.title,"saved!"))
    }
  
  out
}

# appropriate scale limits for different categories

#mature males c(0,10)
#immature males c(0,10)
#mature females c(0,25)
#immature females c(0,10)

### Calculate and save plots (takes a long time)
# For mature males
walk(1982:2017, ~{
  interpolate_data(df=opi_dat,yr=.x,sex="Male",maturity="Mature",scalelimits=c(0,10),save.plot=TRUE)
})
# For mature females
walk(1982:2017, ~{
  interpolate_data(df=opi_dat,yr=.x,sex="Female",maturity="Mature",scalelimits=c(0,25),save.plot=TRUE)
})
# For immature males
walk(1982:2017, ~{
  interpolate_data(df=opi_dat,yr=.x,sex="Male",maturity="Immature",scalelimits=c(0,10),save.plot=TRUE)
})
# For immature females
walk(1982:2017, ~{
  interpolate_data(df=opi_dat,yr=.x,sex="Female",maturity="Immature",scalelimits=c(0,10),save.plot=TRUE)
})

## With gganimate-- we'll need a different function

