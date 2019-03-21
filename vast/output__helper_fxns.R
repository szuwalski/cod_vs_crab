### Helper functions to visualize VAST outputs ####
### specific to Eastern Bering Sea ###

library(tidyverse)
library(sf)
library(VAST)

# basemap
# library(rnaturalearth)
# library(rnaturalearthdata)
# ak <- ne_states(country='United States of America',geounit="Alaska",returnclass = 'sf')

ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(26904)

## data for testing ##
fp = paste0(getwd(),'/vast/output/opilio_gadus_numbers/')
load(paste0(fp,"Record.Rdata"))
load(paste0(fp,"Save.RData"))
load(paste0(fp,"Spatial_List.Rdata"))
load(paste0(fp,"parameter_estimates.Rdata"))
dat <- Save$Data
list2env(Record,envir = environment())

names(dat) <- c("spp","Year","Lon","Lat","area_km2", "Catch_KG","vessel","knot_i")

# ggplot theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_blank(),
        panel.border=element_rect(color='black',fill=NA))
theme_set(plot_theme)

# Plot predicted density
matxt <- plot_maps()
plot_data <- function(Region,Spatial_List, Extrapolation_List,dat) {
  
  Report <- Save$Report
  
  dens <- log(Report$D_xcy)
  
  pts <- Extrapolation_List$Data_Extrap[,c('Lon','Lat','Area_in_survey_km2')] %>% 
    st_as_sf(coords=c('Lon','Lat'),crs=4326) %>% 
    #convert to AK UTM zone
    st_transform(26904)
  pts[,c("E_km","N_km")] <- st_coordinates(pts)
  pts$idx <-as.integer(Spatial_List$PolygonList$NN_Extrap$nn.idx)
  pts <- filter(pts,Area_in_survey_km2>0) %>% select(-Area_in_survey_km2)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$Year),max(dat$Year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$Year)))]
  knots <- dim(dens)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))
  
  for(i in 1:dim(dens)[2]){
    df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=as.numeric(dens[,i,]))
    sp_df<-pts %>% left_join(df,by=c('idx'='knot'))
    # make the facetted plot by year
    cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
    out<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=dens))+
      scale_color_gradientn(colors=cols)+
      facet_wrap(~year)+
      labs(title=cats[i],x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
    # save the plot
    ggsave(plot=out,filename = paste0(fp,"dens_",cats[i],".png"),w=6,h=6)
  }
  
}

# Plot factor maps and factor loadings by category (species)

# Plot species correlations