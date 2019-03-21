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
  st_transform(3571)

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
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)

# Plot predicted density
matxt <- plot_maps()
plot_data <- function(Region,Spatial_List, Extrapolation_List,dat) {
  
  Report <- Save$Report
  
  dens <- log(Report$D_xcy)
  
  pts <- Extrapolation_List$Data_Extrap[,c('E_km','N_km','Area_in_survey_km2')] %>% 
    st_as_sf(coords=c('E_km','N_km'))
  pts$idx <-as.integer(Spatial_List$PolygonList$NN_Extrap$nn.idx)
  pts <- filter(pts,Area_in_survey_km2>0) %>% select(-Area_in_survey_km2)
  Year_Set = seq(min(dat$Year),max(dat$Year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$Year)))]
  knots <- dim(dens)[1]
  
  for(i in 1:dim(dens)[2]){
    df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=as.numeric(dens[,i,]))
    sp_df<-pts %>% left_join(df,by=c('idx'='knot'))
    ggplot(sp_df,aes(E_km,N_km,col=dens))+geom_point()+coord_equal()+scale_color_viridis(limits=c(-5,13))+
      facet_wrap(~year)+
      labs(x="Eastings",y="Northings",col="")+
      theme(axis.text = element_blank())
  }
  
}

# Plot factor maps and factor loadings by category (species)

# Plot species correlations