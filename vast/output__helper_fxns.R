### Helper functions to visualize VAST outputs ####
library(tidyverse)
library(sf)

# basemap
# library(rnaturalearth)
# library(rnaturalearthdata)
# ak <- ne_states(country='United States of America',geounit="Alaska",returnclass = 'sf')

ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(3571)

# ggplot theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        # axis.text=element_blank(),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)

# Plot data
plot_data <- function(Extrapolation_List, Spatial_List, Data, 
                       PlotDir = paste0(getwd(),"/"), 
                       Plot1_name = "Data_and_knots.png", 
                       Plot2_name = "Data_by_year.png") {
  extrapolation <- Extrapolation_List$Data_Extrap %>% st_as_sf(coords=c("Lon","Lat"),crs=4326) %>% 
    st_transform(3571)
  lims=st_bbox(extrapolation)
  ggplot()+
    geom_sf(data=ak)+
    geom_sf(data=extrapolation,size=0.5)+
    xlim(c(lims[1],lims[3]))+ylim(c(lims[2],lims[4]))
}