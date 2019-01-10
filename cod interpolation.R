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