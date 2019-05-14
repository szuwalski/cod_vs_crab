### Goal: Bring together single-species (univariate) results 
### for cod and crab and look at spatial correlation

library(tidyverse)
# import model outputs
load("vast/output/opilio_immature/Save.RData")
opilio_imm_save <- Save$Report
load("vast/output/opilio_mf/Save.RData")
opilio_mf_save <- Save$Report
load("vast/output/gadus/Save.RData")
gadus_save <- Save
rm(Save)

# import spatial knots
load("vast/output/opilio_immature/Spatial_List.RData")
knots <- Spatial_List$loc_x

## Pull out log density estimates by spatial location and year

opilio_imm <- log(opilio_imm_save$Report$D_xcy)
opilio_mf <- log(opilio_mf_save$Report$D_xcy)
gadus <- log(gadus_save$Report$D_xcy)
