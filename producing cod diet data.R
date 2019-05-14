## Clean cod stomach content data ##
library(tidyverse)
library(here)
here <- here::here()
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)


dat <- read_csv("https://access.afsc.noaa.gov/REEM/WebDietData/showdata2.php?NODC=8791030401&Region=BS")
glimpse(dat)
names(dat)

# Notes from Livingston et al 2017
# Use Prey_Name, not Prey_NODC
# Care should be taken when scaling records for ind. fish up to population-level estimates of prey consumption (see Buckley et al 2015)
# per-haul estimates are very noisy and should be treated with caution
# 1981-1984 sample sizes were lower
# Be careful about assumptions made about empty stomachs, because regurgitation is common and varies among samplers
# Temporal aliasing of the survey should be taken into account with respect to seasonal and spatial patterns of prey availability

# Some important variable definitions (all are provided in Livingston et al. 2017 Appendix)
# Prey_cnt is count of a given prey
# Prey_twt is total weight of a given prey; value of 0.0 indicates empty stomach
# Pred_len is fork length of the predator
# Pred_full is a categorical fullness index of the predator stomach
# Prey_lh is prey life history stage (Table 6)
# Pred_dig is how digested a given prey item is

# For now, we select just predator, predator length, prey life history stage, and id codes (haul, lat/lon, etc)

dat_reduced <- dat %>% 
  select(Hauljoin,Year, Month, day, region, Pred_name,Prey_Name,Prey_lh,Rlat,Rlong) %>% 
  filter(Prey_Name=="Opilio Crab")

# Count of prey by life history stage
dat_reduced %>%
  count(Prey_lh) %>% 
  ggplot(aes(Prey_lh,n))+
  geom_bar(stat='identity')+
  labs(title="Number of Observations by Prey Category")

# The grand majority of opilio records are either 7 (Juvenile), C (unknown), or L(egg carrying female)

# Samples by year
dat_reduced %>% 
  count(Year) %>% 
  ggplot(aes(Year,n))+
  geom_point()+geom_line()+
  labs(title="Number of Predation Observations by Year")

# Maps