### Process cod length-frequency data ###
library(tidyverse)
library(here)
here <- here::here()

cod_dat_raw <- read_csv(here("data","race_length_by_haul_PACIFIC_COD.zip"),skip=7,col_types = cols())

# figure theme
plot_theme <-   theme_minimal()+
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=14),
        axis.title=element_text(family="sans",size=14,color="black"),
        axis.text=element_text(family="sans",size=8,color="black"),
        panel.grid.major = element_line(color="gray50",linetype=3))
theme_set(plot_theme)

# remove unnecessary columns and rename
cod_dat <- cod_dat_raw %>%
  select(Year,`Haul Join ID`,`Starting Latitude (dd)`,`Starting Longitude (dd)`,
         `Ending Latitude (dd)`,`Ending Longitude (dd)`, Stratum,
         `Bottom Depth`,`Scientific Name`,`Length (mm)`,`Frequency`,
         `Sex`,`Sample Type`,`Length Type`) %>%
  
  rename(year=Year,hauljoin=`Haul Join ID`,
         startlat=`Starting Latitude (dd)`,
         startlon=`Starting Longitude (dd)`,
         endlat=`Ending Latitude (dd)`,
         endlon=`Ending Longitude (dd)`,
         depth=`Bottom Depth`,
         spp=`Scientific Name`,
         length=`Length (mm)`,
         frequency=Frequency,
         sex=Sex,
         sample_type=`Sample Type`,
         length_type=`Length Type`) %>% 
  mutate(year=as.integer(year)) %>% 
  arrange(year)

sizecomp_by_year <- cod_dat %>% 
  group_by(year,length) %>% 
  summarise(frequency=sum(frequency)) %>%
  ungroup() %>% 
  group_by(year) %>% 
  mutate(prop=frequency/sum(frequency))

# add decade
cod_dat <- cod_dat %>%
  mutate(decade=case_when(
    year<1990 ~ "1980s",
    year<2000&year>1989 ~ "1990s",
    year<2010&year>1999 ~ "2000s",
    year>2009 ~ "2010s"
  ))

sizecomp_by_year <- sizecomp_by_year %>%
  mutate(decade=case_when(
    year<1990 ~ "1980s",
    year<2000&year>1989 ~ "1990s",
    year<2010&year>1999 ~ "2000s",
    year>2009 ~ "2010s"
  ))

## Visualize
library(ggridges)

sizecomp_all<- sizecomp_by_year %>% 
  ggplot(aes(x=length,y=factor(year),size=prop*100))+
  scale_size_continuous(range=c(1,15),breaks=c(1,2,5,10))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length (mm)",y="Year",title="Size Composition",size="Percent")+
  xlim(0,1000)+
  facet_wrap(~decade,scales='free_y')
sizecomp_all
ggsave(plot=sizecomp_all,"plots/cod_size_comp_bubble.png",height=8,width=10)

sizecomp_density_all<-sizecomp_by_year %>% 
  ungroup() %>% 
  ggplot(aes(x=length,y=factor(year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length (mm)",y="Year",title="Size Composition")+
  xlim(0,1000)+
  facet_wrap(~decade,scales='free_y')
sizecomp_density_all
ggsave(plot=sizecomp_density_all,"plots/cod_size_comp_density.png",height=8,width=10)
