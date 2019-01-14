# Compile cod data
library(tidyverse)
library(here)

here <- here::here()

cod_dat_raw <- read_csv(here("data","cod.csv"),col_types='ciiiiiiiddddiiddccciddccddiccddiic')

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
  select(Year,`Haul Join ID`,`Starting Latitude (dd)`,`Starting Longitude (dd)`,`Ending Latitude (dd)`,`Ending Longitude (dd)`, Stratum,
         `Satisfactory Gear Performance`,`Bottom Depth`,`Weight (kg)`,`Number of Fish`) %>%
  rename(hauljoin=`Haul Join ID`,startlat=`Starting Latitude (dd)`,startlon=`Starting Longitude (dd)`,
          endlat=`Ending Latitude (dd)`,endlon=`Ending Longitude (dd)`,gear_satisfactory=`Satisfactory Gear Performance`,depth=`Bottom Depth`,
          weight=`Weight (kg)`,number=`Number of Fish`
          )


# Join appropriate GIS stations (cross-referenced from snow crab data)
load('data/haul_join_key.Rdata')

cod_dat <- cod_dat %>% 
  left_join(haul_join_key,by=c("hauljoin"="HAULJOIN"))

# Aggregate by station and year
cod_dat_clean <- cod_dat %>%
  select(-hauljoin) %>%
  
  #remove all data with no station ID (HAVE TO DEAL WITH THIS LATER)
  filter(!is.na(station)) %>%
  select(station,Year,AreaSwept_km2,midlon,midlat,weight,number) %>%
  
  #density in numbers and weight
  mutate(dens_weight=weight/AreaSwept_km2,dens_num=number/AreaSwept_km2)

save(cod_dat_clean,file="data/cod_dat_clean.Rdata")

# plot total biomass over time

cod_by_yr <-cod_dat_clean %>%
  group_by(Year) %>%
  summarise(n_obs=n(),
            num=mean(dens_num,na.rm = T),
            num_sd=sd(dens_num,na.rm = T),
            weight=mean(dens_weight,na.rm=T),
            weight_sd=sd(dens_weight,na.rm=T)) %>%
  ungroup()

weight_by_yr_plot <-cod_by_yr %>%
  ggplot(aes(Year,weight))+
  geom_line()+
  geom_point()+
  labs(x='year',y='weight (kg/km2)')
num_by_yr_plot <-cod_by_yr %>%
  ggplot(aes(Year,num))+
  geom_line(col='blue')+
  geom_point(col='blue')+
  labs(x='year',y='density (number/km2)')

library(gridExtra)

# compare visually to assessment data
assessment_abundance <- read_csv(here("data","assessment_tot_abun.csv"),col_types = 'dddddd') %>%
  ggplot(aes(Year,Estimate))+
  geom_point(col='darkgreen')+geom_line(col='darkgreen')+
  labs(x='year',y='Abundance (1000s fish)')+
  theme(panel.border = element_rect(fill=NA,color="darkgreen",size=2))


grid.arrange(weight_by_yr_plot,num_by_yr_plot,assessment_abundance)

#size composition (from assessment)
sizecomp <- read_csv('data/assessment_size_comp.csv') %>% 
  gather(cm,prop,-Year,-N) %>%
  mutate(prop=prop/N) %>% 
  filter(cm!="120+") %>% 
  mutate(cm=as.numeric(cm),decade=case_when(
    Year<1990 ~ "1980s",
    Year<2000&Year>1989 ~ "1990s",
    Year<2010&Year>1999 ~ "2000s",
    Year>2009 ~ "2010s"
  ))

#match to bins
# binsize=5
# sizes		 <-seq(27.5,132.5,binsize)
# # size bin matching key
# sizesdf <- tibble(dn=seq(25,130,5),up=seq(30,135,5),size=sizes)

library(ggridges)
sizecomp %>% 
  ungroup() %>% 
  filter(Year<1990) %>% 
  ggplot(aes(x=cm,y=factor(Year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length",y="Year",title="Size Composition, 1980s")+
  xlim(0,100)

sizecomp %>% 
  ungroup() %>% 
  filter(Year<2000,Year>1989) %>% 
  ggplot(aes(x=cm,y=factor(Year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length",y="Year",title="Size Composition, 1990s")+
  xlim(0,100)
sizecomp %>% 
  ungroup() %>% 
  filter(Year<2010,Year>1999) %>% 
  ggplot(aes(x=cm,y=factor(Year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length",y="Year",title="Size Composition, 2000s")+
  xlim(0,100)
sizecomp %>% 
  ungroup() %>% 
  filter(Year>2009) %>% 
  ggplot(aes(x=cm,y=factor(Year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length",y="Year",title="Size Composition, 2010s")+
  xlim(0,100)
#all
sizecomp_density_all<-sizecomp %>% 
  ungroup() %>% 
  ggplot(aes(x=cm,y=factor(Year),height=prop))+
  geom_density_ridges(stat='identity',fill='blue',alpha=0.5,scale=2,color='black')+
  labs(x="Length",y="Year",title="Size Composition")+
  xlim(0,100)+
  facet_wrap(~decade,scales='free_y')
ggsave(plot=sizecomp_density_all,"plots/cod_size_comp_density.png",height=8,width=10)

# bubble?
sizecomp %>% 
  ungroup() %>% 
  filter(Year<1990) %>% 
  ggplot(aes(x=cm,y=factor(Year),size=prop))+
  scale_size_continuous(range=c(1,20))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length",y="Year",title="Size Composition, 1980s")+
  xlim(0,100)

sizecomp %>% 
  ungroup() %>% 
  filter(Year<2000,Year>1989) %>% 
  ggplot(aes(x=cm,y=factor(Year),size=prop))+
  scale_size_continuous(range=c(1,20))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length",y="Year",title="Size Composition, 1990s")+
  xlim(0,100)
sizecomp %>% 
  ungroup() %>% 
  filter(Year<2010,Year>1999) %>% 
  ggplot(aes(x=cm,y=factor(Year),size=prop))+
  scale_size_continuous(range=c(1,20))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length",y="Year",title="Size Composition, 2000s")+
  xlim(0,100)
sizecomp %>% 
  ungroup() %>% 
  filter(Year>2009) %>% 
  ggplot(aes(x=cm,y=factor(Year),size=prop*100))+
  scale_size_continuous(range=c(1,20))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length",y="Year",title="Size Composition, 2010s",size="Percent")+
  xlim(0,100)
#all
sizecomp_bubble_all<-sizecomp %>% 
  ungroup() %>% 
  ggplot(aes(x=cm,y=factor(Year),size=prop*100))+
  scale_size_continuous(range=c(1,20))+
  geom_point(alpha=0.3,color='blue')+
  labs(x="Length",y="Year",title="Size Composition",size="Percent")+
  xlim(0,100)+
  facet_wrap(~decade,scales='free_y')
ggsave(plot=sizecomp_bubble_all,"plots/cod_size_comp_bubble.png",height=8,width=10)
