data<-read.csv("data/opilio.csv",header=T,skip=5)
maturity<-read.csv("data/maturity.csv")
#=========================================================
#==to set up initial conditioons for the projection
#==make spatial distribution of numbers at size by maturity in terminal year
#==============================================================
use_yr<-seq(1982,2017)
sexN		 <-2
binsize  <-5
sizes		 <-seq(27.5,132.5,binsize)
stations<-unique(data$GIS_STATION)

# haul join key (for comparing to other species, e.g. cod, collected in the same survey)
library(dplyr)
haul_join_key <- data %>% select(HAULJOIN,GIS_STATION,MID_LATITUDE,MID_LONGITUDE,AREA_SWEPT_VARIABLE) %>% 
  distinct() %>%
  rename(midlat=MID_LATITUDE,midlon=MID_LONGITUDE,station=GIS_STATION,AreaSwept_km2=AREA_SWEPT_VARIABLE)
save(haul_join_key,file="data/haul_join_key.Rdata")

#==Massage survey data into something more friendly
#==(if you have not done this before, rbind() is probably what you'd start with...don't! it takes forever.
#==making a list of lists then applying rbind() at the end is much faster)
#==I only do this once and then save the files because it takes ~15 min
list_geostat<-list(list())

counter<-1
aintdone<-0
if(aintdone==1)
{
  for(w in 1:length(use_yr))
    for(y in 1:length(stations))
    {
      temp <-  data[which(data$GIS_STATION == stations[y] & data$AKFIN_SURVEY_YEAR==use_yr[w] ),]
      
      if(nrow(temp)>0)
      {  
        for(x in 1:length(sizes))
        {
          temp_df<-data.frame(value = rep(0,8),
                              units =  rep(NA,8),
                              Year =  rep(NA,8),
                              Vessel =  rep(NA,8),
                              AreaSwept_km2 =  rep(NA,8),
                              Lat =  rep(NA,8),
                              Lon =  rep(NA,8),
                              Sex =  rep(NA,8),
                              Maturity =  rep(NA,8),
                              Size =  rep(NA,8),
                              Temp = rep(NA,8))
          
          temp_df$Year<-use_yr[w]
          temp_df$AreaSwept_km2<-temp$AREA_SWEPT_VARIABLE[1]
          temp_df$Lat<-temp$MID_LATITUDE[1]
          temp_df$Lon<-temp$MID_LONGITUDE[1]
          temp_df$Temp<-temp$GEAR_TEMPERATURE[1]
          temp_df$Maturity<-c("Immature","Immature","Mature","Mature","Immature","Immature","Mature","Mature")
          temp_df$Sex<-c("Female","Male","Female","Male","Female","Male","Female","Male")
          temp_df$units<-c(rep("kilos",4),rep("number",4))
          
          size_range_up<-sizes[x]+binsize/2  
          size_range_dn<-sizes[x]-binsize/2 
          temp_df$Size<-mean(c(size_range_dn,size_range_up))
          
          temp2<-temp[(temp$WIDTH>size_range_dn & temp$WIDTH<=size_range_up ),]
          
          #==fill in the numbers and weights by category
          if(nrow(temp2)>0)
          {
            #==mature females
            ugh<-temp2[!is.na(temp2$WIDTH) & temp2$SEX==2 & temp2$CLUTCH_SIZE>0,]
            if(nrow(ugh)>0)
            {
              temp_df$value[3]<-sum(ugh$SAMPLING_FACTOR)
              temp_df$value[7]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)
            }
            #==immature females
            ugh<-temp2[!is.na(temp2$WIDTH) & temp2$SEX==2 & temp2$CLUTCH_SIZE==0,]
            if(nrow(ugh)>0)
            {
              temp_df$value[1]<-sum(ugh$SAMPLING_FACTOR)
              temp_df$value[5]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)
            }
            
            #==males
            ugh<-temp2[!is.na(temp2$WIDTH) & temp2$SEX==1 & temp2$SHELL_CONDITION<=2,]
            ugh2<-temp2[!is.na(temp2$WIDTH) & temp2$SEX==1 & temp2$SHELL_CONDITION>2,]
            if(nrow(ugh)>0 & nrow(ugh2)>0)
            {
              #==mature
              temp_df$value[4]<-sum(ugh$SAMPLING_FACTOR)*maturity[x,2] + sum(ugh2$SAMPLING_FACTOR)*maturity[x,3] 
              temp_df$value[8]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)*maturity[x,2] + sum(ugh2$CALCULATED_WEIGHT*ugh2$SAMPLING_FACTOR)*maturity[x,3]
              #==immature
              temp_df$value[2]<-sum(ugh$SAMPLING_FACTOR)*(1-maturity[x,2]) + sum(ugh2$SAMPLING_FACTOR)*(1-maturity[x,3]) 
              temp_df$value[6]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)*(1-maturity[x,2]) + sum(ugh2$CALCULATED_WEIGHT*ugh2$SAMPLING_FACTOR)*(1-maturity[x,3])
              
            }
            
            if(nrow(ugh)>0 & nrow(ugh2)==0)
            {
              #==mature
              temp_df$value[4]<-sum(ugh$SAMPLING_FACTOR)*maturity[x,2]  
              temp_df$value[8]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)*maturity[x,2] 
              #==immature
              temp_df$value[2]<-sum(ugh$SAMPLING_FACTOR)*(1-maturity[x,2])
              temp_df$value[6]<-sum(ugh$CALCULATED_WEIGHT*ugh$SAMPLING_FACTOR)*(1-maturity[x,2]) 
            }
            
            if(nrow(ugh)==0 & nrow(ugh2)>0)
            {
              #==mature
              temp_df$value[4]<-sum(ugh2$SAMPLING_FACTOR)*maturity[x,3] 
              temp_df$value[8]<-sum(ugh2$CALCULATED_WEIGHT*ugh2$SAMPLING_FACTOR)*maturity[x,3]
              #==immature
              temp_df$value[2]<-sum(ugh2$SAMPLING_FACTOR)*(1-maturity[x,3]) 
              temp_df$value[6]<-sum(ugh2$CALCULATED_WEIGHT*ugh2$SAMPLING_FACTOR)*(1-maturity[x,3])
            }
            
          }
          
          list_geostat[[counter]] <-temp_df
          counter<-counter+1
        }
      }
      print(c(w,y))
    } 
  opi_dat <- do.call("rbind", list_geostat)
  save(opi_dat,file="longform_opilio.Rdata")  
}

#==Load the opilio data
load("data/longform_opilio.Rdata")

#==manipulate opi_dat to produce aggregate time series, then produce spatial datasets for VAST 
#==with the appropriate categories

# calculate the average densities for aggregate time series
# mature male biomass
MMB<- opi_dat[opi_dat$units=="kilos" & opi_dat$Sex=="Male" & opi_dat$Maturity=="Mature",] %>%
  group_by(Year) %>%
  summarise(MMB = mean(value/AreaSwept_km2,na.rm=T),
            MMBsd = sd(value/AreaSwept_km2,na.rm=T))

# mature female biomass
MFB<- opi_dat[opi_dat$units=="kilos" & opi_dat$Sex=="Female" & opi_dat$Maturity=="Mature",] %>%
  group_by(Year) %>%
  summarise(MFB = mean(value/AreaSwept_km2,na.rm=T),
            MFBsd = sd(value/AreaSwept_km2,na.rm=T))

# recruits
rec<- opi_dat[opi_dat$units=="kilos" & opi_dat$Sex=="Male" & opi_dat$Size<38,] %>%
  group_by(Year) %>%
  summarise(rec = mean(value/AreaSwept_km2,na.rm=T),
            recsd = sd(value/AreaSwept_km2,na.rm=T))

#==these need lognormal CIs put around them
par(mfrow=c(3,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
plot(rec$rec~rec$Year,type='l',xaxt='n',las=1)
legend('topright',bty='n',"Recruits")
#==a proxy for exploitable males
plot(MMB$MMB~MMB$Year,type='l',xaxt='n',las=1)
legend('topright',bty='n',"Mature male biomass")
#==a proxy for reproductive potential
plot(MFB$MFB~MFB$Year,type='l',las=1)
legend('topright',bty='n',"Mature female biomass")
mtext(side=2,line=3,"Mean density (kg/km^2",outer=T)


