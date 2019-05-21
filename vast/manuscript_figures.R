# Produce manuscript figures
# 1. Factor Maps
# 2. Factor loadings (PCA)
# 3. 2-dimensional factors
# 4. Correlations between classes
# 5. Lagged correlations (time series)
# 6. Time series indices

library(tidyverse)
library(sf)
library(VAST)
library(abind)
library(viridis)
library(ggsci)
library(RANN)

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

# basemap
# library(rnaturalearth)
# library(rnaturalearthdata)
# ak <- ne_states(country='United States of America',geounit="Alaska",returnclass = 'sf')

ak <- read_sf('data/spatial/cb_2017_02_anrc_500k.shp') %>% 
  st_union() %>% 
  st_transform(26904)

## data for plotting
fp = paste0(getwd(),'/vast/output/const_intercept/')
load(paste0(fp,"derived_quantities.Rdata"))
load(paste0(fp,"Extrapolation_List.Rdata"))
load(paste0(fp,"Spatial_List.Rdata"))
load(paste0(fp,"Save.Rdata"))
Report <- Save$Report

load('data/depth_interpolated.RData')
load('data/nbt_interpolated.RData')

map_points <- function(Extrapolation_List) {
  pts <- Extrapolation_List$Data_Extrap[,c('Lon','Lat','Area_in_survey_km2')] %>% 
    st_as_sf(coords=c('Lon','Lat'),crs=4326) %>% 
    #convert to AK UTM zone
    st_transform(26904)
  pts[,c("E_km","N_km")] <- st_coordinates(pts)
  pts$idx <-as.integer(Spatial_List$PolygonList$NN_Extrap$nn.idx)
  pts <- filter(pts,Area_in_survey_km2>0) %>% select(-Area_in_survey_km2)
  
  return(pts)
}

map_density <- function(dat,Report,Extrapolation_list,category) {
  
  dens <- Report$D_xcy
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  knots <- dim(dens)[1]
  
  cats <- tools::toTitleCase(levels(dat$spp))

  df <- tibble(knot=rep(1:knots,length(Years2Include)),year=rep(Years2Include,each=knots),dens=log(as.numeric(dens[,category,])))
  sp_df<-pts %>% left_join(df,by=c('idx'='knot')) %>% 
    group_by(idx) %>% sample_n(500)
  # make the facetted plot by year
  # cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  out<-ggplot()+
    geom_point(data=sp_df,aes(E_km,N_km,col=dens))+
    scale_color_gradientn(colors=cols)+
    facet_wrap(~year)+
    labs(title=cats[category],x="Eastings",y="Northings",col="Log Numbers\nDensity")+
    geom_sf(data=ak,fill='gray50',col=NA)+
    coord_sf(crs = st_crs(ak), datum = NA) +
    xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
    # scale_color_viridis()+
    theme(axis.text = element_blank(),
          panel.spacing.x = unit(0,"pt"),
          panel.spacing.y=unit(0,"pt"))

  return(out)
}

# Plot factor maps and factor loadings by category (species)
plot_fct_loadings <- function(factors,dat,rotated=TRUE) {
  
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  # collect factor loadings for the different predictors
  if(rotated){
    loadings <- factors$Rotated_loadings
  }
  else{
    loadings <- factors$Loadings
  }
  
  loadings_lng <- map_df(names(loadings),.f = function(x){
    m <- loadings[[x]]
    if(is.matrix(m)){
      tibble(load=as.numeric(m),spp=rep(cats,ncol(m)),fct=x,fct_num=rep(1:ncol(m),each=nrow(m)))
    }
  })
  # # proportion of variance explained
  # round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1)
  fct_labeller <- c(
    "Omega1" = "Spatial Variation\nEncounter",
    "Omega2" = "Spatial Variation\nPositive Abundance",
    "Epsilon1" = "Spatio-Temporal Variation\nEncounter",
    "Epsilon2" = "Spatio-temporal Variation\nPositive Abudance"
  )
  # indicator for pos/neg
  loadings_lng <- loadings_lng %>% 
    # proportion of variance explained
    group_by(fct,fct_num) %>% 
    mutate(loadsq=sum(load^2)) %>%
    ungroup() %>% group_by(fct) %>% 
    mutate(tot=sum(load^2),prop=round(loadsq/tot*100,1)) %>% 
    ungroup() %>% 
    mutate(is_pos=ifelse(load>0,"yes","no"))
  out <-ggplot(loadings_lng,aes(spp,load,fill=is_pos))+
    geom_bar(stat='identity')+
    geom_hline(yintercept=0,col='black')+
    geom_text(aes(x=4.5,y=1.5,label=prop),check_overlap = T)+
    scale_fill_npg()+
    facet_grid(fct_num~fct,labeller = labeller(fct=fct_labeller))+
    guides(fill='none')+
    labs(x="Species",y="Loading",title="Factor Loadings")+
    theme(axis.text.x = element_text(angle=60,hjust=1))
  out
}
# plot 2-dimensional loadings
plot_2d_fct_loadings <- function(dat,factors,rotated=TRUE,which_lpred){
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  # collect factor loadings for the different predictors
  if(rotated){
    loadings <- factors$Rotated_loadings
  }
  else{
    loadings <- factors$Loadings
  }
  
  loadings_lng <- map_df(names(loadings),.f = function(x){
    m <- loadings[[x]]
    if(is.matrix(m)){
      tibble(load=as.numeric(m),spp=rep(cats,ncol(m)),fct=x,fct_num=rep(1:ncol(m),each=nrow(m)))
    }
  })
  # # proportion of variance explained
  # round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1)
  
  # indicator for pos/neg
  loadings_lng <- loadings_lng %>% 
    # proportion of variance explained
    group_by(fct,fct_num) %>% 
    mutate(loadsq=sum(load^2)) %>%
    ungroup() %>% group_by(fct) %>% 
    mutate(tot=sum(load^2),prop=round(loadsq/tot*100,1)) %>% 
    ungroup() %>% 
    mutate(is_pos=ifelse(load>0,"yes","no"))

  pca <- fct_loadings %>% 
    filter(fct_num < 3) %>% 
    select(spp,fct,fct_num,load) %>%
    unite("fct",fct,fct_num) %>% 
    group_by(spp) %>% 
    spread(fct,load) %>% 
    mutate(spp_sym= ifelse(spp %in% c('Opilio Immature','Opilio Spawner'),'crab','cod')) %>% 
    select(spp,spp_sym,contains(which_lpred)) %>% 
    ungroup()
  xmax <- ceiling(max(pca[,3]))
  xmin <- floor(min(pca[,3]))
  ymax <- ceiling(max(pca[,4]))
  ymin <- floor(min(pca[,4]))
  
  varnames <- names(pca)
  
  pca_plot <- pca %>% 
    ggplot(aes_string(varnames[3],varnames[4],col="spp",shape="spp_sym"))+
    geom_point(size=5)+
    coord_equal()+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    xlim(xmin,xmax)+ylim(ymin,ymax)+
    scale_color_npg()+
    guides(shape='none')+
    theme(legend.position = c(0.2,0.2))+
    labs(title="",x="Factor 1",y="Factor 2",col='')+
    theme(legend.text = element_text(size=8),
          legend.position = "bottom")
  pca_plot
}
plot_fct_maps <- function(Extrapolation_List,Report,dat,which_lpred,which_f){
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  Year_Set = seq(min(dat$year),max(dat$year))
  Years2Include = Year_Set[which( Year_Set %in% sort(unique(dat$year)))]
  
  #pull out factor map data, for enc prob and pos values
  fcts <- factors$Rotated_factors %>% pluck(which_lpred)
  n_knots <- dim(fcts)[1]
  fct <- fcts[,which_f,]
  n_knots<-ifelse(grepl('Omega',which_lpred),length(fct),dim(fct)[1])
  n_yrs <- ifelse(grepl('Omega',which_lpred),1,dim(fct)[2])

  df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct))
  if(grepl('Epsilon',which_lpred)) {df <- df %>% mutate(year=rep(Years2Include,each=n_knots))}
  
  cols<-colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))(50)
  
  sp_df<-pts %>% left_join(df,by=c('idx'='knot')) %>% group_by(idx) %>% sample_frac(0.5)
  
  if(grepl('Omega',which_lpred)){
    out<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=value))+
      scale_color_gradientn(colors=cols)+
      labs(title="",x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
  }
  else{
    out<-ggplot()+
      geom_point(data=sp_df,aes(E_km,N_km,col=value))+
      scale_color_gradientn(colors=cols)+
      facet_wrap(~year)+
      labs(title="",x="Eastings",y="Northings",col="")+
      geom_sf(data=ak,fill='gray50',col=NA)+
      coord_sf(crs = st_crs(ak), datum = NA) +
      xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
      theme(axis.text = element_blank(),
            panel.spacing.x = unit(0,"pt"),
            panel.spacing.y=unit(0,"pt"))
  }
  # make the facetted plot by year

    out
}

# Plot species correlations (density correlations)
cor.sig <- function(x,y){
  corr <- cor(x,y,method='spearman')
  test <- suppressWarnings(cor.test(x,y,method='spearman'))
  out <- ifelse(test$p.value<0.05,corr,NA)
  out
}

plot_category_correlations <- function(Report, Extrapolation_List,dat,type="density"){
  
  if(!(type %in% c('density','enc','pos'))){
    stop("type parameter must be one of 'density','enc','pos'")
  }
  
  if(type=='density'){
    est <- log(Report$D_xcy) 
    lab <- "Density"
  }
  if(type=='enc'){
    est <- Report$R1_xcy
    lab <- "Encounter Probability"
  }
  if(type=='pos'){
    est <- Report$R2_xcy
    lab <- "Positive Catch Rate"
  }
  
  # row bind all the years for each species
  report_dims <- dim(Report$D_xcy)
  dim(est) <- c(report_dims[1]*report_dims[3],report_dims[2])
  
  # get correlations
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  spp_cor <- numeric()
  
  for(i in 1:length(cats)){
    for(j in 1:length(cats)){
      spp_cor[(length(cats)*(i-1))+j] = cor.sig(est[,i],est[,j])
    }
  }
  spp_cor[spp_cor==1] <- NA
  
  df <- tibble(spp1=rep(cats,each=length(cats)),spp2=rep(cats,length(cats)),corr=spp_cor)
  
  out <- df %>% 
    ggplot(aes(spp1,spp2,fill=corr))+
    geom_tile()+
    scale_fill_gradient2()+
    geom_text(aes(label=round(corr,2)))+
    coord_equal()+
    labs(x="",y="",title=paste("Correlation in Predicted",lab),fill=expression(rho))+
    theme(axis.text.x = element_text(angle=60,hjust=1),
          legend.title = element_text(size=14))
  out
}

map_category_correlations <- function(Report, Extrapolation_List,dat,type="density",species_1,species_2){
  
  if(!(type %in% c('density','enc','pos'))){
    stop("type parameter must be one of 'density','enc','pos'")
  }
  
  if(type=='density'){
    est <- log(Report$D_xcy) 
    lab <- "Density"
  }
  if(type=='enc'){
    est <- Report$R1_xcy
    lab <- "Encounter Probability"
  }
  if(type=='pos'){
    est <- Report$R2_xcy
    lab <- "Positive Catch Rate"
  }
  
  # row bind all the years for each species
  # report_dims <- dim(Report$D_xcy)
  # dim(est) <- c(report_dims[1]*report_dims[3],report_dims[2])
  
  # get correlations for each knot
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  # spp_cor <- numeric()
  
  spp_corr_by_knot <- array(dim=c(dim(est)[1],dim(est)[2],dim(est)[2]))
  for(k in 1:dim(est)[1]){
    corrmat <- numeric()
    for(i in 1:length(cats)){
      for(j in 1:length(cats)){
        corrmat[(length(cats)*(i-1))+j] <- cor.sig(est[k,i,],est[k,j,])
      }
    }
    spp_corr_by_knot[k,,] <- corrmat
  }
  
  df <- tibble(knot=rep(1:knots,length(cats)*length(cats)),spp1=rep(rep(cats,each=knots),length(cats)),spp2=rep(cats,each=knots*length(cats)),corr=as.numeric(spp_corr_by_knot))
  
  df_cat <- df %>% 
    filter(spp1==species_1,spp2==species_2)
  
  pts <- map_points(Extrapolation_List)
  lims <- st_bbox(pts)
  st_geometry(pts) <- NULL
  
  sp_df<-pts %>% left_join(df_cat,by=c('idx'='knot')) %>% group_by(idx) %>% sample_frac(0.5) %>% ungroup()
  
  out<-ggplot()+
    geom_point(data=sp_df,aes(E_km,N_km,col=corr))+
    scale_color_gradient2(na.value = 'gray80')+
    labs(title="",x="Eastings",y="Northings",col="Correlation")+
    geom_sf(data=ak,fill='gray50',col=NA)+
    coord_sf(crs = st_crs(ak), datum = NA) +
    xlim(lims[1],lims[3])+ylim(lims[2],lims[4])+
    # scale_color_viridis()+
    theme(axis.text = element_blank(),
          panel.spacing.x = unit(0,"pt"),
          panel.spacing.y=unit(0,"pt"))
  out
}

plot_abundance <- function(Index_table,dat){
  abun <- Index_table %>% as.data.frame() %>% 
    mutate(log_abun=log(Estimate_metric_tons),upper=log_abun+SD_log,lower=log_abun-SD_log) %>% 
    mutate(Category=tools::toTitleCase(as.character(Category))) %>% 
    mutate(Category=factor(Category,levels=c("Opilio Immature","Opilio Spawner","Small Cod","Medium Cod","Large Cod")))
  
  ggplot(abun,aes(Year,log_abun,ymax=upper,ymin=lower,fill=Category))+
    geom_ribbon()+
    geom_line()+
    scale_fill_npg()+
    guides(fill='none')+
    facet_wrap(~Category,ncol=1,scales="free_y")+
    labs(x="Year",y="Log Abundance")

}

plot_cog <- function(Report,dat){
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  cog_arr <- Report$mean_Z_cym
  ncats <- dim(cog_arr)[1]
  nyrs <- dim(cog_arr)[2]
  years <- seq(min(dat$year),max(dat$year))
  dim(cog_arr) <- c(ncats*nyrs,dim(cog_arr)[3])
  
  cog_tbl <- tibble(category=rep(cats,nyrs),year=rep(years,each=ncats),East=cog_arr[,1],North=cog_arr[,2]) %>% 
    gather(direction,value,East,North)
  
  ggplot(cog_tbl,aes(year,value,col=category))+
    geom_line(size=2)+
    scale_color_npg()+
    guides(color='none')+
    facet_grid(direction~category,scales="free")+
    labs(x="Year",y="Kilometers")
}

plot_2d_cog <- function(Report,dat) {
  
  cats <- factor(tools::toTitleCase(levels(dat$spp)),levels=tools::toTitleCase(levels(dat$spp)))
  
  cog_arr <- Report$mean_Z_cym
  ncats <- dim(cog_arr)[1]
  nyrs <- dim(cog_arr)[2]
  years <- seq(min(dat$year),max(dat$year))
  dim(cog_arr) <- c(ncats*nyrs,dim(cog_arr)[3])
  
  cog_tbl <- tibble(category=rep(cats,nyrs),year=rep(years,each=ncats),east=cog_arr[,1],north=cog_arr[,2]) %>% 
    group_by(category) %>% 
    mutate(leadeast=lead(east,1),leadnorth=lead(north,1))
  out <- ggplot(cog_tbl,aes(east,north,col=year,xend=leadeast,yend=leadnorth))+
    geom_point()+
    geom_segment(arrow=arrow(angle=20,length = unit(0.2, "cm"),type='closed'))+
    labs(col='Year',x="Eastings",y="Northings")+
    coord_equal()+
    scale_color_viridis_c()+
    scale_size_continuous(range=c(1,4))+
    facet_wrap(~category)
  out
}

# Factor correlations with temperature

find_density_covariates <- function(Spatial_List,dat) {
  # locations of knots from Spatial_List object
  loc_query <- Spatial_List[['loc_x']]
  # years
  yrs <- seq(min(dat$year),max(dat$year))
  # locations of interpolated temperature and depth info
  loc_covars <- nbt.interpolated %>% distinct(x,y)
  # nearest neighbors (knots and interpolated covariates)
  which_nn <- nn2(data=loc_covars,query=loc_query,k=1)$nn.idx[,1]
  # for each knot in each year, find nearest neighbor from the interpolated datasets
  X_xtp <- array(data=NA,dim=c(nrow(loc_query),length(yrs),2))
  # fill in covariates for each knot in each year
  # depth is constant across years
  for(j in 1:length(yrs)) {
    yr_temp <- nbt.interpolated %>% filter(year==yrs[j])
    X_xtp[,j,1] <- yr_temp$temp[which_nn]
    X_xtp[,j,2] <- depth.interpolated$depth[which_nn]
  }
  return(X_xtp)
}
fct_covariate_correlation <- function(Spatial_List,dat,factors,which_lpred,which_f){
  covars <- find_density_covariates(Spatial_List,dat)
  knots <- 1:dim(covars)[1]
  years <- seq(min(dat$year),max(dat$year))
  dim(covars) <- c(dim(covars)[1]*dim(covars)[2],2)
  covars <- as.data.frame(covars) %>% set_names("temp","depth") %>% mutate(knot=rep(knots,length(years)),year=rep(years,each=length(knots)))

  fcts <- factors$Rotated_factors %>% pluck(which_lpred)
  n_knots <- dim(fcts)[1]
  fct <- fcts[,which_f,]
  n_knots<-ifelse(grepl('Omega',which_lpred),length(fct),dim(fct)[1])
  n_yrs <- ifelse(grepl('Omega',which_lpred),1,dim(fct)[2])
  
  if(grepl('Epsilon',which_lpred)) {
    # if looking at spatio-temporal correlations, use temperature anomalies
    covars_anom <- covars %>%
      select(knot,year,temp) %>% 
      group_by(knot) %>% 
      mutate(temp_anom=temp-mean(temp)) %>% 
      ungroup()
    df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct)) %>% 
      mutate(year=rep(Years2Include,each=n_knots)) %>% 
      left_join(covars_anom,by=c('knot','year')) %>% 
      filter(!is.na(temp_anom))
    out <- df %>% 
      group_by(lpred,which_f,year) %>% 
      summarize(temp_anom_corr=cor.sig(value,temp_anom))
  }
  if(grepl('Omega',which_lpred)) {
    df <- tibble(lpred=which_lpred,which_f=which_f,knot=rep(1:n_knots,n_yrs),value=as.numeric(fct)) %>% 
      left_join(covars,by=c('knot'))%>% 
      filter(!is.na(temp),!is.na(depth)) %>% 
      distinct()
    out <- df %>% 
      group_by(lpred,which_f) %>% 
      summarize(tempcorr=cor.sig(value,temp),depthcorr=cor.sig(value,depth))
  }
  out
}
