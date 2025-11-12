#---------------------------------------------------------------------------
### CREATION OF WOODCOCK FEATHER HYDROGEN ISOSCAPE ACROSS EUROPE ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 7 Nov 2025

## following tutorial in: https://bookdown.org/content/782/isoscape.html

## modified based on suggestion by David Soto to create a direct feather isoscape rather than relate to rainfall

rm(list=ls())
library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
library(readxl)
library(terra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(IsoriX)
filter<-dplyr::filter
select<-dplyr::select
library(gridExtra)
library(tidyterra)


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock"),silent=T)



# 1. READ IN WOODCOCK ISOTOPE DATA ----------------

## 1.1. LOADING RAW ISOTOPE DATA ----------------

woco<-fread("data/WOCO_isotopes.csv")
table(woco$JAHR)


## back-convert all measurements to INTERNATIONAL standards following Soto et al 2017
woco$dH<-10.774 + (0.852*woco$dH_scaled)
woco %>% dplyr::filter(is.na(dH))





## 1.2. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
  dplyr::filter(!is.na(dH)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE)
UNK_WC<-woco %>% dplyr::filter(ORIGINE=="UNBEKANNT") %>%
  dplyr::filter(JAGDT==1) %>%
  dplyr::filter(!is.na(DATE)) %>%
  dplyr::filter(!is.na(dH)) %>%
  mutate(Date=lubridate::parse_date_time(x=DATE,orders="dmy", tz = "UTC", drop=T)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,Date,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
  rename(DATE=Date)
UNK_WC$ID<-str_replace_all(UNK_WC$ID, "[^[:alnum:]]", " ")
dim(ORIG_WC)
dim(UNK_WC)

### 1.2.1. add data of known origin provided by andrew Hoodless-------
ORIG_WC<-fread("data/WOCO_known_origin_feathers_Hoodless.csv") %>%
  mutate(ADULTE=ifelse(Age=="Adult",1,0),PROVENANCE_voigt=as.character(LocationCode),
         AGE=ifelse(Age=="Adult", "Adulte","Juvénile")) %>%
  mutate(dH=10.774 + (0.852*DHF)) %>%   ## new correction from Soto et al 2017 inserted after discussion with David Soto on 20 Oct 2025
  rename(ID=RefID, KANTON=Region,
         DATE=Date,LATITUDE=Latitude,LONGITUDE=Longitude) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,LATITUDE,LONGITUDE) %>%
  bind_rows(ORIG_WC)
dim(ORIG_WC)
summary(ORIG_WC)




## 1.3. ADD ELEVATION DATA FOR ISOTOPE RATIOS ----


### 1.3.1. load shapefile of Switzerland and EUROPE -------

SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84

bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)



### 1.3.2. specify origin countries INCLUDING UK AND SWITZERLAND FOR CALIBRATION --------
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium", "Switzerland","France","United Kingdom","Spain","Italy"))



### 1.3.3. load DEM across Europe

#ElevEurope<-terra::rast("data/eurodem.tif")  ## manually downloaded from Copernicus
# ElevEurope<-getelev(
#   file = "~/elevation_world_z5.tif",
#   z = 5,
#   long_min = -10,
#   long_max = 80,
#   lat_min = 35,
#   lat_max = 80,
#   margin_pct = 5,
#   overwrite = TRUE
# )
ElevEurope<- terra::rast('data/elevation_world_z5.tif')
#plot(ElevEurope)



### 1.3.4. create sf object and extract elevation data from DEM --------

woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.vect<-terra::vect(woco.sf)
ORIG_WC$elev<-terra::extract(ElevEurope,woco.vect)$elevation_world_z5

woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)
UNK_WC$elev<-terra::extract(ElevEurope,woco.unk.vect)$elevation_world_z5














# 2. FITTING GEOSTATISTICAL MODELS ----------------


## 2.1. Formatting the input data as they are needed -------
ORIG_WC<-ORIG_WC %>%
  rename(source_ID=ID,lat=LATITUDE, long=LONGITUDE, mean_source_value=dH, age=ADULTE) %>%
  dplyr::select(source_ID,lat,long,elev,age,mean_source_value)
head(ORIG_WC)
hist(ORIG_WC$mean_source_value)

## calculating variance
vardH<-ORIG_WC %>% group_by(lat,long,elev) %>%
  summarise(var_source_value=var(mean_source_value), n_source_value=length(mean_source_value)) %>%
  arrange(desc(var_source_value))

# adding variance to dataset
ORIG_WC<-ORIG_WC %>%
  left_join(vardH, by=c('lat','long','elev')) %>%
  mutate(var_source_value=ifelse(is.na(var_source_value),max(vardH$var_source_value,na.rm=T),var_source_value)) %>%
  mutate(n_source_value=ifelse(n_source_value<5,5,n_source_value)) %>%
  ungroup() %>%
  group_by(source_ID, lat, long, elev)

# ORIG_WC_agg<-ORIG_WC %>%
#   group_by(lat,long,elev) %>%
#   summarise(source_ID=first(source_ID), mean_source_value=mean(source_value),var_source_value=var(source_value), n_source_value=length(source_value)) %>%
#   mutate(var_source_value=ifelse(is.na(var_source_value),max(vardH$var_source_value,na.rm=T),var_source_value)) %>%
#   mutate(n_source_value=ifelse(n_source_value<5,5,n_source_value))


which(is.na(ORIG_WC)==TRUE)



## 2.2. fitting geostatistical model PER AGE CLASS -----------



EuropeFitAD <- isofit(data = ORIG_WC[ORIG_WC$age==1,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE, long_abs=TRUE))

EuropeFitJuv <- isofit(data = ORIG_WC[ORIG_WC$age==0,],
                    mean_model_fix = list(elev = TRUE, lat_abs = TRUE, long_abs=TRUE))




# 3. ADJUSTING THE DEM ----------------

ElevEuropeAD <- prepraster(raster = ElevEurope,
                           isofit = EuropeFitAD,
                           aggregation_factor = 4)

ElevEuropeJuv <- prepraster(raster = ElevEurope,
                         isofit = EuropeFitJuv,
                         aggregation_factor = 4)



# 4. BUILDING WOODCOCK FEATHER ISOSCAPES ----------------


EuropeIsoscapeAD <- isoscape(raster = ElevEuropeAD,
                             isofit = EuropeFitAD)


EuropeIsoscapeJuv <- isoscape(raster = ElevEuropeJuv,
                           isofit = EuropeFitJuv)

plot(EuropeIsoscapeAD)
plot(EuropeIsoscapeJuv)


# 6. EXTRACTING ISOTOPE DISTRIBUTIONS FROM SAMPLED FEATHERS ---------------------

woco.vect<-terra::vect(woco.sf)

woco.sf$d2h_predicted<-terra::extract(EuropeIsoscape$isoscapes,woco.vect)$mean

plot(woco.sf$dH~woco.sf$d2h_predicted)


woco.unk.sf$d2h_local_predicted<-terra::extract(EuropeIsoscape$isoscapes,woco.unk.vect)$mean
woco.unk.sf$d2h_local_sd<-sqrt(terra::extract(EuropeIsoscape$isoscapes,woco.unk.vect)$mean_residVar)
plot(woco.unk.sf$dH~woco.unk.sf$d2h_local_predicted)



## 6.4. load shapefile of Switzerland and EUROPE -------
## OPTION TO CURTAIL TO FOREST AREAS - but only needed for blind assignment

SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84
plot(SUI)

bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)
plot(EUR)



### 6.4.1. curtail woodcock distribution to potential origin countries INCLUDING UK AND SWITZERLAND FOR CALIBRATION
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium", "Switzerland","France","United Kingdom","Spain","Italy"))




## 6.5. curtail woodcock distribution to forest and elevation <2000 m ---------------

## create conversion matrices
forest.mat<-matrix(0, nrow=3, ncol=3)
forest.mat[,1]<-c(11,40,110)  ## from values for conversion matrix
forest.mat[,2]<-c(30,100,230)  ## to values for conversion matrix
forest.mat[,3]<-c(0,1,0)  ## replacement values for conversion matrix

ele.mat<-matrix(0, nrow=3, ncol=3)
ele.mat[,1]<-c(-8000,0,2000)  ## from values for conversion matrix
ele.mat[,2]<-c(0,2000,8000)  ## to values for conversion matrix
ele.mat[,3]<-c(0,1,0)  ## replacement values for conversion matrix


## create elevation raster layer - unnecessary to read in new because ElevEurope already created above
# dem0<-terra::rast("data/eurodem.tif")
dem<-ElevEurope %>%
  terra::project(.,crs(woco.countries)) %>%
  terra::crop(woco.countries) %>%
  terra::classify(rcl=ele.mat,include.lowest=T,right=NA) %>%
  terra::project(.,crs(EuropeIsoscape02$isoscapes))
crs(dem)
summary(dem)
plot(dem)

## create forest raster layer
#globcover<-terra::rast("S:/rasters/landuse/world/globcover2009.tif") %>%
globcover<-terra::rast("data/GLOBCOVER_L4_200901_200912_V2.3.tif") %>%
  terra::crop(woco.countries) %>%
  terra::classify(rcl=forest.mat,include.lowest=T,right=NA) %>%
  terra::project(.,crs(EuropeIsoscape02$isoscapes))
crs(globcover)
summary(globcover)
plot(globcover)


## check whether crs is the same
# crs(globcover)==crs(EuropeIsoscape02$isoscapes)
# origin(globcover)
# origin(EuropeIsoscape02$isoscapes)
# origin(dem)
# 
# ## align extent
# globcover <- terra::resample(globcover, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
# dem <- terra::resample(ElevEurope, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
# 
# 
# ext(EuropeIsoscape05$isoscapes)
# ext(globcover)
# ext(dem)
# 
# 
# res(EuropeIsoscape05$isoscapes)
# res(globcover)
# res(dem)
# 
# compareGeom(EuropeIsoscape02$isoscapes, globcover, dem, stopOnError = FALSE)


## REMOVE NON-FOREST AND HIGH ELEVATION FROM EACH ISOSCAPE
# DO NOT USE terra::crop here because it will change the extent of the raster and the multiplication will fail

globcover02 <- terra::resample(globcover, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
dem02 <- terra::resample(dem, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
WOCO.isoscape02 <- (EuropeIsoscape02$isoscapes %>%
  terra::mask(woco.countries))*globcover02*dem02
rm(globcover02,dem02)
gc()




saveRDS(WOCO.isoscape02,"./data/isoscape02.rds")



## 6.6. extract hydrogen isotope values from that distribution for each year for calibration -------

### 6.6.1 convert known origin woco data to SpatVector to extract values from raster ------
woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.vect<-terra::vect(woco.sf)

yearlist<-c(2002,2003,2005,2007,2008,2009,2010,2013,2014,2015,2016,2017,2018)


for (y in sort(unique(woco.sf$Year))){
  
  woco.sf$d2h_MA[woco.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.vect[woco.vect$Year==y])$mean
  woco.sf$d2h_se_MA[woco.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.vect[woco.vect$Year==y])$mean_predVar
  
}



### 6.6.2 convert unknown origin woco data to SpatVector to extract LOCAL values from raster ------
woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)

for (y in sort(unique(woco.unk.sf$Year))){
  
  woco.unk.sf$d2h_MA[woco.unk.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.unk.vect[woco.unk.vect$Year==y])$mean
  woco.unk.sf$d2h_se_MA[woco.unk.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.unk.vect[woco.unk.vect$Year==y])$mean_predVar
  
}

## check elevation distribution of shot birds
summary(terra::extract(ElevEurope,woco.unk.vect))
hist(terra::extract(ElevEurope,woco.unk.vect)$elevation_world_z5)
length(which(terra::extract(ElevEurope,woco.unk.vect)$elevation_world_z5>2200))


## 6.7. extract hydrogen isotope values from that distribution for each year for calibration -------

### need to curtail countries to potential origin (not including calibration countries like SUI, UK, F, ESP)
woco.countries.orig <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

WOCO.isoscape13 <- (EuropeIsoscape13$isoscapes %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape14 <- (EuropeIsoscape14$isoscapes %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape15 <- (EuropeIsoscape15$isoscapes %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape16 <- (EuropeIsoscape16$isoscapes %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape17 <- (EuropeIsoscape17$isoscapes %>%
                      terra::mask(woco.countries.orig ))


isoscapes.orig<-list(WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17)


## extract mean hydrogen isotope values and sd from that distribution

mean.rain.d2H<-as.numeric()
sd.rain.d2H<-as.numeric()
yearlist<-c(2013,2014,2015,2016,2017)


for (y in sort(unique(UNK_WC$Year))){
  
  rain.d2H<-as.numeric(na.omit(terra::values(isoscapes.orig[[match(y,yearlist)]])[,1]))
  rain.d2H<-rain.d2H[rain.d2H!=0]  ## remove the non-forest values (>10,0000 grid cells are removed)
  mean.rain.d2H[match(y,yearlist)]<-mean(rain.d2H, na.rm=T)
  sd.rain.d2H[match(y,yearlist)]<-sd(rain.d2H, na.rm=T)
  
}





# 7. include phenology data to integrate probability of non-local origin based on date a bird was shot ----------

woco_mig<-readRDS("output/woco_mig_depart_simulation.rds") %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>%
  mutate(prop.remain=1-prop_mig) %>% 
  group_by(week, Date) %>%
  summarise(mig=quantile(prop.remain,0.5),mig.lcl=quantile(prop.remain,0.025),mig.ucl=quantile(prop.remain,0.975))

woco_abundance<-fread("data/AnnualPhenology_abundance_index.csv") %>%
  mutate(Date=lubridate::ymd("2023-12-27") + lubridate::days(Pentade*5)) %>%
  mutate(abund=SOPM/max(SOPM)) %>%
  dplyr::filter(Date>=min(woco_mig$Date))

shot_dates<-hist(lubridate::yday(UNK_WC$DATE), breaks=seq(250,365,7),plot=F)
woco_shot<-tibble(yday=shot_dates$mids, N=shot_dates$counts) %>%
  mutate(Date=parse_date_time(as.integer(yday), orders="j")) %>%
  mutate(abund=N/max(N)) %>%
  mutate(Date=as.Date(Date-years(1)))



# 8. clean up and save workspace -------


## 8.1. prepare the data needed for NIMBLE input ----
woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  #filter(KANTON !="VS") %>% ## remove Valais because only 6 birds from 1 age class, causes imbalance in data
  filter(!is.na(d2h_MA))

woco.sf <- woco.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_MA))

table(woco.unk.sf$AGE,woco.unk.sf$KANTON)


rm(isoscapes, globcover,forest.mat, isoscapes.orig,WOCO.isoscape02,WOCO.isoscape03,WOCO.isoscape05,WOCO.isoscape07,WOCO.isoscape08,WOCO.isoscape09,WOCO.isoscape10,
   WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17,WOCO.isoscape18,ElevEurope,ElevEurope02,ElevEurope03,ElevEurope05,ElevEurope07,ElevEurope08,ElevEurope09,ElevEurope10,ElevEurope13,ElevEurope14,ElevEurope15,ElevEurope16,ElevEurope17,ElevEurope18,      
EuropeFit02,EuropeFit03,EuropeFit05,EuropeFit07,EuropeFit08,EuropeFit09,EuropeFit10,EuropeFit13,EuropeFit14,EuropeFit15,EuropeFit16,EuropeFit17,EuropeFit18,EuropeIsoscape02,EuropeIsoscape03,EuropeIsoscape05,EuropeIsoscape07,EuropeIsoscape08,EuropeIsoscape09,
EuropeIsoscape10,EuropeIsoscape13,EuropeIsoscape14,EuropeIsoscape15,EuropeIsoscape16,EuropeIsoscape17,EuropeIsoscape18)
gc()

save.image("data/woco.input.data.annual.RData")




## 8.2. remove non-SUI calibration data and create reduced input data ----

woco.sf <- woco.sf %>%
  filter(nchar(KANTON)==2) 
dim(woco.sf)
save.image("data/woco.reduced.input.data.annual.RData")





