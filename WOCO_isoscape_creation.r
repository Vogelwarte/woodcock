#---------------------------------------------------------------------------
### CREATION OF HYDROGEN ISOSCAPE ACROSS EUROPE ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 7 Nov 2025

## following tutorial in: https://bookdown.org/content/782/isoscape.html

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



# 1. READ IN GNIP ISOTOPE DATA ----------------

## manually downloaded on 7 Nov 2025 from https://nucleus.iaea.org/wiser/explore/

GNIPData <- read_excel("./data/GNIP_water_isotopes_download.xlsx") %>%   ## changed from GNIP_Wiser_monthly_hydrogen.xlsxafter downloading longer timeseries 
  dplyr::filter(Measurand=="H2") %>%  ## added this because I also downloaded rainfall amount
  dplyr::select(SampleSiteName,Latitude, Longitude,Altitude,SampleDate,dH) %>%
  dplyr::filter(!is.na(dH)) %>%
  dplyr::filter(!is.na(dH)) %>%
  mutate(source_ID=as.factor(SampleSiteName), year=year(SampleDate),month=month(SampleDate)) %>%
  mutate(year=ifelse(month>9,year+1,year)) %>%   ### adjust the year based on the woodcock moult cycle, which is complete in September, so birds with feathers from 2017 will have isotope ratios from Oct 2016 - Sept 2017
  rename(lat=Latitude, long=Longitude, source_value=dH, elev=Altitude) %>%
  dplyr::select(source_ID,lat,long,elev,year,month,source_value)
head(GNIPData)



# GNIPDataEUagg_test <- prepsources(data = GNIPData,
#                                   long_min = -30, long_max = 100,
#                                   lat_min = 30, lat_max = 90)


## 1.1. aggregating data to 1 value per station and year---------
GNIPDataEUagg <- GNIPData %>%
  group_by(source_ID,lat,long,elev,year) %>%
  summarise(mean_source_value=mean(source_value),var_source_value=var(source_value),n_source_value=length(source_value))




# 2. FITTING GEOSTATISTICAL MODELS ----------------


EuropeFit02 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2002,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit05 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2005,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit07 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2007,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit08 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2008,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit09 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2009,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit10 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2010,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))


EuropeFit13 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2013,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit14 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2014,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit15 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2015,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit16 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2016,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit17 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2017,],
                    mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
EuropeFit18 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2018,],
                      mean_model_fix = list(elev = TRUE, lat_abs = TRUE))




# 3. READING IN DEM ----------------
## done once on 10 Nov 2025

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
ElevEurope<- terra::rast('C:/Users/sop/OneDrive - Vogelwarte/Dokumente/elevation_world_z5.tif')
#plot(ElevEurope)


ElevEurope02 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit02,
                           aggregation_factor = 4)

ElevEurope05 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit05,
                           aggregation_factor = 4)

ElevEurope07 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit07,
                           aggregation_factor = 4)

ElevEurope08 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit08,
                           aggregation_factor = 4)

ElevEurope09 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit09,
                           aggregation_factor = 4)

ElevEurope10 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit10,
                           aggregation_factor = 4)



ElevEurope13 <- prepraster(raster = ElevEurope,
                         isofit = EuropeFit13,
                         aggregation_factor = 4)

ElevEurope14 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit14,
                           aggregation_factor = 4)

ElevEurope15 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit15,
                           aggregation_factor = 4)
ElevEurope16 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit16,
                           aggregation_factor = 4)

ElevEurope17 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit17,
                           aggregation_factor = 4)

ElevEurope18 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit18,
                           aggregation_factor = 4)



# 4. BUILDING ANNUAL ISOSCAPES ----------------


EuropeIsoscape02 <- isoscape(raster = ElevEurope02,
                             isofit = EuropeFit02)

EuropeIsoscape05 <- isoscape(raster = ElevEurope05,
                             isofit = EuropeFit05)

EuropeIsoscape07 <- isoscape(raster = ElevEurope07,
                             isofit = EuropeFit07)

EuropeIsoscape08 <- isoscape(raster = ElevEurope08,
                             isofit = EuropeFit08)

EuropeIsoscape09 <- isoscape(raster = ElevEurope09,
                             isofit = EuropeFit09)

EuropeIsoscape10 <- isoscape(raster = ElevEurope10,
                             isofit = EuropeFit10)



EuropeIsoscape13 <- isoscape(raster = ElevEurope13,
                             isofit = EuropeFit13)

EuropeIsoscape14 <- isoscape(raster = ElevEurope14,
                             isofit = EuropeFit14)

EuropeIsoscape15 <- isoscape(raster = ElevEurope15,
                             isofit = EuropeFit15)

EuropeIsoscape16 <- isoscape(raster = ElevEurope16,
                             isofit = EuropeFit16)

EuropeIsoscape17 <- isoscape(raster = ElevEurope17,
                           isofit = EuropeFit17)

EuropeIsoscape18 <- isoscape(raster = ElevEurope18,
                             isofit = EuropeFit18)


# 5. SAVING / RELOADING ANNUAL ISOSCAPES ----------------
## THIS CAUSES A CATASTROPHIC DATA LOSS: https://github.com/rspatial/terra/issues/549
# save.image("./data/isoscapes.RData")
# load("./data/isoscapes.RData")

# saveRDS(EuropeIsoscape02,"./data/isoscape02.rds")
# saveRDS(EuropeIsoscape05,"./data/isoscape05.rds")
# saveRDS(EuropeIsoscape07,"./data/isoscape07.rds")
# saveRDS(EuropeIsoscape08,"./data/isoscape08.rds")
# saveRDS(EuropeIsoscape09,"./data/isoscape09.rds")
# saveRDS(EuropeIsoscape10,"./data/isoscape10.rds")
# saveRDS(EuropeIsoscape13,"./data/isoscape13.rds")
# saveRDS(EuropeIsoscape14,"./data/isoscape14.rds")
# saveRDS(EuropeIsoscape15,"./data/isoscape15.rds")
# saveRDS(EuropeIsoscape16,"./data/isoscape16.rds")
# saveRDS(EuropeIsoscape17,"./data/isoscape17.rds")
# saveRDS(EuropeIsoscape18,"./data/isoscape18.rds")
# 



# 6. EXTRACTING ISOTOPE DISTRIBUTIONS FROM SAMPLED FEATHERS ---------------------

## 6.1. LOADING RAW ISOTOPE DATA ----------------

woco<-fread("data/WOCO_isotopes.csv")
table(woco$JAHR)


## back-convert all measurements to INTERNATIONAL standards following Soto et al 2017
woco$dH<-10.774 + (0.852*woco$dH_scaled)
woco %>% dplyr::filter(is.na(dH))





## 6.2. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

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

### 6.2.1. add data of known origin provided by andrew Hoodless
ORIG_WC<-fread("data/WOCO_known_origin_feathers_Hoodless.csv") %>%
  mutate(ADULTE=ifelse(Age=="Adult",1,0),PROVENANCE_voigt=as.character(LocationCode),
         AGE=ifelse(Age=="Adult", "Adulte","Juvénile")) %>%
  mutate(dH=10.774 + (0.852*DHF)) %>%   ## new correction from Soto et al 2017 inserted after discussion with David Soto on 20 Oct 2025
  rename(ID=RefID, KANTON=Region,
         DATE=Date,LATITUDE=Latitude,LONGITUDE=Longitude) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
  bind_rows(ORIG_WC) %>%
  mutate(Year=year(dmy(DATE)))
dim(ORIG_WC)
table(ORIG_WC$Year)


## 6.3. SPecify when feathers were grown -----------
## sampled feathers are the first secondary in shot birds and greater coverts in live birds
## according to Glutz von Blotzheim adults moult after June until Sept, and juveniles can moult part of their wing coverts and secondaries in Aug/Sept

rain_orig_wc<-ORIG_WC %>% dplyr::select(-PROVENANCE_voigt,-ADULTE) %>%
  mutate(sd=parse_date_time(DATE, orders="dmy")) %>%
  mutate(feather_growth_date=if_else(AGE=="Adulte",
                                      if_else(month(sd)<7,
                                              ymd(paste(year(sd)-1,"-06-15")),
                                     ymd(paste(year(sd),"-06-15"))),
                                     if_else(month(sd)<9,
                                             ymd(paste(year(sd)-1,"-05-15")),
                                             ymd(paste(year(sd),"-08-15"))))) %>%
  mutate(feather_sampling_date=as.Date(sd)) %>%
  dplyr::select(-DATE,-sd,-AGE)
    




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



### 6.4.1. curtail woodcock distribution to potential origin countries OUTSIDE OF SWITZERLAND and FOREST
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

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
  crop(woco.countries) %>%
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





### 6.4.2. multiply isoscape with woodcock distribution and generate mean feather distribution
# WOCO.isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
#   crop(woco.countries)
# plot(WOCO.isoscape)


## check whether crs is the same
crs(globcover)==crs(EuropeIsoscape02$isoscapes)
origin(globcover)
origin(EuropeIsoscape02$isoscapes)
origin(dem)

## align extent
globcover <- terra::resample(globcover, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
dem <- terra::resample(ElevEurope, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data


ext(EuropeIsoscape02$isoscapes)
ext(globcover)
ext(dem)


res(EuropeIsoscape02$isoscapes)
res(globcover)
res(dem)

compareGeom(EuropeIsoscape02$isoscapes, globcover, dem, stopOnError = FALSE)


## REMOVE NON-FOREST AND HIGH ELEVATION FROM EACH ISOSCAPE
# DO NOT USE terra::crop here because it will change the extent of the raster and the multiplication will fail
WOCO.isoscape02 <- (EuropeIsoscape02$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape05 <- (EuropeIsoscape05$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape07 <- (EuropeIsoscape07$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape08 <- (EuropeIsoscape08$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape09 <- (EuropeIsoscape09$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape10 <- (EuropeIsoscape10$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape13 <- (EuropeIsoscape13$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape14 <- (EuropeIsoscape14$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape15 <- (EuropeIsoscape15$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape16 <- (EuropeIsoscape16$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape17 <- (EuropeIsoscape17$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem

WOCO.isoscape18 <- (EuropeIsoscape18$isoscapes %>%
  terra::mask(woco.countries))*globcover*dem


saveRDS(EuropeIsoscape02,"./data/isoscape02.rds")
saveRDS(EuropeIsoscape05,"./data/isoscape05.rds")
saveRDS(EuropeIsoscape07,"./data/isoscape07.rds")
saveRDS(EuropeIsoscape08,"./data/isoscape08.rds")
saveRDS(EuropeIsoscape09,"./data/isoscape09.rds")
saveRDS(EuropeIsoscape10,"./data/isoscape10.rds")
saveRDS(EuropeIsoscape13,"./data/isoscape13.rds")
saveRDS(EuropeIsoscape14,"./data/isoscape14.rds")
saveRDS(EuropeIsoscape15,"./data/isoscape15.rds")
saveRDS(EuropeIsoscape16,"./data/isoscape16.rds")
saveRDS(EuropeIsoscape17,"./data/isoscape17.rds")
saveRDS(EuropeIsoscape18,"./data/isoscape18.rds")




## 6.5. extract hydrogen isotope values from that distribution for each year -------


## extract hydrogen isotope values from that distribution
rain.d2H<-as.numeric(na.omit(terra::values(WOCO.isoscape02[[1]])[,1]))
rain.d2H<-rain.d2H[rain.d2H<0]  ## remove the non-forest values (>10,0000 grid cells are removed)
mean.rain.d2H<-mean(rain.d2H, na.rm=T)
sd.rain.d2H<-sd(rain.d2H, na.rm=T)
hist(rnorm(1000,mean.rain.d2H,sd.rain.d2H))




## 6.6. use known origin data to calibrate rainwater against feather d2H  -------


### 6.6.1 convert woco data to SpatVector to extract values from raster ------
woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.vect<-terra::vect(woco.sf)


## 6.6.2. extract rainwater hydrogen isotopes for the location of known-sample woodcocks -------

woco.sf$d2h_MA<-terra::extract(WOCO.isoscape02,woco.vect)$mean
woco.sf$d2h_se_MA<-terra::extract(isoscape,woco.vect)$mean_predVar

