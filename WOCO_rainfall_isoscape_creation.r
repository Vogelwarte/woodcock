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
  mutate(source_ID=as.factor(SampleSiteName), year=year(SampleDate),month=month(SampleDate)) %>%
  mutate(year=ifelse(month>9,year+1,year)) %>%   ### adjust the year based on the woodcock moult cycle, which is complete in September, so birds with feathers from 2017 will have isotope ratios from Oct 2016 - Sept 2017
  rename(lat=Latitude, long=Longitude, source_value=dH, elev=Altitude) %>%
  dplyr::select(source_ID,lat,long,elev,year,month,source_value)
head(GNIPData)
hist(GNIPData$source_value)


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
EuropeFit03 <- isofit(data = GNIPDataEUagg[GNIPDataEUagg$year==2003,],
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
ElevEurope<- terra::rast('data/elevation_world_z5.tif')
#plot(ElevEurope)


ElevEurope02 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit02,
                           aggregation_factor = 4)

ElevEurope03 <- prepraster(raster = ElevEurope,
                           isofit = EuropeFit03,
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
## commented out as it takes a few hours (~15 min per isoscape)

# EuropeIsoscape02 <- isoscape(raster = ElevEurope02,
#                              isofit = EuropeFit02)
# 
# EuropeIsoscape03 <- isoscape(raster = ElevEurope03,
#                              isofit = EuropeFit03)
# 
# EuropeIsoscape05 <- isoscape(raster = ElevEurope05,
#                              isofit = EuropeFit05)
# 
# EuropeIsoscape07 <- isoscape(raster = ElevEurope07,
#                              isofit = EuropeFit07)
# 
# EuropeIsoscape08 <- isoscape(raster = ElevEurope08,
#                              isofit = EuropeFit08)
# 
# EuropeIsoscape09 <- isoscape(raster = ElevEurope09,
#                              isofit = EuropeFit09)
# 
# EuropeIsoscape10 <- isoscape(raster = ElevEurope10,
#                              isofit = EuropeFit10)
# 
# EuropeIsoscape13 <- isoscape(raster = ElevEurope13,
#                              isofit = EuropeFit13)
# 
# EuropeIsoscape14 <- isoscape(raster = ElevEurope14,
#                              isofit = EuropeFit14)
# 
# EuropeIsoscape15 <- isoscape(raster = ElevEurope15,
#                              isofit = EuropeFit15)
# 
# EuropeIsoscape16 <- isoscape(raster = ElevEurope16,
#                              isofit = EuropeFit16)
# 
# EuropeIsoscape17 <- isoscape(raster = ElevEurope17,
#                            isofit = EuropeFit17)
# 
# EuropeIsoscape18 <- isoscape(raster = ElevEurope18,
#                              isofit = EuropeFit18)


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
  bind_rows(ORIG_WC)
dim(ORIG_WC)
table(ORIG_WC$Year)
summary(ORIG_WC)

## 6.3. SPecify when feathers were grown -----------
## sampled feathers are the first secondary in shot birds and greater coverts in live birds
## according to Glutz von Blotzheim adults moult after June until Sept, and juveniles can moult part of their wing coverts and secondaries in Aug/Sept

ORIG_WC<-ORIG_WC %>% dplyr::select(-PROVENANCE_voigt,-ADULTE) %>%
  mutate(sd=parse_date_time(DATE, orders="dmy")) %>%
  mutate(feather_growth_date=if_else(AGE=="Adulte",
                                      if_else(month(sd)<7,
                                              ymd(paste(year(sd)-1,"-06-15")),
                                     ymd(paste(year(sd),"-06-15"))),
                                     if_else(month(sd)<9,
                                             ymd(paste(year(sd)-1,"-05-15")),
                                             ymd(paste(year(sd),"-08-15"))))) %>%
  mutate(feather_sampling_date=as.Date(sd)) %>%
  mutate(Year=year(feather_growth_date)) %>%
  dplyr::select(-DATE,-sd)


UNK_WC<-UNK_WC %>% 
  mutate(AGE=ifelse(is.na(AGE), "Juvénile",AGE)) %>%  ## fill in ages for 10 birds shot in Ticino in October - very likely juveniles
  mutate(feather_growth_date=if_else(AGE=="Adulte",
                                     if_else(month(DATE)<7,
                                             ymd(paste(year(DATE)-1,"-06-15")),
                                             ymd(paste(year(DATE),"-06-15"))),
                                     if_else(month(DATE)<9,
                                             ymd(paste(year(DATE)-1,"-05-15")),
                                             ymd(paste(year(DATE),"-08-15"))))) %>%
  mutate(feather_sampling_date=as.Date(DATE)) %>%
  mutate(Year=year(feather_growth_date)) %>%
  dplyr::select(-DATE)
    




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
# 
# globcover02 <- terra::resample(globcover, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
# dem02 <- terra::resample(dem, EuropeIsoscape02$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape02 <- (EuropeIsoscape02$isoscapes %>%
#   terra::mask(woco.countries))*globcover02*dem02
# rm(globcover02,dem02)
# gc()
# 
# globcover03 <- terra::resample(globcover, EuropeIsoscape03$isoscapes, method = "max")  # or method = "near" for categorical data
# dem03 <- terra::resample(dem, EuropeIsoscape03$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape03 <- (EuropeIsoscape03$isoscapes %>%
#                       terra::mask(woco.countries))*globcover03*dem03
# rm(globcover03,dem03)
# gc()
# 
# globcover05 <- terra::resample(globcover, EuropeIsoscape05$isoscapes, method = "max")  # or method = "near" for categorical data
# dem05 <- terra::resample(dem, EuropeIsoscape05$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape05 <- (EuropeIsoscape05$isoscapes %>%
#   terra::mask(woco.countries))*globcover05*dem05
# rm(globcover05,dem05)
# gc()
# 
# globcover07 <- terra::resample(globcover, EuropeIsoscape07$isoscapes, method = "max")  # or method = "near" for categorical data
# dem07 <- terra::resample(dem, EuropeIsoscape07$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape07 <- (EuropeIsoscape07$isoscapes %>%
#   terra::mask(woco.countries))*globcover07*dem07
# rm(globcover07,dem07)
# gc()
# 
# globcover08 <- terra::resample(globcover, EuropeIsoscape08$isoscapes, method = "max")  # or method = "near" for categorical data
# dem08 <- terra::resample(dem, EuropeIsoscape08$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape08 <- (EuropeIsoscape08$isoscapes %>%
#   terra::mask(woco.countries))*globcover08*dem08
# rm(globcover08,dem08)
# gc()
# 
# globcover09 <- terra::resample(globcover, EuropeIsoscape09$isoscapes, method = "max")  # or method = "near" for categorical data
# dem09 <- terra::resample(dem, EuropeIsoscape09$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape09 <- (EuropeIsoscape09$isoscapes %>%
#   terra::mask(woco.countries))*globcover09*dem09
# rm(globcover09,dem09)
# gc()
# 
# globcover10 <- terra::resample(globcover, EuropeIsoscape10$isoscapes, method = "max")  # or method = "near" for categorical data
# dem10 <- terra::resample(dem, EuropeIsoscape10$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape10 <- (EuropeIsoscape10$isoscapes %>%
#   terra::mask(woco.countries))*globcover10*dem10
# rm(globcover10,dem10)
# gc()
# 
# globcover13 <- terra::resample(globcover, EuropeIsoscape13$isoscapes, method = "max")  # or method = "near" for categorical data
# dem13 <- terra::resample(dem, EuropeIsoscape13$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape13 <- (EuropeIsoscape13$isoscapes %>%
#   terra::mask(woco.countries))*globcover13*dem13
# rm(globcover13,dem13)
# gc()
# 
# globcover14 <- terra::resample(globcover, EuropeIsoscape14$isoscapes, method = "max")  # or method = "near" for categorical data
# dem14 <- terra::resample(dem, EuropeIsoscape14$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape14 <- (EuropeIsoscape14$isoscapes %>%
#   terra::mask(woco.countries))*globcover14*dem14
# rm(globcover14,dem14)
# gc()
# 
# globcover15 <- terra::resample(globcover, EuropeIsoscape15$isoscapes, method = "max")  # or method = "near" for categorical data
# dem15 <- terra::resample(dem, EuropeIsoscape15$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape15 <- (EuropeIsoscape15$isoscapes %>%
#   terra::mask(woco.countries))*globcover15*dem15
# rm(globcover15,dem15)
# gc()
# 
# globcover16 <- terra::resample(globcover, EuropeIsoscape16$isoscapes, method = "max")  # or method = "near" for categorical data
# dem16 <- terra::resample(dem, EuropeIsoscape16$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape16 <- (EuropeIsoscape16$isoscapes %>%
#   terra::mask(woco.countries))*globcover16*dem16
# rm(globcover16,dem16)
# gc()
# 
# globcover17 <- terra::resample(globcover, EuropeIsoscape17$isoscapes, method = "max")  # or method = "near" for categorical data
# dem17 <- terra::resample(dem, EuropeIsoscape17$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape17 <- (EuropeIsoscape17$isoscapes %>%
#   terra::mask(woco.countries))*globcover17*dem17
# rm(globcover17,dem17)
# gc()
# 
# globcover18 <- terra::resample(globcover, EuropeIsoscape18$isoscapes, method = "max")  # or method = "near" for categorical data
# dem18 <- terra::resample(dem, EuropeIsoscape18$isoscapes, method = "max")  # or method = "near" for categorical data
# WOCO.isoscape18 <- (EuropeIsoscape18$isoscapes %>%
#   terra::mask(woco.countries))*globcover18*dem18
# rm(globcover18,dem18)
# gc()
# 
# 
# saveRDS(WOCO.isoscape02,"./data/isoscape02.rds")
# saveRDS(WOCO.isoscape03,"./data/isoscape03.rds")
# saveRDS(WOCO.isoscape05,"./data/isoscape05.rds")
# saveRDS(WOCO.isoscape07,"./data/isoscape07.rds")
# saveRDS(WOCO.isoscape08,"./data/isoscape08.rds")
# saveRDS(WOCO.isoscape09,"./data/isoscape09.rds")
# saveRDS(WOCO.isoscape10,"./data/isoscape10.rds")
# saveRDS(WOCO.isoscape13,"./data/isoscape13.rds")
# saveRDS(WOCO.isoscape14,"./data/isoscape14.rds")
# saveRDS(WOCO.isoscape15,"./data/isoscape15.rds")
# saveRDS(WOCO.isoscape16,"./data/isoscape16.rds")
# saveRDS(WOCO.isoscape17,"./data/isoscape17.rds")
# saveRDS(WOCO.isoscape18,"./data/isoscape18.rds")

WOCO.isoscape02<-readRDS("./data/isoscape02.rds")
WOCO.isoscape03<-readRDS("./data/isoscape03.rds")
WOCO.isoscape05<-readRDS("./data/isoscape05.rds")
WOCO.isoscape07<-readRDS("./data/isoscape07.rds")
WOCO.isoscape08<-readRDS("./data/isoscape08.rds")
WOCO.isoscape09<-readRDS("./data/isoscape09.rds")
WOCO.isoscape10<-readRDS("./data/isoscape10.rds")
WOCO.isoscape13<-readRDS("./data/isoscape13.rds")
WOCO.isoscape14<-readRDS("./data/isoscape14.rds")
WOCO.isoscape15<-readRDS("./data/isoscape15.rds")
WOCO.isoscape16<-readRDS("./data/isoscape16.rds")
WOCO.isoscape17<-readRDS("./data/isoscape17.rds")
WOCO.isoscape18<-readRDS("./data/isoscape18.rds")


isoscapes<-list(WOCO.isoscape02,WOCO.isoscape03,WOCO.isoscape05,WOCO.isoscape07,WOCO.isoscape08,WOCO.isoscape09,WOCO.isoscape10,
                WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17,WOCO.isoscape18)



## 6.6. extract hydrogen isotope values from that distribution for each year for calibration -------

### 6.6.1 convert known origin woco data to SpatVector to extract values from raster ------
woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.vect<-terra::vect(woco.sf)

yearlist<-c(2002,2003,2005,2007,2008,2009,2010,2013,2014,2015,2016,2017,2018)


for (y in sort(unique(woco.sf$Year))){
  
  woco.sf$d2h_MA[woco.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.vect[woco.vect$Year==y])$mean
  woco.sf$d2h_se_MA[woco.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.vect[woco.vect$Year==y])$mean_respVar
  
}



### 6.6.2 convert unknown origin woco data to SpatVector to extract LOCAL values from raster ------
woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)

for (y in sort(unique(woco.unk.sf$Year))){
  
  woco.unk.sf$d2h_MA[woco.unk.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.unk.vect[woco.unk.vect$Year==y])$mean
  woco.unk.sf$d2h_se_MA[woco.unk.sf$Year==y]<-terra::extract(isoscapes[[match(y,yearlist)]],woco.unk.vect[woco.unk.vect$Year==y])$mean_respVar
  
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

WOCO.isoscape13 <- WOCO.isoscape13 %>%
                      terra::mask(woco.countries.orig)
WOCO.isoscape14 <- (WOCO.isoscape14 %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape15 <- (WOCO.isoscape15 %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape16 <- (WOCO.isoscape16 %>%
                      terra::mask(woco.countries.orig ))
WOCO.isoscape17 <- (WOCO.isoscape17 %>%
                      terra::mask(woco.countries.orig ))


isoscapes.orig<-list(WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17)


## extract mean hydrogen isotope values and sd from that distribution

mean.rain.d2H.Europe<-as.numeric()
sd.rain.d2H.Europe<-as.numeric()
yearlist<-c(2013,2014,2015,2016,2017)


for (y in sort(unique(UNK_WC$Year))){
  
  rain.d2H<-as.numeric(na.omit(terra::values(isoscapes.orig[[match(y,yearlist)]])[,1]))
  rain.d2H<-rain.d2H[rain.d2H!=0]  ## remove the non-forest values (>10,0000 grid cells are removed)
  mean.rain.d2H.Europe[match(y,yearlist)]<-mean(rain.d2H, na.rm=T)
  
  ## calculating pooled variance as sum of sample and other variance: https://stats.stackexchange.com/questions/604824/how-can-i-average-the-variance-extracted-from-different-group-of-samples
  # Sample variance of predictions:
  sample_var.Europe <- as.numeric(na.omit(terra::values(isoscapes.orig[[match(y,yearlist)]])[,5]))
  sample_var.Europe <-sample_var.Europe[sample_var.Europe!=0]  ## remove the non-forest values (>10,0000 grid cells are removed)
  
  # Mean of individual variances:
  mean_var.Europe <- mean(sample_var.Europe, na.rm=T)
  
  # Total variance:
  total_var.Europe <- var(rain.d2H, na.rm=T) + mean_var.Europe
  sd.rain.d2H.Europe[match(y,yearlist)]<-sqrt(total_var.Europe)
  
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





