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

GNIPData <- read_excel("./data/GNIP_Wiser_monthly_hydrogen.xlsx") %>%
  dplyr::select(SampleSiteName,Latitude, Longitude,Altitude,SampleDate,dH) %>%
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


#ElevEurope<-terra::rast("data/eurodem.tif")  ## manually downloaded from Copernicus
# ElevEurope<-getelev(
#   file = "~/elevation_world_z5.tif",
#   z = 5,
#   long_min = -40,
#   long_max = 120,
#   lat_min = 20,
#   lat_max = 90,
#   margin_pct = 5
# )
ElevEurope<- terra::rast('C:/Users/sop/OneDrive - Vogelwarte/Dokumente/elevation_world_z5.tif')
#plot(ElevEurope)


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




# 6. LOADING RAW ISOTOPE DATA ----------------


woco<-fread("data/WOCO_isotopes.csv")
#woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")
table(woco$JAHR)


## back-convert all measurements to INTERNATIONAL standards following Soto et al 2017
# woco$dH<-(woco$dH_scaled - 20.701)/0.979 ## discontinued after discussion with David Soto on 20 Oct 2025
woco$dH<-10.774 + (0.852*woco$dH_scaled)
woco %>% dplyr::filter(is.na(dH))




## 1.3. PLOT HISTOGRAMS FOR SUISSE AND OTHER BIRDS ----

ggplot(woco, aes(x=dH, col=ORIGINE, fill=ORIGINE)) +
  geom_histogram(alpha=0.5,position = position_dodge(width=1))
#ggsave("output/WOCO_isotope_histogram_by_origin.jpg")


# distd2H_knownSUI<-hist(ORIG_WC$dH, breaks=seq(-130,30,5))
# distd2H_knownSUI$density   ## this cannot be used as prior information for the isotope assignment because it is just from a corner of Switzerland




## 1.4. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

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

### 1.4.1. add data of known origin provided by andrew Hoodless
ORIG_WC<-fread("data/WOCO_known_origin_feathers_Hoodless.csv") %>%
  mutate(ADULTE=ifelse(Age=="Adult",1,0),PROVENANCE_voigt=as.character(LocationCode),
         AGE=ifelse(Age=="Adult", "Adulte","Juvénile")) %>%
  mutate(dH=10.774 + (0.852*DHF)) %>%   ## new correction from Soto et al 2017 inserted after discussion with David Soto on 20 Oct 2025
  rename(ID=RefID, KANTON=Region,
         DATE=Date,LATITUDE=Latitude,LONGITUDE=Longitude) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
  bind_rows(ORIG_WC)
dim(ORIG_WC)


## 1.5. EXTRACT GROWTH SEASON RAIN ISOTOPES ----
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
    
#fwrite(rain_orig_wc,"data/woodcock_known_origin_samples.csv")



## report numbers in manuscript
table(ORIG_WC$AGE)
length(UNK_WC$dH)
table(UNK_WC$AGE)
summary(ORIG_WC$dH)
summary(UNK_WC$dH)



# 2. EXTRACTING ISOTOPE DISTRIBUTIONS FROM FOREST MAP OF EUROPE ---------------------

## 2.1. load shapefile of Switzerland and EUROPE -------
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


## 2.2. plot the known origin points of WOCO that were sampled -------

woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)


ggplot(EUR) +
  geom_sf() +
  geom_sf(data=woco.sf,color="red")





## 2.3. extract isoscape map of rainwater dH isotope ratios from across the world ----------

### this function throws a weird error - downloaded and ran function independently without problem
# isoscape <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200) %>%   ## we use MA because that is what Powell's conversion is based on
#   projectRaster(crs = CRS(SRS_string = 'EPSG:4326'))
## go to "C:/STEFFEN/Vogelwarte/assignR" and open assignR.Rproj
## then run the script "DOWNLOAD_ISOSCAPE.R"

### revised on 21 Oct 2025 to adopt David Soto's suggestion that annual mean is better to characterise food web water
# isoscape <- getIsoscapes(isoType = "GlobalPrecipMA", timeout = 1200) %>%   ## we use MA because that is what Powell's conversion is based on
#   projectRaster(crs = CRS(SRS_string = 'EPSG:4326'))

isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
  crop(extent(EUR))
plot(isoscape)
## downloaded from Nelson et al 2021, but poor resolution
# isoscape <- raster("data/Piso.AI_v1.2020_0.5deg_1950-2020.nc") %>%
#   crop(extent(SUI))



### 2.3.1. curtail woodcock distribution to potential origin countries OUTSIDE OF SWITZERLAND and FOREST
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

## create conversion matrices
forest.mat<-matrix(0, nrow=3, ncol=3)
forest.mat[,1]<-c(11,40,110)  ## from values for conversion matrix
forest.mat[,2]<-c(30,100,230)  ## to values for conversion matrix
forest.mat[,3]<-c(0,1,0)  ## replacement values for conversion matrix

ele.mat<-matrix(0, nrow=2, ncol=3)
ele.mat[,1]<-c(-500,2000)  ## from values for conversion matrix
ele.mat[,2]<-c(2000,5000)  ## to values for conversion matrix
ele.mat[,3]<-c(1,0)  ## replacement values for conversion matrix


## create elevation raster layer
dem0<-terra::rast("data/eurodem.tif")
dem<-terra::rast("data/eurodem.tif") %>%
  terra::project(.,crs(woco.countries)) %>%
  crop(woco.countries) %>%
  terra::classify(rcl=ele.mat,include.lowest=T,right=NA) %>%
  terra::project(.,crs(isoscape))
crs(dem)
summary(dem)
plot(dem)

## create forest raster layer
#globcover<-terra::rast("S:/rasters/landuse/world/globcover2009.tif") %>%
globcover<-terra::rast("data/GLOBCOVER_L4_200901_200912_V2.3.tif") %>%
  crop(woco.countries) %>%
  terra::classify(rcl=forest.mat,include.lowest=T,right=NA) %>%
  terra::project(.,crs(isoscape))
crs(globcover)
summary(globcover)
plot(globcover)





### 2.3.2. multiply isoscape with woodcock distribution and generate mean feather distribution
WOCO.isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
  crop(woco.countries)
plot(WOCO.isoscape)

## remove all non-forest areas by multiplying with 0
crs(globcover)==crs(WOCO.isoscape)
origin(globcover)
origin(WOCO.isoscape)

## because the globcover layer has a much finer resolution, we need to resample
globcover <- terra::resample(globcover,WOCO.isoscape, method="max")
WOCO.isoscape <- WOCO.isoscape*globcover
plot(WOCO.isoscape)

## extract hydrogen isotope values from that distribution
rain.d2H<-as.numeric(na.omit(terra::values(WOCO.isoscape)[,1]))
rain.d2H<-rain.d2H[rain.d2H<0]  ## remove the non-forest values (>10,0000 grid cells are removed)
mean.rain.d2H<-mean(rain.d2H, na.rm=T)
sd.rain.d2H<-sd(rain.d2H, na.rm=T)
hist(rnorm(1000,mean.rain.d2H,sd.rain.d2H))




## 2.4. use known origin data to calibrate rainwater against feather d2H  -------


### 2.4.1 convert woco data to SpatVector to extract values from raster ------

woco.vect<-terra::vect(woco.sf)


## 2.4.2. extract rainwater hydrogen isotopes for the location of known-sample woodcocks -------

woco.sf$d2h_MA<-terra::extract(isoscape,woco.vect)$d2h_MA
woco.sf$d2h_se_MA<-terra::extract(isoscape,woco.vect)$d2h_se_MA

