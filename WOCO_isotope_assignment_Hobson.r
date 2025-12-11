#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### ORIGIN ASSIGNMENT OF WOODCOCKS SHOT IN SWITZERLAND ----
## attempt to recreate the assignment method of Hobson et al. 2013
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

## also included annual isoscape suggested by David Soto

rm(list=ls())
library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
library(terra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tictoc)
library(basicMCMCplots) # for trace plots called chainsPlot
filter<-dplyr::filter
select<-dplyr::select
library(gridExtra)
library(tidyterra)


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock"),silent=T)



# 1. READ IN PROCESSED ISOTOPE DATA -------------------------------------------------------------------------------------
# prepared in WOCO_rainfall_isoscape_creation.r
load("data/woco.input.data.annual.RData")
#load("data/woco.reduced.input.data.RData") ## without the calibration feathers from Hoodless and Powell
# try(rm(isoscape, globcover), silent=T)
# ignore the error referring to C++ https://github.com/keblu/MSGARCH/issues/48

## 1.1. modify data to include abundance predictions from Ornitho.ch records ----------
## this approach is intended to address the dilution effect (when local birds are outnumbered by migrants)

woco.unk.abd.prior <- woco.unk.sf %>%
  st_drop_geometry() %>%
  mutate(DAY=yday(feather_sampling_date), DAY_2=(yday(feather_sampling_date)^2))

woco_abundance <- woco_abundance %>%
  mutate(DAY=yday(Date), DAY_2=(yday(Date)^2)) %>%
  dplyr::filter(DAY >= min(woco.unk.abd.prior$DAY)) %>%
  dplyr::filter(DAY <= max(woco.unk.abd.prior$DAY))


m2<-glm(SOPM ~ DAY + DAY_2, data = woco_abundance, family = "poisson")
#m2<-loess(SOPM ~ DAY, data = woco_abundance, family = "poisson")
woco.unk.sf$abd_prior <- predict(m2, newdat=woco.unk.abd.prior, type="response")  ## this is the index of abundance over time - if many woodcocks are there the probability of harvesting a local one is smaller
woco.unk.sf$abd_prior <- woco.unk.sf$abd_prior/ceiling(max(woco.unk.sf$abd_prior)) ## scale to a proportion, but round up maximum because beta prior cannot take values of 1
summary(woco.unk.sf$abd_prior)
hist(woco.unk.sf$abd_prior)

ggplot(data=woco_abundance, aes(x=DAY, y=SOPM/max(SOPM))) +
  geom_point(data=woco_abundance, aes(x=DAY, y=SOPM/max(SOPM))) +
  geom_smooth(method=loess, se=T) +
  geom_line(data=woco.unk.sf, aes(x=yday(feather_sampling_date),y=abd_prior), col="firebrick")


## to get initial values right
ISO_CALIB<-lm(dH~d2h_MA+AGE, data=woco.sf)
summary(ISO_CALIB)




bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)

woco.countries.orig <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

WOCO.isoscape13<-readRDS("./data/isoscape13.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape14<-readRDS("./data/isoscape14.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape15<-readRDS("./data/isoscape15.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape16<-readRDS("./data/isoscape16.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape17<-readRDS("./data/isoscape17.rds") %>%
  terra::mask(woco.countries.orig)


unique(woco.unk.sf$Year)
isoscapes<-list(WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17)



# 2. SPECIFY COMBINED PROBABILITY MODEL TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------

## concept of approach
for (f in 1:842) {  ## loop across all feathers
  
  x<-woco.unk.sf$dH[f]
  iscape<-isoscapes[[woco.unk.sf$Year[f]-2012]]
  for(n in 1:ncells(iscape)) {
    
    xr<- predict(ISO_CALIB,x)
    pnorm(x, mean = iscape$dH[n], sd = sigma, lower.tail = FALSE)
    
    
  }
  
  
  
  
  
  
}



