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


woco.unk <- woco.unk.sf %>%
  st_drop_geometry()


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
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium","Switzerland"))

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
outscape<-WOCO.isoscape17 %>% mutate(N=0)
for (f in 1:dim(woco.unk.sf)[1]) {  ## loop across all feathers
  
  x<-woco.unk[f]
  iscape<-isoscapes[[woco.unk.sf$Year[f]-2012]]
  iscape$mean<-predict(ISO_CALIB,newdat=data.frame(d2h_MA=as.vector(iscape$mean),AGE=x$AGE)) ### convert to feather isotopes with equation above
  
  p_prob <- app(iscape[[c("mean", "mean_predVar")]], fun = function(vals) {
  # vals is a matrix: nrows = number of cells in the chunk, ncols = number of layers (2 here)
  m <- as.matrix(vals)[, 1]
  v <- as.matrix(vals)[, 2]
  
  # Keep variance non-negative, compute sd vectorized
  sd <- sqrt(pmax(v, 0))
  
  # Compute pnorm vectorized; returns a vector of same length as m
  res <- pnorm(x$dH, mean = m, sd = sd)
  
  # Optional: mark invalid sd as NA (vectorized)
  bad <- (sd == 0) | is.na(sd)
  res[bad] <- NA_real_
  
  # Must return a numeric vector with length equal to nrow  # Must return a numeric vector with length equal to nrow(vals)
  res
  
})

  ### APPLY THE 2:1 RATIONALE OF HOBSON AS ORIGIN 1 or 0
  
  thresh<-quantile(as.vector(p_prob),0.333, na.rm=T)
  outscape$N<-as.vector(outscape$N) + as.vector(ifelse(as.vector(p_prob)>thresh,1,0))
}



