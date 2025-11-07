#---------------------------------------------------------------------------
### DATA PREPARATION FOR ISOTOPic ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025

## goal is to estimate proportion of shot birds coming from Switzerland

## USEFUL RESOURCES:
# workshop: https://github.com/JacksonKusack/Isotope-Assignment-Workshop
# paper: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/04-0175
# isocat package: https://github.com/cjcampbell/isocat
# assignR  example: https://cran.r-project.org/web/packages/assignR/vignettes/assignR.html
# isoriX package: https://bookdown.org/content/782/introduction.html [does not download rainfall isotope values]


## received new known origin data from Andrew Hoodless on 28 Sept 2025
## included those data for calibration


## use annual amount-weighted rainwater isotopes (advice by David Soto because it reflects food web much better)

## update after chat with David Soto on 20 Oct 2025
### NEED TO DO: re-scale the isotope metrics for all feathers (from IZW and UK) based on Soto 2017 using that equation
# Correction equation to convert from formerly assigned values for keratins where (a) CBS and KHS or (b) USGS42/43 were
# used as two‐point calibration reference materials (researchers can estimate uncertainties for slope and intercept using FREML,
#                                                    when needed)
# (a) From (7) to (1): Revised δ2
# H VSMOW = 10.774 + 0.852 * old data
## https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/rcm.7893
## revise isoscape and use full annual amount-weighted mean value rather than growing season average


## need to add elevation as prior to isotope data and filter out shot birds occurring above 2200 m (the highest breeding record in Switzerland)



rm(list=ls())
library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
#library(readxl)
library(assignR) ## https://cran.r-project.org/web/packages/assignR/vignettes/assignR.html
# library(isocat) ## https://cran.r-project.org/web/packages/isocat/vignettes/isocat.html
# library(isoAssign) ## https://rdrr.io/github/SMBC-NZP/MigConnectivity/man/isoAssign.html
library(terra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tictoc)
filter<-dplyr::filter
select<-dplyr::select
library(gridExtra)
library(tidyterra)


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock"),silent=T)



# 1. READ IN PROCESSED ISOTOPE DATA ----------------

### data are poorly documented - my reconstruction is as follows based on Voigt 2016 report:
### dH: raw isotope measurements in feather keratin (but only available for half of the samples!!) - based on Berlin lab standard
### dH_reg: corresponding rainfall hydrogen isotope ratio converted with equation in Powell (2012) (age-specific adjustment)
### dH_scaled: is the converted d2H of feather keratin using the formula 0.979*dH + 20.701 to adjust it to Canada lab standard
### dH_correct: is the corrected d2H of feathers IF a value was in column 'bc' (whatever this stands for)


woco<-fread("data/WOCO_isotopes.csv")
#woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")
table(woco$AGE)

## 1.1. PLOT THE DIFFERENT d2H isotope values ----

## raw measurements against rainwater isotopes
ggplot(woco, aes(x=dH,y=dH_reg, col=AGE)) +
  geom_point() +
  geom_smooth()

## raw measurements with Berlin standard compared to Canada standard - this is a straight line because it is based on a predictable equation
ggplot(woco, aes(x=dH_scaled,y=dH, col=AGE)) +
  geom_point() +
  geom_smooth()

## assess how many are missing
length(which(is.na(woco$dH_scaled)))
length(which(is.na(woco$dH)))

## back-convert all measurements to INTERNATIONAL standards following Soto et al 2017
# woco$dH<-(woco$dH_scaled - 20.701)/0.979 ## discontinued after discussion with David Soto on 20 Oct 2025
woco$dH<-10.774 + (0.852*woco$dH_scaled)
woco %>% dplyr::filter(is.na(dH))
# ggplot(woco, aes(x=D2,y=dH, col=AGE)) +
#   geom_point() +
#   geom_smooth()

## USE dH FOR ALL ANALYSES FROM HERE ON!!!




## 1.2. Feather fractionation from rainwater hydrogen isotope ----

## commented out after resolving the issue on 3 June 2025 - dH is the raw value we want to use

# # converting isoscape rainfall values into feather d2H values
# # based on Powell 2012, Table 3.1, page 133 with a sample size of 135 feathers
# n <- 135
# intcept<- rnorm(100,4.50799, (6.34005*sqrt(n)))
# beta.rain<- rnorm(100,0.79760, (0.08049*sqrt(n)))
# beta.age<- rnorm(100,-28.73077, (2.14900*sqrt(n)))
# 
# 
# #woco$dH <- intcept + beta.rain*d2Hrain + beta.age * juvenile  ## the transformation equation of Powell based on mean annual deuterium (MAD)
# woco$dH_Powell <- 4.50799 + 0.79760*woco$dH - 28.73077 * woco$JUVENIL
# 
# 
# ggplot(woco, aes(x=dH_reg,y=dH_Powell, col=AGE)) +
#   geom_point() +
#   geom_smooth()
# 
# ggplot(woco, aes(x=dH_correct,y=dH_Powell, col=AGE)) +
#   geom_point() +
#   geom_smooth()
# 
# 
# ggplot(woco, aes(x=dH_scaled,y=dH_Powell, col=AGE)) +
#   geom_point() +
#   geom_smooth()
# 
# 
# ## try and figure out how values were created
# #summary(lm(dH_correct~dH+AGE, data=woco))
# summary(lm(dH_correct~dH_reg+AGE, data=woco))  ### this seems to be the appropriate equation
# #summary(lm(dH_correct~dH_scaled+AGE, data=woco))
# #summary(lm(dH_Powell~dH_reg+AGE, data=woco))
# 

## USE dH FOR ALL ANALYSES FROM HERE ON!!!



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


## 2.4.3. look at relationship between rainwater to d2H of feather -------
ggplot(woco.sf, aes(x=d2h_MA,y=dH, col=AGE, fill=AGE)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="d2H in rainwater (mean annual average) at sampling location", y="d2H in woodcock feather") +
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_line(linewidth=0.4, colour="grey89", linetype="dashed"),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=18),
        legend.position="inside",
        legend.position.inside=c(0.10,0.90),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16, color="black"))
#ggsave("output/REVISED_known_WOCO_feather_isotope_calibration.jpg", width=9, height=8)

# ISO_CALIB<-lm(dH~d2h_MA+AGE, data=woco.sf)
# summary(ISO_CALIB)




### 2.4.4. apply that equation to all other shot woodcocks -------
## deterministic approach without considering uncertainty
woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)

woco.unk.sf$d2h_MA<-terra::extract(isoscape,woco.unk.vect)$d2h_MA
woco.unk.sf$d2h_se_MA<-terra::extract(isoscape,woco.unk.vect)$d2h_se_MA
woco.unk.sf %>% filter(is.na(AGE))


## 2.5. include phenology data to integrate probability of non-local origin based on date a bird was shot -------

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

# ## draw a plot ##
# colors <- c("All birds" = "darkolivegreen", "Local birds" = "firebrick", "Shot birds" = "gray23")
# 
#   
# FIG_s3<-ggplot()+
#   geom_line(data=woco_mig, aes(x=Date, y=mig, color="Local birds"),linewidth=1) +
#   geom_line(data=woco_abundance, aes(x=Date, y=abund, color="All birds"),linewidth=1) +
#   geom_col(data=woco_shot, aes(x=Date, y=abund, color="Shot birds"),width = 6, alpha=0.5) +
#   labs(color = "Origin") +
#   scale_color_manual(values = colors) +
#   
#   ## format axis ticks
#   scale_x_date(name="Date", date_labels = "%d %b") +
#   scale_y_continuous(name="Relative abundance of woodcocks", limits=c(0,1), breaks=seq(0,1,0.2)) +
#   
#   ## beautification of the axes
#   theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.y=element_text(size=14, color="black"),
#         axis.text.x=element_text(size=14, color="black"), 
#         axis.title=element_text(size=18),
#         legend.direction = "vertical",
#         legend.box = "horizontal",
#         legend.title=element_text(size=16, color="black"),
#         legend.text=element_text(size=14, color="black"),
#         legend.position="inside",
#         legend.key = element_rect(fill = NA, color = NA),
#         legend.background = element_rect(fill = NA, color = NA),
#         legend.position.inside=c(0.1,0.5),
#         strip.background=element_rect(fill="white", colour="black"))
# FIG_s3
# 
# # ggsave(plot=FIG_s3,
# #        filename="output/woco_phenology_comparison.jpg", 
# #        device="jpg",width=11, height=8)
# 





# 3. clean up and save workspace -------


## 3.1. prepare the data needed for NIMBLE input ----
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


rm(isoscape, globcover,forest.mat)
gc()

save.image("data/woco.input.data.RData")




## 3.2. remove non-SUI calibration data and create reduced input data ----

woco.sf <- woco.sf %>%
  filter(PROVENANCE_voigt=="Suisse")
dim(woco.sf)
save.image("data/woco.reduced.input.data.RData")


