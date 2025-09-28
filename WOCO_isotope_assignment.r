### ISOTOPic ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025

## goal is to estimate proportion of shot birds coming from Switzerland

## USEFUL RESOURCES:
# workshop: https://github.com/JacksonKusack/Isotope-Assignment-Workshop
# paper: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/04-0175
# isocat package: https://github.com/cjcampbell/isocat
# assignR  example: https://cran.r-project.org/web/packages/assignR/vignettes/assignR.html


## MAJOR ADDITION ON 5 June 2025: added JAGS model to estimate probability of local origin for shot woodcock
## MAJOR ADDITION ON 11 June 2025: added GlobCover data to extract isotope distribution from forest areas in NW Europe

## UPDATE 17 June 2025:

# Pius: include temporal phenology of abundance and date when bird was shot to alter probability of foreign origin based on date (later more likely to be non-local)
# Murdy: create map for each canton that shows what region is covered by the term 'local' based on equivalent rain d2H
# Elli: get Swiss high resolution isotope data and create calibration for rainwater-feather isotope ratios based on actual local simultaneous measurements. This calibration then cannot be used outside of Switzerland, so only purpose would be to show potential differences in calibration?
# Urs: provide overall mean probability 

## UPDATE 3 July 2025:
## tried to add a time distribution for monitoring data

## received new known origin data from Andrew Hoodless on 28 Sept 2025
## included those data for calibration


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
library(coda)
library(MCMCvis)
library(nimble)
library(tictoc)
library(basicMCMCplots) # for trace plots called chainsPlot
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
### dH_correct: is the corrected d2H of feathers IF a value was in column 'bc' (whatever the fuck this stands for)


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

## back-convert all measurements to IZW (Berlin) standards
woco$dH<-(woco$dH_scaled - 20.701)/0.979
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
  rename(ID=RefID, KANTON=Region,
         DATE=Date,dH=DHF,LATITUDE=Latitude,LONGITUDE=Longitude) %>%
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



# 2. ESTIMATING PROPORTION OF LOCAL BIRDS IN HUNTING BAG ---------------------

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

isoscape <- readRDS("data/global_d2H_GS_isoscape.rds") %>%
  crop(extent(SUI))

## downloaded from Nelson et al 2021, but poor resolution
# isoscape <- raster("data/Piso.AI_v1.2020_0.5deg_1950-2020.nc") %>%
#   crop(extent(SUI))



### 2.3.1. curtail woodcock distribution to potential origin countries and FOREST
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Switzerland","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

## create a conversion matrix
forest.mat<-matrix(0, nrow=3, ncol=3)
forest.mat[,1]<-c(11,40,110)  ## from values for conversion matrix
forest.mat[,2]<-c(30,100,230)  ## to values for conversion matrix
forest.mat[,3]<-c(0,1,0)  ## replacement values for conversion matrix

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
WOCO.isoscape <- readRDS("data/global_d2H_GS_isoscape.rds") %>%
  crop(woco.countries)
plot(WOCO.isoscape)

## remove all non-forest areas by multiplying with 0
crs(globcover)==crs(WOCO.isoscape)
origin(globcover)
origin(WOCO.isoscape)

## because the globcover layer has a much finer resolution, we need to resample
globcover <- terra::resample(globcover,WOCO.isoscape, method="max")
WOCO.isoscape <- WOCO.isoscape*globcover

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

woco.sf$d2h_GS<-terra::extract(isoscape,woco.vect)$d2h
woco.sf$d2h_se_GS<-terra::extract(isoscape,woco.vect)$d2h.se


## 2.4.3. create equation to link rainwater to d2H of feather -------
# SUPERSEDED BY NIMBLE MODEL BELOW
ggplot(woco.sf, aes(x=d2h_GS,y=dH, col=AGE, fill=AGE)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="d2H in rainwater (growing season average) at sampling location", y="d2H in Swiss woodcock feather") +
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
#ggsave("output/SUI_WOCO_feather_isotope_calibration.jpg", width=9, height=8)

ISO_CALIB<-lm(dH~d2h_GS+AGE, data=woco.sf)
summary(ISO_CALIB)




### 2.4.4. apply that equation to all other shot woodcocks -------
## deterministic approach without considering uncertainty
woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)

woco.unk.sf$d2h_GS<-terra::extract(isoscape,woco.unk.vect)$d2h
woco.unk.sf$d2h_se_GS<-terra::extract(isoscape,woco.unk.vect)$d2h.se
woco.unk.sf %>% filter(is.na(AGE))
# 
# d2H_pred<-predict(ISO_CALIB,newdata=woco.unk.sf,se.fit=TRUE,interval="prediction",level=0.95, type = "response")$fit
# 
# 
# 
# ## 2.4.5. test whether each shot feather came from within location-specific prediction interval -------
# 
# OUT_SUMMARY<-bind_cols(woco.unk.sf,d2H_pred) %>%
#   filter(!is.na(AGE)) %>%
#   filter(!is.na(dH)) %>%
#   st_drop_geometry() %>%
#   mutate(LOCAL=ifelse(dH>lwr & dH<upr,1,0)) %>%
#   group_by(AGE, KANTON) %>%
#   summarise(prop.local=mean(LOCAL, na.rm=T)) %>%
#   spread(key=AGE, value=prop.local) %>%
#   rename(prop_adult_local=Adulte, prop_juvenile_local=Juvénile)
# 
# OUT_SUMMARY<-bind_cols(woco.unk.sf,d2H_pred) %>%
#   filter(!is.na(AGE)) %>%
#   filter(!is.na(dH)) %>%
#   st_drop_geometry() %>%
#   mutate(LOCAL=ifelse(dH>lwr & dH<upr,1,0)) %>%
#   group_by(AGE, KANTON) %>%
#   summarise(n=length(LOCAL)) %>%
#   spread(key=AGE, value=n, fill=0) %>%
#   rename(N_adult=Adulte, N_juvenile=Juvénile) %>%
#   left_join(OUT_SUMMARY, by="KANTON")
# 
# OUT_SUMMARY
# 
# #fwrite(OUT_SUMMARY, "output/prop_local_WOCO_shot_isotopes.csv")




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

## draw a plot ##
colors <- c("All birds" = "darkolivegreen", "Local birds" = "firebrick", "Shot birds" = "gray23")

  
FIG_s3<-ggplot()+
  geom_line(data=woco_mig, aes(x=Date, y=mig, color="Local birds"),linewidth=1) +
  geom_line(data=woco_abundance, aes(x=Date, y=abund, color="All birds"),linewidth=1) +
  geom_col(data=woco_shot, aes(x=Date, y=abund, color="Shot birds"),width = 6, alpha=0.5) +
  labs(color = "Origin") +
  scale_color_manual(values = colors) +
  
  ## format axis ticks
  scale_x_date(name="Date", date_labels = "%d %b") +
  scale_y_continuous(name="Relative abundance of woodcocks", limits=c(0,1), breaks=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title=element_text(size=16, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.position="inside",
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.position.inside=c(0.1,0.5),
        strip.background=element_rect(fill="white", colour="black"))
FIG_s3

# ggsave(plot=FIG_s3,
#        filename="output/woco_phenology_comparison.jpg", 
#        device="jpg",width=11, height=8)










# 3. COMBINED PROBABILITY MODEL TO ESTIMATE PROB OF LOCAL ORIGIN IN NIMBLE ----

## 3.1. specify model in NIMBLE ----

woco.orig.model<-nimbleCode({
  
  # Parameters:
  # p.loc: probability of local origin (that shot woodcock was a local bird) - varies by age and canton
  
  # b.age: regression parameter for age that translates rainfall d2H to feather d2H 
  # b.rain: regression parameter for rainfall d2H that translates rainfall d2H to feather d2H
  # int.rain: regression intercept that translates rainfall d2H to feather d2H 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CALIBRATION REGRESSION FOR KNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Likelihood:  
  for (i in 1:nind.known){
        d2H_feather.known[i] ~ dnorm(mu.known[i], sd=sigma.calib)
        mu.known[i] <- int.rain + b.age*age.known[i] + b.rain*d2H_rain.known[i]
      }
    
  # Priors informed by known origin woodcock feathers in UK (Powell thesis 2012)
  int.rain ~ dnorm(4.5, sd=5) # informative prior based on Powell
  b.age ~ dnorm(-28.7, sd=5) # informative prior based on Powell
  b.rain ~ dnorm(0.8, sd=1) # informative prior based on Powell
  sigma.calib ~ dunif(0,20) # standard deviation
  
  
  
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # TIME DISTRIBUTION FOR OVERALL ABUNDANCE
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # negative binomial distribution:  
  # for (i in 1:n.count.weeks){
  #   abund.known[i] ~ dnegbin(mu.abd.known[i], r.abd)
  #   mu.abd.known[i] <- int.abd + b.countday*count.day[i] + b2.countday*count.day.sq[i]
  # }
  # 
  # # Priors not informed by anything
  # int.abd ~ dnorm(0.2, sd=5) # uninformative prior
  # b.countday ~ dnorm(1, sd=5) # uninformative prior
  # b2.countday ~ dnorm(1, sd=5) # uninformative prior for quadratic effect
  # r.abd ~ dgamma(0.01,0.01) # dispersion parameter for neg bin distribution
  
  
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # PROBABILITY DISTRIBUTION FROM DEPARTURE MODEL
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # for (t in 1:(nweeks)){
  #   logit.mig[t] <- lm.mean +      ### intercept from migration model
  #     b.mig.week*(week[t])     ### week effect from migration model
  #   mig.prob[t] <-ilogit(logit.mig[t])
  # }    
  #     


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #### RANDOM EFFECT FOR EACH CANTON
  # for (ct in 1:ncanton){
  #   cant.reff[ct] ~ dnorm(0,sd=1)   ##  random normal effect centred on zero to get canton-specific probabilities
  # } #a

  
  #### OVERALL PRIOR PROBABILITY OF NON-LOCAL ORIGIN BASED ON AGE CLASS
  # #for (ct in 1:ncanton){
  #   for (a in 1:2){
  #     mean.p.nonlocal.iso[a] ~ dbeta(2,2)   ##  hist(rbeta(1000,2,2)) very low and very high probabilities are less likely a priori
  #     logit.mean.p.nonlocal[a]<-log(mean.p.nonlocal.iso[a] / (1-mean.p.nonlocal.iso[a]))
  #   } #a
  # #} #ct
  
  
  #### LINEAR PREDICTOR OF AGE AND CANTON-SPECIFIC PROBABILITY OF NON-LOCAL ORIGIN
  # moved to individual loop below
  # for (ct in 1:ncanton){
  #   for (a in 1:2){
  #     logit.p.nonlocal[ct,a] <- logit.mean.p.nonlocal[a] + cant.reff[ct]  ##  combination of overall mean and random effect
  #     p.nonlocal.iso[ct,a]<-ilogit(logit.p.nonlocal[ct,a])
  #   } #a
  # } #ct
  
  
  # Likelihood 
  sd.unknown[2]<-sigma.calib+sd.rain.d2H  ### overall distribution across Europe
  sd.unknown[1]<-sigma.calib
  
  for (i in 1:nind.unkn){
    
    # linear predictor of isotope non-local probability
    p.nonlocal.iso[i] ~ dunif(0,1)   ##  hist(rbeta(1000,2,2)) very low and very high probabilities are less likely a priori
    # logit.mean.p.nonlocal[a]<-log(mean.p.nonlocal.iso[a] / (1-mean.p.nonlocal.iso[a]))
    # 
    # logit.p.nonlocal.iso[i] <- logit.mean.p.nonlocal[age.unknown[i]+1] + cant.reff[canton[i]]  ##  combination of overall mean and random effect
    # p.nonlocal.iso[i]<-ilogit(logit.p.nonlocal.iso[i])
    
    # Latent indicator for isotope distribution
    z[i] ~ dbern(p.nonlocal.iso[i])  ## age is specified as 0 and 1 in data, need to add 1 for index to work
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H    ### overall distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*d2H_rain.unknown[i]  ### local expected feather isotope distribution

    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[z[i]+1])
    
    # Potential local distribution based on timing
    # logit.mig.unknown[i,2] <- int.abd + b.countday*unk.day[i] + b2.countday*unk.day.sq[i]    ### overall abundance distribution from abundance data
    # mig.prob[i,2] <-ilogit(abd.unknown[i,2])   ### predicted probability from departure model
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal.time[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from departure model

    
    # combine the probabilities of migration and isotope into a single probability using log odds
    # Convert each probability to odds:
    odds_iso[i] <- p.nonlocal.iso[i] / (1 - p.nonlocal.iso[i])
    odds_time[i] <- p.nonlocal.time[i]  / (1 - p.nonlocal.time[i] )
    
    # Combine odds
    combined_odds[i] <- odds_iso[i] * odds_time[i]
    
    # Convert back to probability
    p.nonlocal[i] <- combined_odds[i] / (1 + combined_odds[i])

    
    
    # #### CALCULATE SUMS PER AGE AND CANTON FOR DERIVED MEAN PROBABILITY ######
    # 
    # for (a in 1:2){
    #   
    #   age.count[i,a] <- equals(age.unknown[i], a-1)
    #   age.sum.p[i,a] <- p.nonlocal[i] * age.count[i,a]
    #   
    #   for (ct in 1:ncanton){
    #     
    #     is_canton[i,a,ct] <- equals(canton[i], ct)
    #     ct.sum.p[i,a,ct] <- p.nonlocal[i] * ct.count[i,a,ct]
    #     ct.count[i,a,ct] <- (age.count[i,a] * is_canton[i, a,ct])
    #     
    #   } #ct
    # } #a
    
  } #i
  
  
  # #### DERIVED OVERALL PROBABILITY BY AGE CLASS AND CANTON
  # this does not work in Nimble because it cannot sum across indices
  # 
  # for (a in 1:2){
  #   mean.p.nonlocal[a]<-sum(age.sum.p[1:nind.unkn,a]) / sum(age.count[1:nind.unkn,a])
  #   
  #   for (ct in 1:ncanton){
  #     p.nonlocal[ct,a]<-sum(ct.sum.p[1:nind.unkn,a,ct]) / sum(ct.count[1:nind.unkn,a,ct])
  #     
  #   } #ct
  # } #a

}) ## end of nimble code chunk






## 3.2. prepare the data needed for NIMBLE input ----
woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  #filter(KANTON !="VS") %>% ## remove Valais because only 6 birds from 1 age class, causes imbalance in data
  filter(!is.na(d2h_GS))

woco.sf <- woco.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_GS))

table(woco.unk.sf$AGE,woco.unk.sf$KANTON)

#### DISTINGUISH CONSTANTS AND DATA----
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      nind.known = dim(woco.sf)[1],
                      mean.rain.d2H = mean.rain.d2H,
                      sd.rain.d2H = sd.rain.d2H,
                      d2H_rain.known=woco.sf$d2h_GS,
                      d2H_rain.unknown=woco.unk.sf$d2h_GS,
                      #ncanton=length(unique(woco.unk.sf$KANTON)),
                      age.unknown = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      age.known = ifelse(woco.sf$AGE=="Adulte",0,1),
                      #canton = as.numeric(as.factor(woco.unk.sf$KANTON)),
                      lm.mean.mig=logit(readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[2,1]),
                      b.mig.week=readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[1,1],
                      #unk.day=yday(UNK_WC$DATE),
                      unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26"))   ## migration weeks start only in August
                      #count.day=yday(woco_abundance$Date),
                      #count.day.sq=yday(woco_abundance$Date)^2,
                      #n.count.weeks=length(woco_abundance$SOPM)
                      )

iso.data <- list(d2H_feather.known = woco.unk.sf$dH,
                 #abund.known=woco_abundance$SOPM,
                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.3. specify NIMBLE run settings ----

# Parameters monitored
parameters.iso <- c("b.rain","b.age","int.rain","sigma.calib","p.nonlocal","p.nonlocal.iso","p.nonlocal.time")
                    #"int.abd","b.countday","b2.countday","r.abd","mean.p.nonlocal","p.nonlocal",)

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_GS-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,4.5, sd=5), # informative prior based on Powell
                  b.age = rnorm(1,-28.7, sd=5), # informative prior based on Powell
                  b.rain = rnorm(1,0.8, sd=1), # informative prior based on Powell
                  sigma.calib = runif(1,0,20), # standard deviation
                  
                  ## INITS FOR ABUNDANCE REGRESSION
                  # int.abd = rnorm(1,0.2, sd=5), # uninformative prior
                  # b.countday = rnorm(1,1, sd=5), # uninformative prior
                  # b2.countday = rnorm(1,1, sd=5), # uninformative prior for quadratic effect
                  # r.abd  = rgamma(1,0.01,0.01), # dispersion parameter for neg bin distribution
                  
                  #### FLAT PRIOR PROBABILITY FOR BEING LOCAL BIRD
                  p.nonlocal.iso=runif(dim(woco.unk.sf)[1],0,1)  ###  hist(rbeta(1000,2,4.5))
                  #cant.reff=rnorm(iso.constants$ncanton,0,1)
                  #p.loc=matrix(rbeta(length(unique(woco.unk.sf$KANTON))*2,2,4.5),ncol=2)  ###  hist(rbeta(1000,2,4.5))
)

#iso.inits$p.nonlocal<-matrix(rbeta(iso.constants$ncanton*2,2,2),ncol=2)



# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
# test <- nimbleModel(code = woco.orig.model,
#                     constants=iso.constants,
#                     data = iso.data,
#                     inits = iso.inits,
#                     calculate=TRUE)
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
# test$initializeInfo()
#help(modelInitialization)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 50000
n.burnin <- 25000
n.chains <- 4



## 3.4. run model in NIMBLE ------------------------------------------
## this takes 5 min

tic()
### this takes 4000 sec (2 hrs) for 20000 iterations and converges in that time
woco.iso <- nimbleMCMC(code = woco.orig.model,
                        constants=iso.constants,
                        data = iso.data,
                        inits = iso.inits,,
                        monitors = parameters.iso,
                        thin=4,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains,
                        progressBar = getNimbleOption("MCMCprogressBar"),
                        summary=T)
toc() 



saveRDS(woco.iso,"output/woco_iso_origin_model.rds")
#woco_surv<-readRDS("output/woco_iso_origin_model.rds")





## 3.5. assess model and examine convergence ------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("p.nonlocal","b.rain","b.age","int.rain","sigma.calib"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_iso_origin_parm_estimates.csv")
#out<-fread("output/woco_iso_origin_parm_estimates.csv")


#MCMCplot(woco.iso$samples, params=c("mean.p.nonlocal"))


## look at the chains and whether they mixed well
chainsPlot(woco.iso$samples,
           var=c("b.rain","b.age","int.rain","sigma.calib"))



## 3.5.1 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples <- rbind(woco.iso$samples$chain1,woco.iso$samples$chain2,woco.iso$samples$chain3,woco.iso$samples$chain4)
mean.p.nonlocal <- as_tibble(samples[,grep("p.nonlocal", colnames(samples))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) 

hist(mean.p.nonlocal$p.nonlocal)




## 3.5.2 compare probabilities between isotope and time ------------------------------------------

p.comp<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("p.nonlocal.iso","p.nonlocal.time")))
p.comp$parameter<-row.names(p.comp)
names(p.comp)[c(3,4,5)]<-c('lcl','median', 'ucl')
p.comp<- p.comp %>% mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) %>%
  mutate(info=if_else(grepl("time", parameter)==T,"time","isotope")) %>%
  select(ind,age,ctn,info,median) %>%
  pivot_wider(.,names_from = info,values_from = median)


## plot the probabilities
p.comp %>% 
ggplot(aes(x=isotope, y=time,col=ctn))+
  geom_point(size=1)+
  
  ## format axis ticks
  xlab("p.nonlocal isotope")+
  ylab("p.nonlocal time") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"),
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"),
        strip.text.y=element_text(size=18, color="black"),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        strip.background=element_rect(fill="white", colour="black"))




## 3.6. summarise output for manuscript text ------------------------------------------
## changed after model includes mean and random effect
# out %>%
#   filter(!parameter %in% c("b.rain","b.age","int.rain","sigma.calib")) %>%
#   mutate(Age=as.numeric(substr(parameter,nchar(parameter)-2,nchar(parameter)-1))) %>%
#   mutate(Age=ifelse(Age==1,"Adult","Juvenile")) %>%
#   group_by(Age) %>%
#   summarise(m=mean(mean),l=min(lcl), u=max(ucl))

agemeans<-out %>%
  filter(parameter %in% c("mean.p.nonlocal[1]","mean.p.nonlocal[2]")) %>%
  mutate(Age=as.numeric(substr(parameter,nchar(parameter)-1,nchar(parameter)-1))) %>%
  mutate(Age=ifelse(Age==1,"Adult","Juvenile"))
agemeans  


## 3.7. summarise output in graphical form ------------------------------------------

# LOAD AND MANIPULATE ICONS TO REMOVE BACKGROUND
require(png)
library(grid)
library(gtable)
require(jpeg)
library(magick)
imgGun<-readPNG("manuscript/rifleicon.png")
gunicon <- rasterGrob(imgGun, interpolate=TRUE)
imgWOCO<-image_read("manuscript/woodcock.jpg") %>% image_transparent("white", fuzz=5)
wocoicon <- rasterGrob(imgWOCO, interpolate=TRUE)


FIGURE<-mean.p.nonlocal %>%
  group_by(age,ctn,ind) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(age,ctn) %>%
  summarise(for.med=median(p.nonlocal.mean),for.ucl=quantile(p.nonlocal.mean,0.025), for.lcl=quantile(p.nonlocal.mean,0.975)) %>%
  
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  #mutate(Kanton=levels(as.factor(woco.unk.sf$KANTON))[ctn]) %>%
  
  ggplot(aes(x=ctn, y=for.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=for.lcl, ymax=for.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  
  ## format axis ticks
  labs(y="Prop. of shot woodcock that are not local",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  # scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  # scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  # 
  ### add the bird icons
  annotation_custom(gunicon, xmin=5.2, xmax=5.8, ymin=0.88, ymax=0.96)+
  annotation_custom(wocoicon, xmin=5.5, xmax=6.9, ymin=0.85, ymax=1.02)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major.y = element_line(linewidth=0.5, colour="grey59", linetype="dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14, color="black"),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title=element_text(size=14, color="black"),
        legend.position="inside",
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.position.inside=c(0.65,0.85),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE

ggsave(plot=FIGURE,
       filename="output/woco_iso_time_origin_probability_estimates.jpg", 
       device="jpg",width=11, height=8)


# FIGURE<-out %>%
#   filter(!parameter %in% c("b.rain","b.age","int.rain","sigma.calib","mean.p.nonlocal[1]","mean.p.nonlocal[2]")) %>%
#   mutate(Age=as.numeric(substr(parameter,nchar(parameter)-2,nchar(parameter)-1))) %>%
#   mutate(Age=ifelse(Age==1,"Adult","Juvenile")) %>%
#   mutate(Canton=readr::parse_number(parameter,locale=locale(grouping_mark=". ", decimal_mark=","))) %>%  ### parse_number doesn't deal with dots: https://stackoverflow.com/questions/61328339/r-parse-number-fails-if-the-string-contains-a-dot
#   mutate(Kanton=levels(as.factor(woco.unk.sf$KANTON))[Canton]) %>%
#   mutate(for.med=(median),for.ucl=lcl, for.lcl=ucl) %>%
#   
#   ggplot(aes(x=Canton, y=for.med))+
#   geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
#   geom_errorbar(aes(ymin=for.lcl, ymax=for.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
#   
#   ## format axis ticks
#   labs(y="Prop. of shot woodcock that are not local",x="Swiss Canton",col="") +
#   scale_x_continuous(limits=c(0.5,6.5), breaks=c(1:6), labels=levels(as.factor(woco.unk.sf$KANTON))) +
#   # scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
#   # scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
#   # 
#   ### add the bird icons
#   annotation_custom(gunicon, xmin=5.2, xmax=5.8, ymin=0.88, ymax=0.96)+
#   annotation_custom(wocoicon, xmin=5.5, xmax=6.9, ymin=0.85, ymax=1.02)+
#   
#   ## beautification of the axes
#   theme(panel.background=element_rect(fill="white", colour="black"),
#         panel.grid.major.y = element_line(linewidth=0.5, colour="grey59", linetype="dashed"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.y=element_text(size=14, color="black"),
#         axis.text.x=element_text(size=14, color="black"),
#         axis.title=element_text(size=16),
#         legend.text=element_text(size=14, color="black"),
#         legend.direction = "vertical",
#         legend.box = "horizontal",
#         legend.title=element_text(size=14, color="black"),
#         legend.position="inside",
#         legend.key = element_rect(fill = NA, color = NA),
#         legend.background = element_rect(fill = NA, color = NA),
#         legend.position.inside=c(0.75,0.15),
#         strip.text=element_text(size=18, color="black"),
#         strip.background=element_rect(fill="white", colour="black"))
# FIGURE
# 
# ggsave(plot=FIGURE,
#        filename="output/woco_iso_origin_probability_estimates.jpg", 
#        device="jpg",width=11, height=8)
# 


# 4. CREATE MAPS FOR EACH CANTON WHAT LOCAL RAINFALL ENCOMPASSES -----------
plot_list <- list()
for (ct in 1:length(unique(woco.unk.sf$KANTON))){

  ## get canton-wise distribution
  woco.cnt <- woco.unk.sf %>% dplyr::filter(KANTON==unique(woco.unk.sf$KANTON)[ct]) 
  cnt.iso <-  woco.cnt %>% st_drop_geometry() %>%
    dplyr::select(ID,KANTON,d2h_GS,d2h_se_GS) %>%
    rowwise() %>%
    mutate(pot.orig.d2H=rnorm(1,d2h_GS,d2h_se_GS)) %>%
    ungroup() %>%
    group_by(KANTON) %>%
    summarise(min=min(pot.orig.d2H), max=max(pot.orig.d2H))
  
  ## EXTRACT ISOSCAPE CELLS FALLING INTO THIS RANGE
  CNT.isoscape <- ifel(WOCO.isoscape >= cnt.iso$min & WOCO.isoscape <= cnt.iso$max, 1, 0)
  
  ## create plot
  plot_list[[ct]] <- ggplot(EUR) +
    geom_sf() +
    tidyterra::geom_spatraster(data=CNT.isoscape, aes(fill=d2h))+
    geom_sf(data=woco.cnt,color="red") +
    geom_sf(data=EUR, colour="grey12", fill=NA) +
    #ggtitle(unique(woco.unk.sf$KANTON)[ct]) +
    
    annotation_custom(wocoicon, xmin=-10, xmax=5, ymin=65, ymax=78) +
    annotate("text", x = 35, y = 75, label = unique(woco.unk.sf$KANTON)[ct], size=8, colour="darkolivegreen") +
    scale_fill_gradient(low = 'lightgray', high = 'lightgreen', na.value=NA) +
    
    ## beautification of the axes
    theme(panel.background=element_rect(fill="white", colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=12, color="black"),
          axis.title=element_blank(),
          legend.position="none",
          plot.margin= unit(rep(.5, 4), "lines"),
          strip.text=element_text(size=18, color="black"),
          strip.background=element_rect(fill="white", colour="black"))

} #ct

grid.arrange(grobs=plot_list,ncol=3)
ggsave(filename="output/woco_iso_origin_local_maps.jpg", 
       device="jpg",width=8, height=12)

# # 5. Geographic assignment using assignR ----
# 
# ## 5.1. load shapefile of Switzerland and EUROPE ----
# ### OPTION TO CURTAIL TO FOREST AREAS - but only needed for blind assignment
# 
# SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
#   st_transform(st_crs('EPSG:4326')) # Project to WGS84
# plot(SUI)
# 
# bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
# EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
#   st_crop(bbox) %>%
#   st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
#   dplyr::select(admin,name,adm0_a3,geometry)
# plot(EUR)
# 
# 
# 
# 
# ## 5.2. plot the known origin points of WOCO that were sampled ----
# 
# woco.sf <- ORIG_WC %>%
#   st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
# 
# 
# ggplot(SUI) +
#   geom_sf() +
#   geom_sf(data=woco.sf,color="red")
# 
# 
# 
# 
# 
# ## 5.3. extract isoscape ----
# 
# ### this function throws a weird error - downloaded and ran function independently without problem
# # isoscape <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200) %>%   ## we use MA because that is what Powell's conversion is based on
# #   projectRaster(crs = CRS(SRS_string = 'EPSG:4326'))
# ## go to "C:/STEFFEN/Vogelwarte/assignR" and open assignR.Rproj
# ## then run the script "DOWNLOAD_ISOSCAPE.R"
# 
# isoscape <- readRDS("data/global_d2H_GS_isoscape.rds") %>%
#   crop(extent(SUI))
# 
# ## downloaded from Nelson et al 2021, but poor resolution
# # isoscape <- raster("data/Piso.AI_v1.2020_0.5deg_1950-2020.nc") %>%
# #   crop(extent(SUI))
# 
# plot(isoscape, xlab="Longitude", ylab="Latitude")
# 
# 
# 
# 
# 
# ## 5.4. get canned calibration data and calibrate rainwater against feather d2H  ----
# 
# # data("knownOrig")
# # knownOrig$sources[1:2]
# # knownOrig$samples %>% filter(Group=="Ground bird")
# # 
# # cal.grd.bird <- subOrigData(marker = 'd2H', dataset = c(5), ref_scale = NULL) # Extract data
# # str(cal.grd.bird$data)
# # 
# # d2Hcalib <- calRaster(known = cal.grd.bird, isoscape = isoscape, interpMethod = 1, verboseLM = F, genplot = F) # Calibrate the mean and sd values from our isoscape
# # summary(d2Hcalib$lm.model) # Extract model results
# #  
# # ggplot(data = d2Hcalib$lm.data, aes(y = tissue.iso, x = isoscape.iso)) + 
# #   geom_point()  + 
# #   stat_smooth(method = "lm", formula = 'y ~ x') + 
# #   theme_classic()
# 
# ## manually rescale the isoscape with the Powell equation
# 
# isoscape.AD<- 4.50799 + 0.79760*isoscape
# isoscape.JUV<- 4.50799 + 0.79760*isoscape - 28.73077
# 
# plot(isoscape.AD, xlab="Longitude", ylab="Latitude")
# 
# 
# 
# 
# ## 5.4.1 building our own equation using the known origin feather samples ----
# ## this looks weird like somebody had calculated the dH_reg
# # summary(lm(dH~dH_reg+AGE, data=ORIG_WC))
# 
# ## 5.4.1.1. convert data to SpatVector ------
# 
# woco.vect<-terra::vect(woco.sf)
# 
# 
# ## 5.4.1.2. extract rainwater hydrogen isotopes for the location of known-sample woodcocks -------
# 
# woco.sf$d2h_GS<-terra::extract(isoscape,woco.vect)$d2h
# woco.sf$d2h_se_GS<-terra::extract(isoscape,woco.vect)$d2h.se
# 
# 
# ## 5.4.1.3. create equation to link rainwater to d2H of feather -------
# ggplot(woco.sf, aes(x=d2h_GS,y=dH_reg, col=AGE, fill=AGE)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   labs(x="d2H in rainwater (growing season average) at sampling location", y="d2H in Swiss woodcock feather") +
#   ## beautification of the axes
#   theme(panel.background=element_rect(fill="white", colour="black"),
#         panel.grid.major = element_line(linewidth=0.4, colour="grey89", linetype="dashed"),
#         panel.grid.minor = element_blank(),
#         axis.text=element_text(size=16, color="black"), 
#         axis.title=element_text(size=18),
#         legend.position="inside",
#         legend.position.inside=c(0.10,0.90),
#         legend.box.background = element_blank(),
#         legend.background = element_blank(),
#         legend.title=element_text(size=18),
#         legend.text=element_text(size=16, color="black"))
# 
# ggsave("output/SUI_WOCO_feather_isotope_calibration.jpg", width=9, height=8)
# 
# ISO_CALIB<-lm(dH_reg~d2h_GS, data=woco.sf)
# summary(ISO_CALIB)
# 
# 
# 
# 
# ## 5.4.1.4. apply that equation to all other shot woodcocks -------
# woco.unk.sf <- UNK_WC %>%
#   st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
# woco.unk.vect<-terra::vect(woco.unk.sf)
# 
# woco.unk.sf$d2h_GS<-terra::extract(isoscape,woco.unk.vect)$d2h
# woco.unk.sf$d2h_se_GS<-terra::extract(isoscape,woco.unk.vect)$d2h.se
# woco.unk.sf %>% filter(is.na(AGE))
# 
# d2H_pred<-predict(ISO_CALIB,newdata=woco.unk.sf,se.fit=TRUE,interval="prediction",level=0.95, type = "response")$fit
# 
# 
# 
# ## 5.4.1.5. test whether each shot feather came from within location-specific prediction interval -------
# 
# OUT_SUMMARY<-bind_cols(woco.unk.sf,d2H_pred) %>%
#   filter(!is.na(dH_reg)) %>%
#   st_drop_geometry() %>%
#   mutate(LOCAL=ifelse(dH_reg>lwr & dH_reg<upr,1,0)) %>%
#   group_by(AGE, KANTON) %>%
#   summarise(prop.local=mean(LOCAL, na.rm=T)) %>%
#   spread(key=AGE, value=prop.local)
# 
# fwrite(OUT_SUMMARY, "output/prop_local_WOCO_shot_isotopes.csv")
# 
# 
# 
# 
# 
# ## 5.5. PRIOR OF DISTRIBUTION  ----
# min(values(isoscape)[,1],na.rm=T)
# max(values(isoscape)[,1],na.rm=T)
# distd2H_knownSUI<-hist(ORIG_WC$dH_correct, breaks=seq(-135,-5,5))
# prior.probs<-tibble(d2H=distd2H_knownSUI$breaks[-1], prob=distd2H_knownSUI$density)
# prior.probs$prob[findInterval(x=woco$dH_correct, vec=prior.probs$d2H)]
# 
# priorscape <- isoscape
# length(values(isoscape)[,1])
# length(values(priorscape)[,1])
# length(findInterval(x=values(isoscape)[,1], vec=prior.probs$d2H))
# values(priorscape)[,1]<-prior.probs$prob[findInterval(x=values(isoscape)[,1], vec=prior.probs$d2H)]
# values(priorscape)[,2]<-2  ## arbitrary value for se of prior probability
# 
# 
# 
# 
# ## 5.6. ASSIGNMENT TO ORIGIN  ----
# ## this is computationally very intensive!!
# ### because the isotope data are already converted into rainfall isotopes, we do not need the calibrated isoscape
# origins <- pdRaster(isoscape,
#                     unknown=UNK_WC[!(is.na(UNK_WC$dH_correct)),c(1,9)],
#                     #prior=priorscape,
#                     mask = as(EUR, 'Spatial'),
#                     genplot = F)
#                     #outDir="C:/Users/sop/OneDrive - Vogelwarte/Woodcock/output")
# 
# 
# saveRDS(origins,"output/WOCO_origin_assignment.rds")
# origins<-readRDS("output/WOCO_origin_assignment.rds")
# 
# ### check whether values sum to 1 - they do NOT sum to 1!!??
# ALLout<-values(origins)
# ALLprob<-apply(ALLout,2,sum,na.rm=T)
# summary(ALLprob)
# global(origins[[1]], 'sum', na.rm = TRUE)
# 
# ## 4.7. CROP ASSIGNMENT TO SWITZERLAND  ----
# 
# SUIorigins<-terra::crop(origins,SUI)
# dim(SUIorigins)
# out<-values(SUIorigins)
# 
# 
# ## extract the cumulative probability of Switzerland for all of the birds shot in Switzerland
# SUIprob<-apply(out,2,sum)
# UNK_WC$Swissorigin_prob<-as.numeric(SUIprob)
# any(names(SUIprob)!=UNK_WC$ID)
# 
# ## plot output against isotopes and date
# ggplot(UNK_WC, aes(x=Swissorigin_prob, col=PROVENANCE_voigt, fill=PROVENANCE_voigt)) +
#   geom_histogram(alpha=0.5,position = position_dodge(width=1))
# 
# 
# ggplot(UNK_WC, aes(x=Swissorigin_prob,y=dH_correct, col=AGE)) +
#   geom_point()
# 


############# TEMPLATE CODE PROVIDED BY COPILOT FOR COMBINING BOTH TIME AND d2H IN ASSIGNMENT PROBABILITY ##########

library(nimble)

# Define the model
code <- nimbleCode({
  # Prior for class membership
  z ~ dcat(pi[1:2])  # z = 1 for A, z = 2 for B
  
  # Priors for class probabilities
  pi[1] <- 0.5
  pi[2] <- 0.5
  
  # Conditional distributions
  d2H ~ dnorm(mu_d2H[z], sd = sigma_d2H[z])
  time ~ dcustom(time_likelihood(time, z))
})

# Custom likelihood for time variable
time_likelihood <- nimbleFunction(
  run = function(x = double(0), z = integer(0)) {
    returnType(double(0))
    if (z == 1) {
      # Logistic distribution for Population A
      mu <- 0
      s <- 1
      loglik <- -x + 2 * log(1 / (1 + exp(-x)))
    } else {
      # Inverse distribution for Population B (e.g., 1/x)
      if (x <= 0) return(-1e6)  # Penalize invalid values
      loglik <- -log(x^2)
    }
    return(loglik)
  }
)

# Constants and data
constants <- list(
  mu_d2H = c(-50, -100),
  sigma_d2H = c(10, 15)
)

# Example data point
data <- list(
  d2H = -60,
  time = 2
)

# Initial values
inits <- list(z = 1)

# Build and compile the model
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

# Configure and run MCMC
conf <- configureMCMC(model, monitors = "z")
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)
samples <- runMCMC(cmcmc, niter = 1000)

# Posterior probability of class membership
table(samples) / length(samples)

