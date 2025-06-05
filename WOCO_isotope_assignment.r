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
library(nimble)
library(tictoc)
filter<-dplyr::filter
select<-dplyr::select


## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")



# 1. READ IN PROCESSED ISOTOPE DATA ----------------

### data are poorly documented - my reconstruction is as follows based on Voigt 2016 report:
### dH: raw isotope measurements in feather keratin (but only available for half of the samples!!) - based on Berlin lab standard
### dH_reg: corresponding rainfall hydrogen isotope ratio converted with equation in Powell (2012) (age-specific adjustment)
### dH_scaled: is the converted d2H of feather keratin using the formula 0.979*dH + 20.701 to adjust it to Canada lab standard
### dH_correct: is the corrected d2H of feathers IF a value was in column 'bc' (whatever the fuck this stands for)


woco<-fread("data/WOCO_isotopes.csv")
#woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")


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
ggsave("output/WOCO_isotope_histogram_by_origin.jpg")


# distd2H_knownSUI<-hist(ORIG_WC$dH, breaks=seq(-130,30,5))
# distd2H_knownSUI$density   ## this cannot be used as prior information for the isotope assignment because it is just from a corner of Switzerland




## 1.4. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
  dplyr::filter(!is.na(dH)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE)
UNK_WC<-woco %>% filter(ORIGINE=="UNBEKANNT") %>%
  dplyr::filter(!is.na(dH)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE)
UNK_WC$ID<-str_replace_all(UNK_WC$ID, "[^[:alnum:]]", " ")
dim(ORIG_WC)
dim(UNK_WC)









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


ggplot(SUI) +
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



## overlay isoscape with woodcock distribution and generate mean feather distribution
woco.countries <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Switzerland","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))
WOCO.isoscape <- readRDS("data/global_d2H_GS_isoscape.rds") %>%
  crop(woco.countries)
plot(WOCO.isoscape)

mean.rain.d2H<-mean(na.omit(values(WOCO.isoscape)[,1]), na.rm=T)
sd.rain.d2H<-sd(na.omit(values(WOCO.isoscape)[,1]), na.rm=T)
hist(rnorm(1000,mean.rain.d2H,sd.rain.d2H))


## 2.4. use known origin data to calibrate rainwater against feather d2H  -------


### 2.4.1 convert woco data to SpatVector to extract values from raster ------

woco.vect<-terra::vect(woco.sf)


## 2.4.2. extract rainwater hydrogen isotopes for the location of known-sample woodcocks -------

woco.sf$d2h_GS<-terra::extract(isoscape,woco.vect)$d2h
woco.sf$d2h_se_GS<-terra::extract(isoscape,woco.vect)$d2h.se


## 2.4.3. create equation to link rainwater to d2H of feather -------
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
woco.unk.sf <- UNK_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
woco.unk.vect<-terra::vect(woco.unk.sf)

woco.unk.sf$d2h_GS<-terra::extract(isoscape,woco.unk.vect)$d2h
woco.unk.sf$d2h_se_GS<-terra::extract(isoscape,woco.unk.vect)$d2h.se
woco.unk.sf %>% filter(is.na(AGE))

d2H_pred<-predict(ISO_CALIB,newdata=woco.unk.sf,se.fit=TRUE,interval="prediction",level=0.95, type = "response")$fit



## 2.4.5. test whether each shot feather came from within location-specific prediction interval -------

OUT_SUMMARY<-bind_cols(woco.unk.sf,d2H_pred) %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  st_drop_geometry() %>%
  mutate(LOCAL=ifelse(dH>lwr & dH<upr,1,0)) %>%
  group_by(AGE, KANTON) %>%
  summarise(prop.local=mean(LOCAL, na.rm=T)) %>%
  spread(key=AGE, value=prop.local) %>%
  rename(prop_adult_local=Adulte, prop_juvenile_local=Juvénile)

OUT_SUMMARY<-bind_cols(woco.unk.sf,d2H_pred) %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  st_drop_geometry() %>%
  mutate(LOCAL=ifelse(dH>lwr & dH<upr,1,0)) %>%
  group_by(AGE, KANTON) %>%
  summarise(n=length(LOCAL)) %>%
  spread(key=AGE, value=n, fill=0) %>%
  rename(N_adult=Adulte, N_juvenile=Juvénile) %>%
  left_join(OUT_SUMMARY, by="KANTON")

OUT_SUMMARY

fwrite(OUT_SUMMARY, "output/prop_local_WOCO_shot_isotopes.csv")









# 3. COMBINED PROBABILITY MODEL TO ESTIMATE PROB OF LOCAL ORIGIN IN NIMBLE ----

## 3.1. specify model in NIMBLE ----

woco.orig.model<-nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # p.loc: probability of local origin (that shot woodcock was a local bird) - varies by age and canton
  
  # b.age: regression parameter for age that translates rainfall d2H to feather d2H 
  # b.rain: regression parameter for rainfall d2H that translates rainfall d2H to feather d2H
  # int.rain: regression intercept that translates rainfall d2H to feather d2H 
  
  # -------------------------------------------------
  # CALIBRATION REGRESSION FOR KNOWN ORIGIN BIRDS 
  # -------------------------------------------------
  # Likelihood:  
  for (i in 1:nind.known){
        d2H_feather.known[i] ~ dnorm(mu.known[i], sd=sigma.calib) # tau is precision (1 / variance)
        mu.known[i] <- int.rain + b.age*age.known[i] + b.rain*d2H_rain.known[i]
      }
    
  # Priors informed by known origin woodcock feathers in UK (Powell thesis 2012)
  int.rain ~ dnorm(4.5, sd=50) # informative prior based on Powell
  b.age ~ dnorm(-28.7, sd=25) # informative prior based on Powell
  b.rain ~ dnorm(0.8, sd=1) # informative prior based on Powell
  sigma.calib ~ dunif(0,10) # standard deviation
      
      


# -------------------------------------------------
# PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
# -------------------------------------------------

#### PRIOR PROBABILITIES OF LOCAL ORIGIN
for (ct in 1:ncanton){
  for (a in 1:2){
    p.loc[ct,a] ~ dbeta(2,4.5)   ##  hist(rbeta(1000,2,4.5)) very low and very high probabilities are less likely a priori
  } #a
} #ct

  # Likelihood 
  sd.unknown[1]<-sigma.calib+sd.rain.d2H  ### overall distribution across Europe
  sd.unknown[2]<-sigma.calib
  
  for (i in 1:nind.unkn){
    
    # Latent indicator 
    z[i] ~ dbern(p.loc[canton[i],age.unknown[i]+1])  ## age is specified as 0 and 1 in data, need to add 1 for index to work
    
    # Potential local distribution
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H    ### overall distribution across Europe
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*d2H_rain.unknown[i]

    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[z[i]+1])
  } #i
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
                      age.unknown = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      age.known = ifelse(woco.sf$AGE=="Adulte",0,1),
                      canton = as.numeric(as.factor(woco.unk.sf$KANTON)),
                      ncanton=length(unique(woco.unk.sf$KANTON))
                      )

iso.data <- list(d2H_feather.known = woco.unk.sf$dH,
                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.3. specify NIMBLE run settings ----

# Parameters monitored
parameters.iso <- c("p.loc","b.rain","b.age","int.rain","sigma.calib")

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_GS-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,4.5, sd=50), # informative prior based on Powell
                  b.age = rnorm(1,-28.7, sd=25), # informative prior based on Powell
                  b.rain = rnorm(1,0.8, sd=1), # informative prior based on Powell
                  sigma.calib = runif(1,0,10), # standard deviation
                  
                  #### FLAT PRIOR PROBABILITY FOR BEING LOCAL BIRD
                  p.loc=matrix(rbeta(length(unique(woco.unk.sf$KANTON))*2,2,4.5),ncol=2)  ###  hist(rbeta(1000,2,4.5))
 
)



# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.orig.model,
                    constants=iso.constants,
                    data = iso.data,
                    inits = iso.inits,
                    calculate=TRUE)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
test$initializeInfo()
#help(modelInitialization)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 50000
n.burnin <- 25000
n.chains <- 4



## 3.4. run model in NIMBLE ------------------------------------------



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
#woco_surv<-readRDS("output/woco_mig_depart_output_nimble.rds")






















# # 4. Geographic assignment using assignR ----
# 
# ## 4.1. load shapefile of Switzerland and EUROPE ----
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
# ## 4.2. plot the known origin points of WOCO that were sampled ----
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
# ## 4.3. extract isoscape ----
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
# ## 4.4. get canned calibration data and calibrate rainwater against feather d2H  ----
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
# ## 4.4.1 building our own equation using the known origin feather samples ----
# ## this looks weird like somebody had calculated the dH_reg
# # summary(lm(dH~dH_reg+AGE, data=ORIG_WC))
# 
# ## 4.4.1.1. convert data to SpatVector ------
# 
# woco.vect<-terra::vect(woco.sf)
# 
# 
# ## 4.4.1.2. extract rainwater hydrogen isotopes for the location of known-sample woodcocks -------
# 
# woco.sf$d2h_GS<-terra::extract(isoscape,woco.vect)$d2h
# woco.sf$d2h_se_GS<-terra::extract(isoscape,woco.vect)$d2h.se
# 
# 
# ## 4.4.1.3. create equation to link rainwater to d2H of feather -------
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
# ## 4.4.1.4. apply that equation to all other shot woodcocks -------
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
# ## 4.4.1.5. test whether each shot feather came from within location-specific prediction interval -------
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
# ## 4.5. PRIOR OF DISTRIBUTION  ----
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
# ## 4.6. ASSIGNMENT TO ORIGIN  ----
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

