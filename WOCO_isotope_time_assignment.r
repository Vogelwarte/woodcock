### ISOTOPic ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025




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
load("data/woco.input.data.RData")

# 2. COMBINED PROBABILITY MODEL TO ESTIMATE PROB OF LOCAL ORIGIN IN NIMBLE ----

## 3.1. specify model in NIMBLE ----

woco.orig.model<-nimbleCode({
  
  # Parameters:
  # p.nonlocal: probability of local origin (that shot woodcock was a local bird) - varies by age and canton
  
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

  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # Likelihood 
  sd.unknown[2]<-sigma.calib+sd.rain.d2H  ### overall distribution across Europe
  sd.unknown[1]<-sigma.calib
  
  for (i in 1:nind.unkn){
    
    # Potential local distribution based on timing used as prior
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from departure model
    
    # Latent indicator for isotope distribution
    z[i] ~ dbern(p.nonlocal[i])  ## 
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H    ### overall distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*d2H_rain.unknown[i]  ### local expected feather isotope distribution

    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[z[i]+1])
    
    
  } #i
  

}) ## end of nimble code chunk






## 3.2. prepare the data needed for NIMBLE input ----
woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
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
                      lm.mean.mig=logit(readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[2,1]),
                      b.mig.week=readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[1,1],
                      unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26"))   ## migration weeks start only in August
                      )

iso.data <- list(d2H_feather.known = woco.unk.sf$dH,

                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.3. specify NIMBLE run settings ----

# Parameters monitored
parameters.iso <- c("b.rain","b.age","int.rain","sigma.calib","p.nonlocal")

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_GS-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,4.5, sd=5), # informative prior based on Powell
                  b.age = rnorm(1,-28.7, sd=5), # informative prior based on Powell
                  b.rain = rnorm(1,0.8, sd=1), # informative prior based on Powell
                  sigma.calib = runif(1,0,20) # standard deviation

)


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



saveRDS(woco.iso,"output/woco_iso_time_origin_model.rds")
#woco_surv<-readRDS("output/woco_iso_origin_model.rds")





## 3.5. assess model and examine convergence ------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("p.nonlocal","b.rain","b.age","int.rain","sigma.calib"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_iso_time_origin_parm_estimates.csv")
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

grid.arrange(grobs=plot_list,ncol=2)
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

