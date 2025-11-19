#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### ORIGIN ASSIGNMENT OF WOODCOCKS SHOT IN SWITZERLAND ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025
## includes feather isotope ratios and time of harvest to inform probability whether woodcock was of local (Swiss) origin
## completely revised on 17 Nov 2025 to include 3 models with different priors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

## minor adjustments based on Pius Korner's suggestions from 25 July 2025
## 3 MODELS INCLUDE (1) no priors; (2) both timing and abundance prior; (3) only abundance prior
## this leads to vastly different estimates, and it will be necessary to report that difference

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



# 2. SPECIFY COMBINED PROBABILITY MODEL TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------

woco.t.abd.model<-nimbleCode({
  
  # Parameters:
  # p.nonlocal: probability of non-local origin (that shot woodcock was a not a local bird) - varies by individual
  
  # b.age: regression parameter for age that translates rainfall d2H to feather d2H 
  # b.rain: regression parameter for rainfall d2H that translates rainfall d2H to feather d2H
  # int.rain: regression intercept that translates rainfall d2H to feather d2H 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CALIBRATION REGRESSION FOR KNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Likelihood:  
  for (i in 1:nind.known){
        d2H_feather.known[i] ~ dnorm(mu.known[i], sd=d2H_rain.known.sd[i])
        mu.known[i] <- int.rain + b.age*age.known[i] + b.rain*d2H_rain.known[i]
      }
    
  # Priors informed by known origin woodcock feathers in UK see code above
  int.rain ~ dnorm(-13, sd=56) # 
  b.age ~ dnorm(-24, sd=35) # 
  b.rain ~ dnorm(0.2, sd=1) # 
  dispersion ~ dnorm(25,sd=1) # dispersion parameter to convert prior probability into beta distribution - almost fixed quantity
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (i in 1:nind.unkn){
    
    # Potential local distribution based on timing used as prior (from WOCO_migration_model.R - could be integrated into this model)
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal.prior1[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from migration model
    p.nonlocal.prior[i] <- max(p.nonlocal.prior1[i],p.nonlocal.prior2[i])  ## Option 3: combining priors by taking the minimum (once birds have migrated the lowest prob will do)
    
    # Latent indicator for isotope distribution based on beta distribution of shooting time prior
    alpha[i] <- max(1e-3, p.nonlocal.prior[i] * dispersion)   ## with safeguard to avoid values <0
    beta[i]  <- max(1e-3, (1 - p.nonlocal.prior[i]) * dispersion)   ## with safeguard to avoid values <0
    p.nonlocal[i] ~ dbeta(alpha[i], beta[i])
    z[i] ~ dbern(p.nonlocal[i])  ## indicator variable from binomial draw of which isotope distribution fits better
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.Europe[year[i]]    ### overall year-specific distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.local[i]  ### local expected feather isotope distribution
    
    sd.unknown[i,2] <- sd.rain.d2H.Europe[year[i]]    ### overall year-specific variance across Europe
    sd.unknown[i,1] <- sd.rain.d2H.local[i]  ### local expected feather isotope distribution variation
    
    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[i,z[i]+1])
    
  } #i
  

}) ## end of nimble code chunk






# 3. prepare the data needed for NIMBLE input ------------------------------------------------------------
woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_MA)) %>%
  mutate(year=year(feather_growth_date))

woco.sf <- woco.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_MA)) %>%
  mutate(year=year(feather_growth_date))

table(woco.unk.sf$AGE,woco.unk.sf$KANTON)

## 3.1. DISTINGUISH CONSTANTS AND DATA-----------------
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      nind.known = dim(woco.sf)[1],
                      mean.rain.d2H.Europe = mean.rain.d2H.Europe,
                      mean.rain.d2H.local=woco.unk.sf$d2h_MA,
                      sd.rain.d2H.Europe = sd.rain.d2H.Europe,
                      sd.rain.d2H.local = ifelse(sqrt(woco.unk.sf$d2h_se_MA)==0,mean(woco.unk.sf$d2h_MA),sqrt(woco.unk.sf$d2h_se_MA)),  ## given as variance in data frame, needs to be sqrt for sd
                      d2H_rain.known=woco.sf$d2h_MA,
                      d2H_rain.known.sd=sqrt(woco.sf$d2h_se_MA),   ## given as variance in data frame, needs to be sqrt for sd
                      age.unknown = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      age.known = ifelse(woco.sf$AGE=="Adulte",0,1),
                      p.nonlocal.prior2 = woco.unk.sf$abd_prior,
                      year=as.numeric(as.factor(woco.unk.sf$year)),
                      # lm.mean.mig=logit(readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[2,1]),
                      # b.mig.week=readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[1,1],
                      # lm.mean.mig=logit(MCMCsummary(readRDS("output/woco_mig_depart_output_nimble.rds"))[2,1]),  ## became necessary after switching to stepwise run which produces different output than nimbleMCMC
                      # b.mig.week=MCMCsummary(readRDS("output/woco_mig_depart_output_nimble.rds"))[1,1],  ## became necessary after switching to stepwise run which produces different output than nimbleMCMC
                      lm.mean.mig=logit(as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[1,1])),
                      b.mig.week=as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[2,1]),
                      unk.week=week(UNK_WC$feather_sampling_date)-week(ymd("2024-07-26"))   ## migration weeks start only in August
                      )

iso.data <- list(d2H_feather.known = woco.sf$dH,
                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.2. specify NIMBLE run settings --------------------

# Parameters monitored
parameters.iso <- c("b.rain","b.age","int.rain","dispersion","p.nonlocal") #,"p.nonlocal.prior","p.nonlocal.prior1") including these bloats the output object

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8
iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_MA-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,-13, sd=20), # informative prior based on Powell
                  b.age = rnorm(1,-24, sd=35), # informative prior based on Powell
                  b.rain = rnorm(1,0.2, sd=1), # informative prior based on Powell
                  dispersion = 25, # dispersion parameter for beta distribution set to sensible value
                  p.nonlocal = 1-woco.unk.sf$abd_prior # initial start value will be replaced after test run
)


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 100000
n.burnin <- 50000
n.chains <- 4



## 3.3. PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.t.abd.model,
                    constants=iso.constants,
                    data = iso.data,
                    inits = iso.inits,
                    calculate=TRUE)
# cModel <- compileNimble(test)
# 
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
# test$initializeInfo()
# #help(modelInitialization)
# test$p.nonlocal.prior1
# test$p.nonlocal.prior
# summary(test$p.nonlocal)
# test$b.age
# test$b.rain
# test$int.rain
# 
# hist(rbeta(1000,test$alpha[17],test$beta[17]))
# 
# # Calculate all log probabilities
# cModel$calculate()
# 
# # Get all stochastic nodes
# nodes <- cModel$getNodeNames(stochOnly = TRUE)
# 
# # Find nodes with -Inf logProb
# badNodes <- nodes[sapply(nodes, function(n) cModel$getLogProb(n)) == -Inf]
# cat("Nodes with -Inf logProb:\n")
# print(badNodes)
# 
# cModel$logProb_d2H_feather.unknown[17]  ### why are there no -Inf?
# iso.data$d2H_feather.unknown[17]
# cModel$logit.mig.unknown[16:17]
# cModel$p.nonlocal[16:17]
# cModel$mu.unknown[16:17,1]
# cModel$mu.unknown[16:17,2]
# cModel$sd.unknown[16:17,1]
# cModel$sd.unknown[16:17,2]
  


### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits

## update initial values
iso.inits$p.nonlocal<-test$p.nonlocal.prior


## 3.4. run model in NIMBLE ------------------------------------------
### this takes 600 sec for 100000 iterations and converges in that time
tic()

woco.iso <- nimbleMCMC(code = woco.t.abd.model,
                        constants=iso.constants,
                        data = iso.data,
                        inits = iso.inits,,
                        monitors = parameters.iso,
                        thin=5,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains,
                        progressBar = getNimbleOption("MCMCprogressBar"),
                        summary=T)
toc() 



# saveRDS(woco.iso,"output/woco_iso_time_origin_model_comb_prior.rds") # no point saving - it only runs 10 min and the file is too large for GitHub
# woco.iso<-readRDS("output/woco_iso_time_origin_model_comb_prior.rds")





# 4. SAVE MODEL OUTPUT AND ASSESS CONVERGENCE --------------------------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("b.rain","b.age","int.rain","dispersion"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
#fwrite(out,"output/woco_iso_time_origin_parm_estimates_comb_prior.csv")
#out<-fread("output/woco_iso_time_origin_parm_estimates_comb_prior.csv")


## for comparing the p.nonlocal estimates for 842 shot birds
# out<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("p.nonlocal"))) #"int.abd","b.countday","b2.countday","r.abd")))
# out$parameter<-row.names(out)
# names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
# out
# #fwrite(out,"output/woco_p_nonlocal_comb_prior_SUI_calib.csv")
# comp<-fread("output/woco_p_nonlocal_comb_prior_SUI_calib.csv") %>%
#   dplyr::select(parameter,median) %>%
#   left_join(out, by="parameter") %>%
#   mutate(diff=median.x-median.y)
# hist(comp$diff)
# summary(comp)
# 
# hist(rnorm(1000,mean=(5.53-31.3*1+1.14*mean.rain.d2H), sd=(13.19)))


#MCMCplot(woco.iso$samples, params=c("p.nonlocal"))


## look at the chains and whether they mixed well
chainsPlot(woco.iso$samples,
           var=c("b.rain","b.age","int.rain","dispersion"))



## 4.1 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples <- rbind(woco.iso$samples$chain1,woco.iso$samples$chain2,woco.iso$samples$chain3,woco.iso$samples$chain4)
mean.p.nonlocal <- as_tibble(samples[,grep("p.nonlocal\\[", colnames(samples))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) %>%
  mutate(prior="combined abundance and migration")
##fwrite(mean.p.nonlocal,"output/WOCO_nonlocal_probs_comb_prior.csv")

mean.p.nonlocal %>% filter(is.na(age))

# summarise across SUI

out.sui<- mean.p.nonlocal %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.sui  
#fwrite(out.sui,"output/woco_nonlocal_origin_estimates_SUI_comb_prior.csv")


# summarise by Canton

out.ctn<- mean.p.nonlocal %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age,ctn) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.ctn
#fwrite(out.ctn,"output/woco_nonlocal_origin_estimates_CANTON_comb_prior.csv")





# 5. RUN MIGRATION ONLY PROBABILITY MODEL TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------

woco.t.model<-nimbleCode({
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CALIBRATION REGRESSION FOR KNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Likelihood:  
  for (i in 1:nind.known){
    d2H_feather.known[i] ~ dnorm(mu.known[i], sd=d2H_rain.known.sd[i])
    mu.known[i] <- int.rain + b.age*age.known[i] + b.rain*d2H_rain.known[i]
  }
  
  # Priors informed by known origin woodcock feathers in UK see code above
  int.rain ~ dnorm(-13, sd=56) # 
  b.age ~ dnorm(-24, sd=35) # 
  b.rain ~ dnorm(0.2, sd=1) # 
  dispersion ~ dnorm(25,sd=1) # dispersion parameter to convert prior probability into beta distribution - almost fixed quantity
  

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (i in 1:nind.unkn){
    
    # Potential local distribution based on timing used as prior (from WOCO_migration_model.R - could be integrated into this model)
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal.prior[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from migration model

    # Latent indicator for isotope distribution based on beta distribution of shooting time prior
    alpha[i] <- max(1e-3, p.nonlocal.prior[i] * dispersion)   ## with safeguard to avoid values <0
    beta[i]  <- max(1e-3, (1 - p.nonlocal.prior[i]) * dispersion)   ## with safeguard to avoid values <0
    p.nonlocal[i] ~ dbeta(alpha[i], beta[i])
    z[i] ~ dbern(p.nonlocal[i])  ## indicator variable from binomial draw of which isotope distribution fits better
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.Europe[year[i]]    ### overall year-specific distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.local[i]  ### local expected feather isotope distribution
    
    sd.unknown[i,2] <- sd.rain.d2H.Europe[year[i]]    ### overall year-specific variance across Europe
    sd.unknown[i,1] <- sd.rain.d2H.local[i]  ### local expected feather isotope distribution variation
    
    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[i,z[i]+1])
    
  } #i
  
  
}) ## end of nimble code chunk




## 5.2. run model in NIMBLE ------------------------------------------
### this takes 600 sec for 100000 iterations and converges in that time
tic()

woco.iso.migprior <- nimbleMCMC(code = woco.t.model,
                       constants=iso.constants,
                       data = iso.data,
                       inits = iso.inits,
                       monitors = parameters.iso,
                       thin=5,
                       niter = n.iter,
                       nburnin = n.burnin,
                       nchains = n.chains,
                       progressBar = getNimbleOption("MCMCprogressBar"),
                       summary=T)
toc() 





## 5.3. SAVE MODEL OUTPUT AND ASSESS CONVERGENCE --------------------------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso.migprior$samples, params=c("b.rain","b.age","int.rain","sigma.calib","dispersion"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_iso_time_origin_parm_estimates_mig_prior.csv")



## 5.4 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples.migprior <- rbind(woco.iso.migprior$samples$chain1,woco.iso.migprior$samples$chain2,woco.iso.migprior$samples$chain3,woco.iso.migprior$samples$chain4)
mean.p.nonlocal.migprior <- as_tibble(samples.migprior[,grep("p.nonlocal\\[", colnames(samples.migprior))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) %>%
  mutate(prior="only migration")
fwrite(mean.p.nonlocal.migprior,"output/WOCO_nonlocal_probs_mig_prior.csv")

# summarise across SUI

out.sui<- mean.p.nonlocal.migprior %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  select(-age)  %>%
  mutate(prior="only migration") %>%
  bind_rows(out.sui)
out.sui  
fwrite(out.sui,"output/woco_nonlocal_origin_estimates_SUI.csv")


# summarise by Canton

out.ctn<- mean.p.nonlocal.migprior %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age,ctn) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age)  %>%
  mutate(prior="only migration") %>%
  bind_rows(out.ctn)
out.ctn
fwrite(out.ctn,"output/woco_nonlocal_origin_estimates_CANTON.csv")








# 6. RUN ASSIGNMENT MODEL WITHOUT PRIOR TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------

woco.null.model<-nimbleCode({
  
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

  # Standard deviation for isotope ratios in rainwater 
  sd.unknown[2]<-max(sigma.calib,sd.rain.d2H)  ### overall distribution across Europe
  sd.unknown[1]<-sigma.calib
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SPECIFY AGE AND CANTON SPECIFIC PRIORS
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (ag in c(1,2)){
    for(ctn in 1:ncantons){
      p.nonlocal[ag,ctn] ~dbeta(2,2)
    }
  }
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (i in 1:nind.unkn){
    
    # Latent indicator for isotope distribution based on beta distribution of shooting time prior
    z[i] ~ dbern(p.nonlocal[age[i],canton[i]])  ## indicator variable from binomial draw of which isotope distribution fits better
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.Europe[year[i]]    ### overall year-specific distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H.local[i]  ### local expected feather isotope distribution
    
    sd.unknown[i,2] <- sd.rain.d2H.Europe[year[i]]    ### overall year-specific variance across Europe
    sd.unknown[i,1] <- sd.rain.d2H.local[i]  ### local expected feather isotope distribution variation
    
    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[i,z[i]+1])
    
    
  } #i
  
  
}) ## end of nimble code chunk



## 6.1. DISTINGUISH CONSTANTS AND DATA-----------------
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      nind.known = dim(woco.sf)[1],
                      mean.rain.d2H.Europe = mean.rain.d2H,
                      sd.rain.d2H.Europe = sd.rain.d2H,
                      sd.rain.d2H.local = sd.rain.d2H,
                      d2H_rain.known=woco.sf$d2h_MA,
                      mean.rain.d2H.local=woco.unk.sf$d2h_MA,
                      age.unknown = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      age.known = ifelse(woco.sf$AGE=="Adulte",0,1),
                      canton=as.numeric(as.factor(woco.unk.sf$CANTON)),
                      year=as.numeric(as.factor(woco.unk.sf$year)),
                      ncantons=length(unique(woco.unk.sf$CANTON)),
                      lm.mean.mig=logit(as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[1,1])),
                      b.mig.week=as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[2,1]),
                      unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26"))   ## migration weeks start only in August
)

iso.data <- list(d2H_feather.known = woco.sf$dH,
                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.2. specify NIMBLE run settings --------------------

# Parameters monitored
parameters.iso <- c("b.rain","b.age","int.rain","p.nonlocal") #,"p.nonlocal.prior","p.nonlocal.prior1") including these bloats the output object

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_MA-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,4.5, sd=5), # informative prior based on Powell
                  b.age = rnorm(1,-28.7, sd=5), # informative prior based on Powell
                  b.rain = rnorm(1,0.8, sd=1), # informative prior based on Powell
                  p.nonlocal = matrix(rbeta(16,2,2),ncol=8,nrow=2)
)

## 6.2. run model in NIMBLE ------------------------------------------
### this takes 600 sec for 100000 iterations and converges in that time
tic()

woco.iso.null <- nimbleMCMC(code = woco.null.model,
                                constants=iso.constants,
                                data = iso.data,
                                inits = iso.inits,
                                monitors = parameters.iso,
                                thin=5,
                                niter = n.iter,
                                nburnin = n.burnin,
                                nchains = n.chains,
                                progressBar = getNimbleOption("MCMCprogressBar"),
                                summary=T)
toc() 





## 5.3. SAVE MODEL OUTPUT AND ASSESS CONVERGENCE --------------------------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso.migprior$samples, params=c("b.rain","b.age","int.rain","sigma.calib","dispersion"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_iso_time_origin_parm_estimates_mig_prior.csv")



## 5.4 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples.migprior <- rbind(woco.iso.migprior$samples$chain1,woco.iso.migprior$samples$chain2,woco.iso.migprior$samples$chain3,woco.iso.migprior$samples$chain4)
mean.p.nonlocal.migprior <- as_tibble(samples.migprior[,grep("p.nonlocal\\[", colnames(samples.migprior))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) %>%
  mutate(prior="only migration")
fwrite(mean.p.nonlocal.migprior,"output/WOCO_nonlocal_probs_mig_prior.csv")

# summarise across SUI

out.sui<- mean.p.nonlocal.migprior %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  select(-age)  %>%
  mutate(prior="only migration") %>%
  bind_rows(out.sui)
out.sui  
fwrite(out.sui,"output/woco_nonlocal_origin_estimates_SUI.csv")


# summarise by Canton

out.ctn<- mean.p.nonlocal.migprior %>%
  #group_by(age,ctn,ind) %>%
  #summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  #ungroup() %>%
  group_by(age,ctn) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age)  %>%
  mutate(prior="only migration") %>%
  bind_rows(out.ctn)
out.ctn
fwrite(out.ctn,"output/woco_nonlocal_origin_estimates_CANTON.csv")









# 7. summarise output in graphical form ------------------------------------------

# mean.p.nonlocal.migprior<- fread("output/WOCO_nonlocal_probs_mig_prior.csv")
# mean.p.nonlocal<- fread("output/WOCO_nonlocal_probs_comb_prior.csv")

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



FIGURE2<-bind_rows(mean.p.nonlocal,mean.p.nonlocal.migprior) %>%
  group_by(age,ctn,ind, prior) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(age,ctn, prior) %>%
  summarise(for.med=median(p.nonlocal.mean),for.ucl=quantile(p.nonlocal.mean,0.025), for.lcl=quantile(p.nonlocal.mean,0.975)) %>%
  
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  #mutate(Kanton=levels(as.factor(woco.unk.sf$KANTON))[ctn]) %>%
  
  ggplot(aes(x=ctn, y=for.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=for.lcl, ymax=for.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  facet_wrap(~prior, ncol = 1) +
  
  # annotation_custom(grob=gunicon, xmin=0.5, xmax=1.5, ymin=0.05, ymax=0.18) +
  # annotation_custom(wocoicon, xmin=0.5, xmax=2.9, ymin=0.10, ymax=0.35) +
  
  ## format axis ticks
  labs(y="Proportion of shot woodcocks of non-local origin",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +

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
        legend.position.inside=c(0.90,0.65),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE2


ggsave(plot=FIGURE2,
       filename="output/woco_iso_time_origin_probability_estimates_comb_prior.jpg", 
       device="jpg",width=9, height=12)






# 6. SUMMARY FIGURES FOR MANUSCRIPT ----------------------------------

## 6.1. PLOT HISTOGRAMS FOR SUISSE AND OTHER BIRDS ----
woco$ORIGINE<-ifelse(woco$ORIGINE=="SCHWEIZ","Switzerland","unknown")
FIGURES2<-ggplot(woco, aes(x=dH, col=ORIGINE, fill=ORIGINE)) +
  geom_histogram(alpha=0.5,position = position_dodge(width=1)) +
  
  ## format axis ticks
  labs(y="Number of woodcock feathers",
       x=expression(paste(delta^{2}, "H (\u2030)")),
       col="Origin", fill="Origin") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_line(linewidth=0.5, colour="grey59", linetype="dashed"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14, color="black"),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title=element_text(size=14, color="black"),
        legend.position="inside",
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.position.inside=c(0.85,0.85),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURES2
#ggsave("output/WOCO_isotope_histogram_by_origin.jpg")


## report numbers in manuscript
table(ORIG_WC$AGE)
length(UNK_WC$dH)/9681
table(UNK_WC$AGE)
summary(ORIG_WC$dH[ORIG_WC$PROVENANCE_voigt=="Suisse"])
summary(UNK_WC$dH)














# 7. CREATE MAPS FOR EACH CANTON WHAT LOCAL RAINFALL ENCOMPASSES ----------------------------------
plot_list <- list()
WOCO.isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
  crop(woco.countries)
for (ct in 1:length(unique(woco.unk.sf$KANTON))){

  ## get canton-wise distribution
  woco.cnt <- woco.unk.sf %>% dplyr::filter(KANTON==unique(woco.unk.sf$KANTON)[ct]) 
  cnt.iso <-  woco.cnt %>% st_drop_geometry() %>%
    dplyr::select(ID,KANTON,d2h_MA,d2h_se_MA) %>%
    rowwise() %>%
    mutate(pot.orig.d2H=rnorm(1,d2h_MA,d2h_se_MA)) %>%
    ungroup() %>%
    group_by(KANTON) %>%
    summarise(min=min(pot.orig.d2H), max=max(pot.orig.d2H))
  
  ## EXTRACT ISOSCAPE CELLS FALLING INTO THIS RANGE
  CNT.isoscape <- ifel(WOCO.isoscape >= cnt.iso$min & WOCO.isoscape <= cnt.iso$max, 1, 0)
  
  ## create plot
  plot_list[[ct]] <- ggplot(EUR) +
    geom_sf() +
    tidyterra::geom_spatraster(data=CNT.isoscape, aes(fill=d2h_MA))+
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

