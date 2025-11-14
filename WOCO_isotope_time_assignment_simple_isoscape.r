#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### ORIGIN ASSIGNMENT OF WOODCOCKS SHOT IN SWITZERLAND ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025
## includes feather isotope ratios and time of harvest to inform probability whether woodcock was of local (Swiss) origin
## completed on 13 Nov 2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

## adjusted with SIMPLE isoscape - fitted with feather isotope ratios directly, hence model needs no fractionation equation

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
# prepared in WOCO_feather_isoscape_creation.r
load("data/woco.input.data.simple.RData")


## 1.1. modify data to include abundance predictions from Ornitho.ch records ----------
## this approach is intended to address the dilution effect (when local birds are outnumbered by migrants)

woco.unk.abd.prior <- woco.unk.sf %>%
  st_drop_geometry() %>%
  mutate(DAY=yday(DATE), DAY_2=(yday(DATE)^2))

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


# 2. SPECIFY COMBINED PROBABILITY MODEL TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------
## NOTE: ensure you uncheck the desired prior option

woco.orig.model<-nimbleCode({
  
  # Parameters:
  # p.nonlocal: probability of non-local origin (that shot woodcock was a not a local bird) - varies by individual
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (i in 1:nind.unkn){
    
    # Potential local distribution based on timing used as prior (from WOCO_migration_model.R - could be integrated into this model)
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal.prior1[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from migration model
    p.nonlocal.prior[i] <- max(p.nonlocal.prior1[i],p.nonlocal.prior2[i])  ## combining priors by taking the maximum (once birds have migrated the highest prob will do)
    
    # Latent indicator for isotope distribution based on beta distribution of shooting time prior
    alpha[i] <- max(1e-3, p.nonlocal.prior[i] * dispersion)   ## with safeguard to avoid values <0
    beta[i]  <- max(1e-3, (1 - p.nonlocal.prior[i]) * dispersion)   ## with safeguard to avoid values <0
    p.nonlocal[i] ~ dbeta(alpha[i], beta[i])
    z[i] ~ dbern(p.nonlocal[i])  ## indicator variable from binomial draw of which isotope distribution fits better

    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- mean.feather.d2H[age[i]+1]    ### overall distribution across Europe
    mu.unknown[i,1] <- d2H_feather.local[i]  ### local expected feather isotope distribution
    
    sd.unknown[i,2] <- sd.feather.d2H[age[i]+1]    ### sd of overall distribution across Europe
    sd.unknown[i,1] <- d2H_feather.local.sd[i]  ### sd of local expected feather isotope distribution
    
    
    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i]+1], sd=sd.unknown[i,z[i]+1])
    
    
  } #i
  

}) ## end of nimble code chunk






# 3. prepare the data needed for NIMBLE input ------------------------------------------------------------

table(woco.unk.sf$AGE,woco.unk.sf$KANTON)

## 3.1. DISTINGUISH CONSTANTS AND DATA-----------------
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      mean.feather.d2H = c(mean.rain.d2H.ad,mean.rain.d2H.juv),
                      sd.feather.d2H = c(sd.rain.d2H.ad,sd.rain.d2H.juv),
                      d2H_feather.local=woco.unk.sf$d2h_local_predicted,
                      d2H_feather.local.sd=woco.unk.sf$d2h_local_sd,
                      age = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      p.nonlocal.prior2 = woco.unk.sf$abd_prior,
                      lm.mean.mig=logit(as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[1,1])),
                      b.mig.week=as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[2,1]),
                      unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26"))   ## migration weeks start only in August
                      )

iso.data <- list(d2H_feather.unknown = woco.unk.sf$dH)


## 3.2. specify NIMBLE run settings --------------------

# Parameters monitored
parameters.iso <- c("dispersion","p.nonlocal") #,"p.nonlocal.prior","p.nonlocal.prior1") including these bloats the output object

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < woco.unk.sf$d2h_local_predicted,0,1),
                  dispersion = 25, # dispersion parameter for beta distribution set to sensible value
                  p.nonlocal = 1-woco.unk.sf$abd_prior # initial start value will be replaced after test run
)


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 100000
n.burnin <- 50000
n.chains <- 4



## 3.3. PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.orig.model,
                    constants=iso.constants,
                    data = iso.data,
                    inits = iso.inits,
                    calculate=TRUE)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
test$initializeInfo()
#help(modelInitialization)
test$p.nonlocal.prior1
test$p.nonlocal.prior
test$p.nonlocal
hist(rbeta(1000,test$alpha[17],test$beta[17]))

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits

## update initial values
iso.inits$p.nonlocal<-test$p.nonlocal.prior


## 3.4. run model in NIMBLE ------------------------------------------
### this takes 600 sec for 100000 iterations and converges in that time
tic()

woco.iso <- nimbleMCMC(code = woco.orig.model,
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





# 4. SAVE MODEL OUTPUT AND ASSESS CONVERGENCE --------------------------------------------------------------


## 4.1 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples <- rbind(woco.iso$samples$chain1,woco.iso$samples$chain2,woco.iso$samples$chain3,woco.iso$samples$chain4)
mean.p.nonlocal <- as_tibble(samples[,grep("p.nonlocal\\[", colnames(samples))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) %>%
  mutate(prior="combined abundance and migration")


# summarise across SUI

out.sui<- mean.p.nonlocal %>%
  group_by(age) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.sui  
#fwrite(out.sui,"output/woco_nonlocal_origin_estimates_SUI_comb_prior.csv")


# summarise by Canton

out.ctn<- mean.p.nonlocal %>%
  group_by(age,ctn) %>%
  summarise(foreign.med=median(p.nonlocal),foreign.lcl=quantile(p.nonlocal,0.025), foreign.ucl=quantile(p.nonlocal,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.ctn
#fwrite(out.ctn,"output/woco_nonlocal_origin_estimates_CANTON_comb_prior.csv")




## 4.2. summarise output in graphical form ------------------------------------------



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



mean.p.nonlocal %>%
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
        legend.position.inside=c(0.90,0.35),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))





# 5. SPECIFY SIMPLE PROBABILITY MODEL TO ESTIMATE PROBABILITY OF LOCAL ORIGIN IN NIMBLE ----------------------------------------------

# uncertainty explodes when you move from a single global prior on p to individual priors p[i] inside the loop is due to how hierarchical structure and parameterization affect identifiability and shrinkage.
# 
# When you use one global p:
#   
#   All birds share the same prior probability of being local.
# The model pools information across individuals → strong shrinkage → posterior uncertainty is reduced.
# 
# 
# When you switch to p[i] for each bird:
#   
#   Each bird now has its own Bernoulli probability parameter.
# You’ve introduced 800 extra parameters with only one observation each.
# There is almost no data to inform π[i] beyond its prior, so the posterior for each π[i] remains close to the prior (Beta(1,1) = uniform).
# Result: extreme uncertainty (posterior ranges 0–1) because the model cannot learn much about π[i] from a single observation.


woco.orig.simple<-nimbleCode({
  
  # Parameters:
  # p.local: probability of local origin (that shot woodcock was a local bird) - varies by individual
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  p.local ~ dunif(0, 1)  ## completely uninformative prior
  for (i in 1:nind.unkn){
    

    z[i] ~ dbern(p.local)  ## indicator variable from binomial draw of which isotope distribution fits better
    
    # # Potential local distribution based on isotope ratio
    # mu.unknown[i,1] <- mean.feather.d2H[age[i]+1]    ### overall distribution across Europe
    # mu.unknown[i,] <- d2H_feather.local[i]  ### local expected feather isotope distribution
    # 
    # sd.unknown[i,2] <- sd.feather.d2H[age[i]+1]    ### sd of overall distribution across Europe
    # sd.unknown[i,1] <- d2H_feather.local.sd[i]  ### sd of local expected feather isotope distribution
    # 
    
    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mean = d2H_feather.local[i]*z[i] + mean.feather.d2H[age[i]+1]*(1 - z[i]),
                                    sd = d2H_feather.local.sd[i]*z[i] + sd.feather.d2H[age[i]+1]*(1 - z[i]))
    
    
    
  } #i
  
  
}) ## end of nimble code chunk


iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      mean.feather.d2H = c(mean.rain.d2H.ad,mean.rain.d2H.juv),
                      sd.feather.d2H = c(sd.rain.d2H.ad,sd.rain.d2H.juv),
                      d2H_feather.local=woco.unk.sf$d2h_local_predicted,
                      d2H_feather.local.sd=woco.unk.sf$d2h_local_sd,
                      age = ifelse(woco.unk.sf$AGE=="Adulte",0,1))   ## migration weeks start only in August

iso.data <- list(d2H_feather.unknown = woco.unk.sf$dH)
parameters.iso <- c("p.local") #,"p.nonlocal.prior","p.nonlocal.prior1") including these bloats the output object
iso.inits <- list(z = ifelse(woco.unk.sf$dH < woco.unk.sf$d2h_local_predicted,0,1),
                  p.local = 1 # initial start value will be replaced after test run
)

## explore whether precise assignment was possible if variation wasn't so large
iso.constants$d2H_feather.local<-iso.constants$d2H_feather.local+50
iso.constants$sd.feather.d2H<-iso.constants$sd.feather.d2H*0.1
iso.constants$d2H_feather.local.sd<-iso.constants$d2H_feather.local.sd*0.1
iso.data$d2H_feather.unknown<-iso.constants$mean.feather.d2H[iso.constants$age+1]


woco.iso.simple <- nimbleMCMC(code = woco.orig.simple,
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

# compile all the samples
samples <- rbind(woco.iso.simple$samples$chain1,woco.iso.simple$samples$chain2,woco.iso.simple$samples$chain3,woco.iso.simple$samples$chain4)
mean.p.local <- as_tibble(samples[,grep("p.local\\[", colnames(samples))]) %>%
  gather(key="parameter",value="p.local") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind])


# summarise by Canton

out.ctn<- mean.p.local %>%
  group_by(age,ctn) %>%
  summarise(foreign.med=median(p.local),foreign.lcl=quantile(p.local,0.025), foreign.ucl=quantile(p.local,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.ctn






# 6. SPECIFY DETERMINISTIC ASSIGNMENT MODEL WITH BAYES FACTOR  ----------------------------------------------


woco.orig.bf<-nimbleCode({
  
  
  for (i in 1:nind.unkn){
  # Parameters:
  # p.nonlocal: probability of nonlocal origin (that shot woodcock was a local bird) - varies by individual
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
        # Potential local distribution based on timing used as prior (from WOCO_migration_model.R - could be integrated into this model)
        logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
          b.mig.week*(unk.week[i])     ### week effect from migration model
        p.local.prior1[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from migration model
        p.local.prior[i] <- max((1-p.local.prior1[i]),p.local.prior2[i])  ## combining priors by taking the maximum (once birds have migrated the highest prob will do)

    # Compute densities
    local_density[i] <- dnorm(d2H_feather.unknown[i], mean = d2H_feather.local[i], sd = d2H_feather.local.sd[i])
    global_density[i] <- dnorm(d2H_feather.unknown[i], mean = mean.feather.d2H[age[i]+1], sd = sd.feather.d2H[age[i]+1])
    
    # Posterior probability using Bayes formula
    p.local[i] <- (p.local.prior[i] * local_density[i]) / 
      (p.local.prior[i] * local_density[i] + (1-p.local.prior[i]) * global_density[i])
    
  } #i
  
  
}) ## end of nimble code chunk


    iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                          mean.feather.d2H = c(mean.rain.d2H.ad,mean.rain.d2H.juv),
                          sd.feather.d2H = c(sd.rain.d2H.ad,sd.rain.d2H.juv),
                          d2H_feather.local=woco.unk.sf$d2h_local_predicted,
                          d2H_feather.local.sd=woco.unk.sf$d2h_local_sd,
                          age = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                          p.local.prior2 = (1-woco.unk.sf$abd_prior),
                          lm.mean.mig=logit(as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[1,1])),
                          b.mig.week=as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[2,1]),
                          unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26"))   ## migration weeks start only in August
    )
    
    iso.data <- list(d2H_feather.unknown = woco.unk.sf$dH)

    parameters.iso <- c("p.local") #,"p.nonlocal.prior","p.nonlocal.prior1") including these bloats the output object
    
    # Initial values  FOR ALL PARAMETERS
    ## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8
    
    iso.inits <- list(p.local = (1-woco.unk.sf$abd_prior)) # initial start value will be replaced after test run

    ## explore whether precise assignment was possible if variation wasn't so large
    #iso.constants$d2H_feather.local<-iso.constants$d2H_feather.local+50
    iso.constants$sd.feather.d2H<-iso.constants$sd.feather.d2H*0.1
    iso.constants$d2H_feather.local.sd<-iso.constants$d2H_feather.local.sd*0.1
    #iso.data$d2H_feather.unknown<-iso.constants$mean.feather.d2H[iso.constants$age+1]
    
    
    woco.iso <- nimbleMCMC(code = woco.orig.bf,
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
    
    
    
    
    
    # compile all the samples
    samples <- rbind(woco.iso$samples$chain1,woco.iso$samples$chain2,woco.iso$samples$chain3,woco.iso$samples$chain4)
    mean.p.local <- as_tibble(samples[,grep("p.local\\[", colnames(samples))]) %>%
      gather(key="parameter",value="p.local") %>%
      mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
      mutate(age=iso.constants$age[ind]) %>%
      mutate(ctn=woco.unk.sf$KANTON[ind]) 
    
    out.ctn<- mean.p.local %>%
      group_by(age,ctn) %>%
      summarise(foreign.med=median(p.local),foreign.lcl=quantile(p.local,0.025), foreign.ucl=quantile(p.local,0.975)) %>%
      mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
      ungroup() %>%
      select(-age) %>%
      mutate(prior="combined abundance and migration")
    out.ctn

ggplot(data=out.ctn,aes(x=ctn, y=foreign.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=foreign.lcl, ymax=foreign.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  
  ## format axis ticks
  labs(y="Proportion of shot woodcocks of local origin",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
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
        legend.position.inside=c(0.90,0.35))




# 7. CALCULATE LOCAL PROBABILITY IN R FROM LOCAL DISTRIBUTION ONLY ------------------------

local_density <- dnorm(woco.unk.sf$dH, mean = woco.unk.sf$d2h_local_predicted, sd = woco.unk.sf$d2h_local_sd)
global_density.ad <- dnorm(woco.unk.sf$dH[woco.unk.sf$AGE=="Adulte"], mean = mean.rain.d2H.ad, sd = sd.rain.d2H.ad)
global_density.juv <- dnorm(woco.unk.sf$dH[woco.unk.sf$AGE=="Adulte"], mean = mean.rain.d2H.juv, sd = sd.rain.d2H.juv)
global.density<-c(global_density.ad,global_density.juv)

# Posterior probability using Bayes formula
p.local.naive <- local_density / (local_density + global_density)
p.local <- ((1-woco.unk.sf$abd_prior) * local_density) / 
  ((1-woco.unk.sf$abd_prior) * local_density + woco.unk.sf$abd_prior * global_density)

hist(p.local.naive)

out.ctn<-as_tibble(p.local.naive) %>%
  rename(p.local=value) %>%
  mutate(age=iso.constants$age) %>%
  mutate(ctn=woco.unk.sf$KANTON) %>%
  group_by(age,ctn) %>%
  summarise(local.med=median(p.local),local.lcl=quantile(p.local,0.025), local.ucl=quantile(p.local,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.ctn

ggplot(data=out.ctn,aes(x=ctn, y=local.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=local.lcl, ymax=local.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  
  ## format axis ticks
  labs(y="Proportion of shot woodcocks of local origin",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
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
        legend.position.inside=c(0.90,0.85))




# 8. SIMULATE LOCAL PROBABILITY IN R WHEN UNCERTAINTY IS MUCH REDUCED ------------------------

local_density <- dnorm(woco.unk.sf$dH, mean = woco.unk.sf$d2h_local_predicted, sd = woco.unk.sf$d2h_local_sd)
global_density.ad <- dnorm(woco.unk.sf$dH[woco.unk.sf$AGE=="Adulte"], mean = mean.rain.d2H.ad, sd = sd.rain.d2H.ad*0.1)
global_density.juv <- dnorm(woco.unk.sf$dH[woco.unk.sf$AGE=="Adulte"], mean = mean.rain.d2H.juv, sd = sd.rain.d2H.juv*0.1)
global.density<-c(global_density.ad,global_density.juv)

# Posterior probability using Bayes formula
p.local.mig.prior<-plogis(logit(as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[1,1])) +
                          as.numeric(fread("output/woco_telemetry_seasonal_surv_parm.csv")[2,1])*(week(UNK_WC$DATE)-week(ymd("2024-07-26"))))    ### week effect from migration model
p.local.naive <- local_density / (local_density + global_density)
p.local.abd <- ((1-woco.unk.sf$abd_prior) * local_density) / 
  ((1-woco.unk.sf$abd_prior) * local_density + woco.unk.sf$abd_prior * global_density)
p.local.mig <- ((1-p.local.mig.prior) * local_density) / 
  ((1-p.local.mig.prior) * local_density + p.local.mig.prior * global_density)

hist(p.local.naive)
hist(p.local.abd)
hist(p.local.mig)


out.ctn<-as_tibble(p.local.abd) %>%
  rename(p.local=value) %>%
  mutate(age=iso.constants$age) %>%
  mutate(ctn=woco.unk.sf$KANTON) %>%
  group_by(age,ctn) %>%
  summarise(n=length(p.local),local.med=median(p.local),local.lcl=quantile(p.local,0.025), local.ucl=quantile(p.local,0.975)) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  ungroup() %>%
  select(-age) %>%
  mutate(prior="combined abundance and migration")
out.ctn

ggplot(data=out.ctn,aes(x=ctn, y=local.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=local.lcl, ymax=local.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  
  ## format axis ticks
  labs(y="Proportion of shot woodcocks of local origin",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
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
        legend.position.inside=c(0.90,0.85))
