###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS --------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds
## background: https://www.ffw.ch/de/news/waldschnepfe/

## for reporting: estimate proportion that has departed at time when hunting season opens

# Clear workspace ---------------------------------------------------------

rm(list=ls()) 
Sys.setenv(LANG = "en") ##change language of error messages to english
# install.packages(c('tictoc','coda','MCMCvis','basicMCMCplots','lubridate','tidyverse','dplyr','data.table','png','magick','jpeg','grid','gtable','igraph', 'coda', 'R6', 'pracma', 'numDeriv'), dependencies=T)
# install.packages('nimble', dependencies=T)
# Load libraries, functions etc -------------------------------------------
library(lubridate)
library(tidyverse)
library(dplyr)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select
library(coda)
library(MCMCvis)
library(nimble)
library(tictoc)
library(basicMCMCplots) # for trace plots called chainsPlot

## set root folder for project
#setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")
setwd("C:/woodcock")  ## for HPC


# ## set PATH environment for nimble
# path <- Sys.getenv('PATH')
# newPath <- paste("C:\\rtools44\\mingw32\\bin;C:\\rtools44\\mingw64\\bin;",
#              path, sep = "")
# Sys.setenv(PATH = newPath) 


# Load data ---------------------------------------------------------------
# prepared in script WOCO_survival_data_compilation.R
load("data/woco_mig_input.RData")
nimbleOptions(allowDynamicIndexing = TRUE)



## ADJUST INPUT DATA (provided by Jaume Badia Boher)

# remove NAs in the encounter histories to remove the number of NAs when assembling the model with nimbleModel, and will make error tracing easier.

for(i in 1:nrow(y.telemetry)){
  
  if(f.telemetry[i] > 1){
    
    y.telemetry[i, 1:(f.telemetry[i]-1)] <- 5
    
  }
  
}

# add an absorbing state to the model (long dead)
# probably more correct because otherwise a dead observation is always under the estimation of a recovery probability
# Thus, recode the observation and state matrix

ww3 <- which(y.telemetry == 3, arr.ind = T)

for(i in 1:nrow(ww3)){
  
  if(ww3[i,2] != nweeks){
    
    y.telemetry[ww3[i,1], (ww3[i,2] + 1):nweeks] <- 5 ## state 5 for not observed
    z.telemetry[ww3[i,1], (ww3[i,2] + 1):nweeks] <- 5 ## state 5 for long dead
    
  }
  
}

ww4 <- which(y.telemetry == 4, arr.ind = T)

for(i in 1:nrow(ww4)){
  
  if(ww4[i,2] != nweeks){
    
    y.telemetry[ww4[i,1], (ww4[i,2] + 1):nweeks] <- 5 ## state 5 for not observed
    z.telemetry[ww4[i,1], (ww4[i,2] + 1):nweeks] <- 5 ## state 5 for long dead
    
  }
  
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP NIMBLE CODE FOR SURVIVAL MODEL ---------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify model in NIMBLE format
woco.mig.model<-nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi: weekly survival probability intercept
  # mig: weekly migration probability
  
  # p.obs: probability to be observed alive even through capture or transmitter
  # p.found.dead: probability for carcass to be recovered
  
  # -------------------------------------------------
  # True States (S) - these are often unknown and cannot be observed
  # 1 alive in study area
  # 2 dead in study area
  # 3 alive outside study area (after migration)
  # 4 dead outside study area (shot after migration)
  # 5 absorbing state / long-dead
  
  # Observed States (O) - these are based on the actual transmission history
  # 1 Bird (re)captured alive or recorded via transmitter in study area
  # 2 Bird recorded via transmitter outside study area
  # 3 Dead bird recovered  in study area
  # 4 Dead bird recovered  outside study area
  # 5 No signal (=not observed)
  
  # -------------------------------------------------
  
  # Priors and constraints
  
  
  #### WEEKLY SURVIVAL, MIGRATION and DETECTION PROBABILITIES
  for (i in 1:nind){
    for (t in f[i]:(nweeks)){
      phi[i,t] <- mean.phi[year[i]] #+      ###
      p.dead.in[i,t] <- mean.p.dead.in[year[i]] #+      ###
      p.dead.out[i,t] <- mean.p.dead.out[year[i]] #+      ###
      
      logit.mig[i,t] <- lm.mean +      ### intercept for mean survival
        b.mig.week*(week[t])     ### migration probability dependent on week
      mig[i,t] <-ilogit(logit.mig[i,t])
      
      logit.p.obs.in[i,t] <- lpin.mean +      ### intercept for mean survival
        b.obs.effort*(effort[i,t]) +    ### observation dependent on effort in that week and year
        b.obs.tag*(tag[i])    ### observation dependent on whether animal was tagged
      p.obs.in[i,t] <-ilogit(logit.p.obs.in[i,t])
      
      logit.p.obs.out[i,t] <- lpout.mean +      ### intercept for mean survival
        b.obs.tag*(tag[i])    ### observation dependent on whether animal was tagged
      p.obs.out[i,t] <-ilogit(logit.p.obs.out[i,t])
      
    } #t
  } #i
  
  #### BASELINE FOR SURVIVAL PROBABILITY (varies by year)
  for(s in 1:nyears) {
    mean.phi[s] ~ dunif(0.9,1)   # fairly uninformative prior for weekly survival probabilities
    #lp.mean[s] <- log(mean.phi[s]/(1 - mean.phi[s]))    # logit transformed survival intercept
    
    #### DEAD RECOVERY PROBABILITY VARIES BY YEAR
    mean.p.dead.in[s] ~ dunif(0,1)
    mean.p.dead.out[s] ~ dunif(0,0.3)
  }
  
  #### BASELINE FOR MIGRATION PROBABILITY (varies by year)
  mean.mig ~ dnorm(0,1)   # fairly uninformative prior for intercept of weekly departure probability
  lm.mean <- log(mean.mig/(1 - mean.mig))    # logit transformed migration intercept
  
  #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
  mean.p.in ~ dunif(0,1)   # fairly uninformative prior for weekly detection probabilities
  lpin.mean <- log(mean.p.in/(1 - mean.p.in))    # logit transformed detection intercept
  
  mean.p.out ~ dunif(0,0.3)   # fairly uninformative prior for weekly detection probabilities
  lpout.mean <- log(mean.p.out/(1 - mean.p.out))    # logit transformed detection intercept
  
  #### SLOPE PARAMETERS FOR PROBABILITIES ON LOGIT SCALE
  b.mig.week ~ dnorm(1,sd=2)         # Prior for week effect on migration probability on logit scale - must be positive
  b.obs.effort ~ dnorm(1, sd=2)         # Prior for effort effect on observation probability on logit scale  - must be positive
  b.obs.tag ~ dnorm(1, sd=2)         # Prior for tag effect on observation probability on logit scale  - must be positive 
  
  # #### DEAD RECOVERY PROBABILITY VARIES BY INDIVIDUAL
  # for (i in 1:nind){
  #   p.found.dead.in[i] ~ dunif(0,1)
  #   p.found.dead.out[i] ~ dunif(0,0.3)
  # } #i
  
  
  
  # -------------------------------------------------
  # Define state-transition and observation matrices 
  # -------------------------------------------------
  
  for (i in 1:nind){
    
    for (t in f[i]:(nweeks-1)){
      
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      
      ps[1,i,t,1]<-phi[i,t]*(1-mig[i,t])
      ps[1,i,t,2]<-(1-phi[i,t])*(1-mig[i,t])
      ps[1,i,t,3]<-phi[i,t]*mig[i,t]
      ps[1,i,t,4]<-(1-phi[i,t])*mig[i,t]
      ps[1,i,t,5] <- 0
      
      ps[2,i,t,1]<-0
      ps[2,i,t,2]<-0
      ps[2,i,t,3]<-0
      ps[2,i,t,4]<-0
      ps[2,i,t,5]<-1
      
      ps[3,i,t,1]<-0
      ps[3,i,t,2]<-0
      ps[3,i,t,3]<-phi[i,t]
      ps[3,i,t,4]<-(1-phi[i,t])
      ps[3,i,t,5] <- 0
      
      ps[4,i,t,1]<-0
      ps[4,i,t,2]<-0
      ps[4,i,t,3]<-0
      ps[4,i,t,4]<-0
      ps[4,i,t,5] <- 1
      
      #Add an absorbing state, otherwise the recovery probability is applied all over and over I think
      
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0
      ps[5,i,t,3] <- 0
      ps[5,i,t,4] <- 0
      ps[5,i,t,5] <- 1 #Probability (dead to remain dead)
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      
      po[1,i,t,1]<-p.obs.in[i,t]
      po[1,i,t,2]<-0
      po[1,i,t,3]<-0
      po[1,i,t,4]<-0
      po[1,i,t,5]<-(1-p.obs.in[i,t])
      
      po[2,i,t,1]<-0
      po[2,i,t,2]<-0
      po[2,i,t,3]<-p.dead.in[i,t]
      po[2,i,t,4]<-0
      po[2,i,t,5]<-(1-p.dead.in[i,t])
      
      po[3,i,t,1]<-0
      po[3,i,t,2]<-p.obs.out[i,t]
      po[3,i,t,3]<-0
      po[3,i,t,4]<-0
      po[3,i,t,5]<-(1-p.obs.out[i,t])
      
      po[4,i,t,1]<-0
      po[4,i,t,2]<-0
      po[4,i,t,3]<-0
      po[4,i,t,4]<-p.dead.out[i,t]
      po[4,i,t,5]<-(1-p.dead.out[i,t])
      
      #Add match event - absorbing state: long-dead individuals will always be unobservable with prob. 1
      
      po[5,i,t,1] <- 0
      po[5,i,t,2] <- 0
      po[5,i,t,3] <- 0
      po[5,i,t,4] <- 0
      po[5,i,t,5] <- 1
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1 ## alive in study area on first record - all other individuals and observations should have been eliminated because they are not relevant to the model
    for (t in (f[i]+1):nweeks){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:5])     ### check different CJS distribution in nimbleEcology dcjs
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])
      #rep.states[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])    ## include this as goodness of fit measure
    } #t
  } #i
}) ## end of nimble code chunk




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE NIMBLE RUN WITH DATA INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### BUNDLE DATA INTO A LIST

#### DISTINGUISH CONSTANTS AND DATA----
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
telemetry.constants <- list(f = f.telemetry,
                            effort = effort,  # must be same dimension as obs matrix
                            week = seq(1:nweeks),
                            year = year,
                            tag = tag,
                            nind = nind,
                            nweeks=nweeks,
                            nyears = nyears)
telemetry.data <- list(y = y.telemetry)


### SPECIFY DIMENSIONS for arrays used with empty indices
telemetry.dims <- list(ps = c(5, dim(y.telemetry)[1], dim(y.telemetry)[2], 5),
                       po = c(5, dim(y.telemetry)[1], dim(y.telemetry)[2], 5))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters monitored
parameters.telemetry <- c("mean.mig","b.mig.week","mean.phi")

# Initial values for some parameters
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

# Missing values (NAs) or non-finite values were found in model variables: phi, b.phi.age, b.phi.mig, b.phi.feed, b.phi.sex, b.phi.FRA, b.phi.SUI, b.phi.ESP, p.obs, tag.fail, p.found.dead, ps, po, z, y, rep.states.
smartInit1 <- list(z = z.telemetry,
                  # z = ifelse(is.na(z.telemetry),1,z.telemetry),
                  # y = ifelse(is.na(y.telemetry),1,y.telemetry),
                  #  y = y.telemetry,
                   
                   mean.phi = runif(nyears,0.95,0.999),
                   mean.p.dead.in = runif(nyears,0.2,0.9),
                   mean.p.dead.out = runif(nyears,0,0.3),
                   
                   # phi = array(runif(nweeks*nind,0.95,0.999),dim=c(nind,nweeks)),
                   # p.dead.in = array(runif(nweeks*nind,0.2,0.9),dim=c(nind,nweeks)),
                   # p.dead.out = array(runif(nweeks*nind,0,0.3),dim=c(nind,nweeks)),
                   
                   #### BASELINE FOR MIGRATION PROBABILITY (varies by year)
                   mean.mig = 0.1,   # fairly uninformative prior for weekly migration probabilities

                   #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
                   mean.p.in = 0.2,   # fairly uninformative prior for weekly detection probabilities
                   mean.p.out = 0.05,   # fairly uninformative prior for weekly detection probabilities
                   
                   #### SLOPE PARAMETERS FOR PROBABILITIES ON LOGIT SCALE
                   b.mig.week = 0.5,         # Prior for week effect on migration probability on logit scale - must be positive
                   b.obs.effort = 0.5,         # Prior for effort effect on observation probability on logit scale  - must be positive
                   b.obs.tag = 0.5         # Prior for tag effect on observation probability on logit scale  - must be positive 
                   
)



# # MCMC settings
# # number of posterior samples per chain is n.iter - n.burnin
# n.iter <- 200
# n.burnin <- 100
# n.chains <- 4
# 
# 
# # PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
# test <- nimbleModel(code = woco.mig.model,
#                     constants=telemetry.constants,
#                     data = telemetry.data,
#                     inits = smartInit1,
#                     calculate=TRUE)
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
# test$initializeInfo()
# #help(modelInitialization)
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits
# 
# 
# # use test output as starting values or check where the NA comes from
# test$logProb_b.mig.week
# test$logProb_mean.phi
# test$b.obs.effort
# test$mean.phi
# test$phi[1:3,]
# y.telemetry[1:3,15:16]
# z.telemetry[1:3,15:16]
# test$logProb_y ### there should not be any -Inf in this matrix
# y.telemetry[17,15:18]
# z.telemetry[17,]



### FIT FULL MODEL TO ALL DATA ############################
# RUN NIMBLE MODEL --------------------

# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 50000
n.burnin <- 25000
n.chains <- 4

tic()
### this takes 4000 sec (2 hrs) for 20000 iterations and converges in that time
woco_surv <- nimbleMCMC(code = woco.mig.model,
                        constants=telemetry.constants,
                        data = telemetry.data,
                        inits = smartInit1,
                        dimensions=telemetry.dims,
                        monitors = parameters.telemetry,
                        thin=4,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains,
                        progressBar = getNimbleOption("MCMCprogressBar"),
                        summary=T)
toc() 


saveRDS(woco_surv,"output/woco_mig_depart_output_nimble.rds")
#woco_surv<-readRDS("output/woco_seasonal_survival_output_nimble.rds")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMINE OUTPUT AND DIAGNOSTICS WITH MCMCvis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### MODEL ASSESSMENT ####
out<- as.data.frame(MCMCsummary(woco_surv$samples, params=c("mean.mig","b.mig.week","mean.phi")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_telemetry_seasonal_surv_parm.csv")
#out<-fread("output/woco_telemetry_surv_parm_v1.csv")


MCMCplot(woco_surv$samples, params=c("mean.mig","b.mig.week","mean.phi"))
ggsave("output/woco_seasonal_survival_parameter_estimates.jpg", height=9, width=9)
# ggsave("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/wocopopmod/output/Fig_S1_parameter_estimates.jpg", height=11, width=8)

## look at the chains and whether they mixed well
chainsPlot(woco_surv$samples,
           var=c("mean.mig","b.mig.week"
           ))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ MADE UP GOF TEST  #############################------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### suggested by Pius Korner in July 2023
### simulate data from model
### compare frequency of states from simulated prediction to frequency in observed data
### run once, then commented out as rep.states slowed down convergence

# OBS<-table(as.factor(y.telemetry))
# REPraw<-MCMCpstr(woco_surv$samples, params=c("rep.states"), type="chains")
# REP<-table(as.factor(as.numeric(REPraw$rep.states[REPraw$rep.states>0])))  ## remove the 0s because they are NA in the y.telemetry matrix from occ before an animal was tagged
# chisq.test(OBS,REP)  ### should have a p-value >> 0.05 otherwise there would be disconcerting lack of fit




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ ESTIMATE DEPARTURE TIME OF WOODCOCKS  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### PREPARE RAW MCMC OUTPUT
parmcols<-dimnames(woco_surv$samples[[1]])[[2]]

# ### COMBINE SAMPLES ACROSS CHAINS
MCMCout<-rbind(woco_surv$samples[[1]],woco_surv$samples[[2]],woco_surv$samples[[3]],woco_surv$samples[[4]])
# str(MCMCout)

### SET UP ANNUAL TABLE FOR PLOTTING THE ANNUAL SURVIVAL GRAPH
AnnTab<-data.frame(week = seq(1:nweeks))

### CALCULATE PREDICTED VALUE FOR EACH SAMPLE
MCMCpred<-data.frame()
for(s in 1:nrow(MCMCout)) {
  
  X<-  AnnTab %>%
    
    ##CALCULATE MONTHLY SURVIVAL
    mutate(logit.mig=logit(as.numeric(MCMCout[s,match("mean.mig",parmcols)]))+
             as.numeric(MCMCout[s,match("b.mig.week",parmcols)])*week)  %>%
    
    ##BACKTRANSFORM TO NORMAL SCALE
    mutate(mig=plogis(logit.mig)) %>%
    mutate(simul=s)              
  
  MCMCpred<-rbind(MCMCpred,as.data.frame(X)) 
}

### CREATE PLOT

FIGURE<-MCMCpred %>% rename(raw.mig=mig) %>%
  group_by(week) %>%
  summarise(mig=quantile(raw.mig,0.5),mig.lcl=quantile(raw.mig,0.025),mig.ucl=quantile(raw.mig,0.975)) %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>%
  
  ggplot()+
  geom_ribbon(aes(x=Date, ymin=mig.lcl, ymax=mig.ucl), alpha=0.2, fill="firebrick") +   ##
  geom_line(aes(x=Date, y=mig),linewidth=1, col="firebrick")+     ##
  
  ## format axis ticks
  scale_x_date(name="Week of the year", date_labels = "%b %d") +
  scale_y_continuous(name="Probability to leave study area", limits=c(0,1), breaks=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE

ggsave(plot=FIGURE,
       filename="output/woco_seasonal_mig.jpg", 
       device="jpg",width=11, height=8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ CUMULATIVE FIGURE: PROP POPULATION THAT HAS DEPARTED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## first calculate what proportion of a population of 1000 birds has departed by a given date
woco.mig<-MCMCpred %>% rename(raw.mig=mig) %>%
  mutate(pop=1000, left=0, prop_mig=0)

for(s in 1: max(woco.mig$simul)){
  for(w in 1: nweeks){
    woco.mig$left[woco.mig$week==w & woco.mig$simul==s]<-rbinom(1,
                                                                woco.mig$pop[woco.mig$week==w & woco.mig$simul==s],
                                                                woco.mig$raw.mig[woco.mig$week==w & woco.mig$simul==s])
    woco.mig$prop_mig[woco.mig$week==w & woco.mig$simul==s]<-sum(woco.mig$left[woco.mig$week<=w & woco.mig$simul==s])/1000
    woco.mig$pop[woco.mig$week==w+1 & woco.mig$simul==s]<-woco.mig$pop[woco.mig$week==w & woco.mig$simul==s]-
                                                                woco.mig$left[woco.mig$week==w & woco.mig$simul==s]
  }
}

saveRDS(woco.mig,"output/woco_mig_depart_simulation.rds")
#woco_mig<-readRDS("output/woco_mig_depart_simulation.rds")



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

FIGURE2<-woco.mig %>% 
  group_by(week) %>%
  summarise(mig=quantile(prop_mig,0.5),mig.lcl=quantile(prop_mig,0.025),mig.ucl=quantile(prop_mig,0.975)) %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>%
  
  ggplot()+
  geom_ribbon(aes(x=Date, ymin=mig.lcl, ymax=mig.ucl), alpha=0.2, fill="firebrick") +   ##
  geom_line(aes(x=Date, y=mig),linewidth=1, col="firebrick")+     ##
  
  ### add vertical lines to specify key dates
  geom_vline(aes(xintercept=min(Date[mig>0.95])), linetype="dashed", col="forestgreen", linewidth=1.5) +
  geom_segment(x=lubridate::ymd("2024-09-15"),y=0,yend=0.6, linetype="dashed", col="grey36", linewidth=2) +
  geom_segment(x=lubridate::ymd("2024-10-01"),y=0,yend=0.6, linetype="dashed", col="grey27", linewidth=2) +
  geom_segment(x=lubridate::ymd("2024-10-15"),y=0,yend=0.6, linetype="dashed", col="grey18", linewidth=2) +
  geom_text(x=lubridate::ymd("2024-09-15"),y=0.65,label = "JU\nNE", size=6,col="grey36", vjust = 'bottom')+
  geom_text(x=lubridate::ymd("2024-10-01"),y=0.65,label = "BE\nVS", size=6,col="grey27", vjust = 'bottom')+
  geom_text(x=lubridate::ymd("2024-10-15"),y=0.65,label = "FR\nTI\nVD", size=6,col="grey18", vjust = 'bottom')+
  
  
  ### add the bird icons
  annotation_custom(gunicon, xmin=lubridate::ymd("2024-08-15"), xmax=lubridate::ymd("2024-09-07"), ymin=0.6, ymax=0.8)+
  annotation_custom(wocoicon, xmin=lubridate::ymd("2024-08-01"), xmax=lubridate::ymd("2024-09-01"), ymin=0.8, ymax=1)+
  
  
  ## format axis ticks
  scale_x_date(name="Week of the year", date_labels = "%d %b") +
  scale_y_continuous(name="Cumulative proportion of woodcocks that have left study area", limits=c(0,1), breaks=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE2

ggsave(plot=FIGURE2,
       filename="output/woco_cumulative_mig_prop.jpg", 
       device="jpg",width=11, height=8)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ CALCULATE PROP POPULATION THAT HAS DEPARTED BY CERAIN DATES -------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Nationale Ebene 15 September: https://www.fedlex.admin.ch/eli/cc/1988/506_506_506/de
# Kanton BE und VS 1 Okt
# Kanton FR, TI, and VD 15 Okt
# Kanton JU und NE 15 Sept

# most woodcock are shot between 15 Oct and 15 Nov, so include 15 Nov as well

OUTPUT<-woco.mig %>%
  group_by(week) %>%
  summarise(mig=quantile(prop_mig,0.5),mig.lcl=quantile(prop_mig,0.025),mig.ucl=quantile(prop_mig,0.975)) %>%
  mutate(rem=1-mig,rem.lcl=1-mig.lcl,rem.ucl=1-mig.ucl) %>%
  mutate(WeekStartDate=lubridate::ymd("2024-07-28") + lubridate::weeks(week - 1)) %>%
  filter(week %in% c(8,10,12,17)) %>%
  mutate(Percent_departed=paste(round(mig*100,1)," (",round(mig.lcl*100,1)," - ",round(mig.ucl*100,1),")", sep="")) %>%
  mutate(Percent_remaining=paste(round(rem*100,1)," (",round(rem.lcl*100,1)," - ",round(rem.ucl*100,1),")", sep="")) %>%
  select(WeekStartDate,Percent_departed,Percent_remaining)

fwrite(OUTPUT,"output/WOCO_prop_departed_key_hunting_dates.csv")
