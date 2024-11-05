###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS --------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds
## background: https://www.ffw.ch/de/news/waldschnepfe/

#Checked by Jaume A. Badia-Boher, 05.11.24

#Fixed errors (see code as well):

#First encounters in the observed data took values different than 1 - I changed them to 1s
#Model formulation retained dead individuals in a state where recovery probabilities were applied at all years. This would probably lead to a bias in the estimation of recovery probabilities. I added an absorbing state (the code now has 5 states)
#There were two individuals showing forbidden transitions (detected in Switzerland, then abroad, then again in Switzerland)
#There were errors in the initial values for the latent state matrix z - I rewrote it

#There is still an annoying error: the likelihood for the observations (matrix y) is still providing "Inf" values at some locations. Interestingly, these are just for some of the individuals first marked at t=1 (not all of them), and the loglikelihood issue is always at t = 2. I actually do not get what is going on, because there are some individuals first marked at t=1 and with similar capture histories for which the loglikelihood can be successfully computed. There is an "effort" covariate for detection probability, can this be generating this issue? See list of individuals with loglikelihood issues here:


#warning: logProb of data node y[32, 2]: logProb is -Inf.
#warning: logProb of data node y[66, 2]: logProb is -Inf.
#warning: logProb of data node y[67, 2]: logProb is -Inf.
#warning: logProb of data node y[73, 2]: logProb is -Inf.
#warning: logProb of data node y[74, 2]: logProb is -Inf.
#warning: logProb of data node y[77, 2]: logProb is -Inf.
#warning: logProb of data node y[78, 2]: logProb is -Inf.
#warning: logProb of data node y[80, 2]: logProb is -Inf.
#warning: logProb of data node y[82, 2]: logProb is -Inf.
#warning: logProb of data node y[87, 2]: logProb is -Inf.
#warning: logProb of data node y[91, 2]: logProb is -Inf.
#warning: logProb of data node y[92, 2]: logProb is -Inf.
#warning: logProb of data node y[95, 2]: logProb is -Inf.
#warning: logProb of data node y[96, 2]: logProb is -Inf.
#warning: logProb of data node y[105, 2]: logProb is -Inf.
#warning: logProb of data node y[106, 2]: logProb is -Inf.
#warning: logProb of data node y[109, 2]: logProb is -Inf.
#warning: logProb of data node y[113, 2]: logProb is -Inf.
#warning: logProb of data node y[114, 2]: logProb is -Inf.
#warning: logProb of data node y[119, 2]: logProb is -Inf.
#warning: logProb of data node y[120, 2]: logProb is -Inf.
#warning: logProb of data node y[122, 2]: logProb is -Inf.
#warning: logProb of data node y[123, 2]: logProb is -Inf.
#warning: logProb of data node y[153, 2]: logProb is -Inf.
#warning: logProb of data node y[154, 2]: logProb is -Inf.
#warning: logProb of data node y[155, 2]: logProb is -Inf.
#warning: logProb of data node y[156, 2]: logProb is -Inf.
#warning: logProb of data node y[157, 2]: logProb is -Inf.
#warning: logProb of data node y[158, 2]: logProb is -Inf.
#warning: logProb of data node y[159, 2]: logProb is -Inf.
#warning: logProb of data node y[160, 2]: logProb is -Inf.
#warning: logProb of data node y[161, 2]: logProb is -Inf.
#warning: logProb of data node y[162, 2]: logProb is -Inf.
#warning: logProb of data node y[163, 2]: logProb is -Inf.
#warning: logProb of data node y[164, 2]: logProb is -Inf.
#warning: logProb of data node y[165, 2]: logProb is -Inf.
#warning: logProb of data node y[166, 2]: logProb is -Inf.
#warning: logProb of data node y[167, 2]: logProb is -Inf.
#warning: logProb of data node y[168, 2]: logProb is -Inf.
#warning: logProb of data node y[169, 2]: logProb is -Inf.
#warning: logProb of data node y[170, 2]: logProb is -Inf.
#warning: logProb of data node y[171, 2]: logProb is -Inf.
#warning: logProb of data node y[172, 2]: logProb is -Inf.
#warning: logProb of data node y[173, 2]: logProb is -Inf.
#warning: logProb of data node y[174, 2]: logProb is -Inf.
#warning: logProb of data node y[176, 2]: logProb is -Inf.
#warning: logProb of data node y[179, 2]: logProb is -Inf.
#warning: logProb of data node y[184, 2]: logProb is -Inf.
#warning: logProb of data node y[215, 2]: logProb is -Inf.
#warning: logProb of data node y[248, 2]: logProb is -Inf.
#warning: logProb of data node y[260, 2]: logProb is -Inf.
#warning: logProb of data node y[261, 2]: logProb is -Inf.
#warning: logProb of data node y[296, 2]: logProb is -Inf.
#warning: logProb of data node y[323, 2]: logProb is -Inf.
#warning: logProb of data node y[333, 2]: logProb is -Inf.
#warning: logProb of data node y[334, 2]: logProb is -Inf.
#warning: logProb of data node y[335, 2]: logProb is -Inf.
#warning: logProb of data node y[336, 2]: logProb is -Inf.
#warning: logProb of data node y[337, 2]: logProb is -Inf.
#warning: logProb of data node y[338, 2]: logProb is -Inf.
#warning: logProb of data node y[339, 2]: logProb is -Inf.
#warning: logProb of data node y[340, 2]: logProb is -Inf.
#warning: logProb of data node y[341, 2]: logProb is -Inf.
#warning: logProb of data node y[342, 2]: logProb is -Inf.
#warning: logProb of data node y[343, 2]: logProb is -Inf.
#warning: logProb of data node y[346, 2]: logProb is -Inf.
#warning: logProb of data node y[347, 2]: logProb is -Inf.
#warning: logProb of data node y[348, 2]: logProb is -Inf.
#warning: logProb of data node y[350, 2]: logProb is -Inf.
#warning: logProb of data node y[351, 2]: logProb is -Inf.
#warning: logProb of data node y[352, 2]: logProb is -Inf.
#warning: logProb of data node y[353, 2]: logProb is -Inf.
#warning: logProb of data node y[355, 2]: logProb is -Inf.
#warning: logProb of data node y[356, 2]: logProb is -Inf.
#warning: logProb of data node y[357, 2]: logProb is -Inf.
#warning: logProb of data node y[358, 2]: logProb is -Inf.
#warning: logProb of data node y[482, 2]: logProb is -Inf.
#warning: logProb of data node y[486, 2]: logProb is -Inf.
#warning: logProb of data node y[491, 2]: logProb is -Inf


# Clear workspace ---------------------------------------------------------

rm(list=ls()) 
Sys.setenv(LANG = "en") ##change language of error messages to english

# Load libraries, functions etc -------------------------------------------
library(lubridate)
library(tidyverse)
library(dplyr)
filter<-dplyr::filter
select<-dplyr::select
library(coda)
library(MCMCvis)
library(nimble)
library(tictoc)
library(basicMCMCplots) # for trace plots called chainsPlot

## set root folder for project

setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Small projects/Woodcocks Steffen/woodcock")
#setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")

# Load data ---------------------------------------------------------------
# prepared in script WOCO_survival_data_compilation.R
load("data/woco_mig_input.RData")
nimbleOptions(allowDynamicIndexing = TRUE)

#Let's see what we have here

length(f.telemetry) #498 individuals
dim(effort) #A covariate for detection probability, linking to field effort
length(week) #1:24, sequence indicating the number of the week
nweeks #I guess this will be used in the loops (n.occasions)
tag #Indicating which individual is tagged
nind #Num. of individuals
nyears #Ah, I thought data were only for one year
y.telemetry #Actually, can NAs be specified at y here? On the one hand, values before first encounter should not intervene in the likelihood calculation. On the other hand, we usually add the event value for "missed".
#Why are there individuals that start with a 5 in the event matrix, and why do others start with NA?

#In the encounter history matrix, individuals not first encountered at t = 1 always have their first detection with evend code 1, but individuals first encountered at occ. 1 start with an event code 5. Why is that? For instance, individuals 12, 14, and 21
f.telemetry[12]
#Some individuals start with a 1, some with a 5, some with a 2.
y.telemetry[,1]
#While these values mean...
#1: Bird captured or recorded alive inside the study ara
#2: Bird recorded via transmitter outside the study area. How can this be true as a first encounter? Did the individuals get tagged outside of it? Well that is not possible - the individual must be Swiss, so if it is first encountered outside of Switzerland, the "Swiss nationality" cannot be assigned.
#5: Not observed / No signal. If the individual is first encountered at t=1, but the event is "not seen", that must lead to an error somewhere!

#One problem is here for sure: The likelihood does not allow first encounters to be anywhere else than state "Alive in the study area".
# Likelihood 
#for (i in 1:nind){
  # Define latent state at first capture
#  z[i,f[i]] <- 1 ## alive when first marked

#That initial state is inconsistent with being observed abroad.

#Let's check not only y = 1, but the state at which every individual is first encountered.

y.first <- numeric(length = length(f.telemetry))
for(i in 1:length(y.first))
  y.first[i] <- y.telemetry[i, f.telemetry[i]]

y.first[i]
#There's a lot of individuals starting at events other than 1.

#let's set all individuals as 1 at t = 1

for(i in 1:length(f.telemetry))
  y.telemetry[i, f.telemetry[i]] <- 1

#Do a check
y.first2 <- numeric(length = length(f.telemetry))
for(i in 1:length(y.first2))
  y.first2[i] <- y.telemetry[i, f.telemetry[i]]


y.first2
#Everything looks fine now

#I do not know whether NAs in the encounter histories could lead to an error, but removing them will for sure remove the number of NAs when assembling the model with nimbleModel, and will make error tracing easier.

for(i in 1:nrow(y.telemetry)){
  
  if(f.telemetry[i] > 1){
    
    y.telemetry[i, 1:(f.telemetry[i]-1)] <- 5
    
  }
  
}

#I modified the formulation of the model to add an absorbing state, it is probably more correct because otherwise a dead observation is always under the estimation of a recovery probability, which is probably not what we want. Thus, let's also recode the observation matrix

ww3 <- which(y.telemetry == 3, arr.ind = T)

for(i in 1:nrow(ww3)){
  
  if(ww3[i,2] != nweeks){
    
   y.telemetry[ww3[i,1], (ww3[i,2] + 1):nweeks] <- 5
    
  }
  
}

ww4 <- which(y.telemetry == 4, arr.ind = T)

for(i in 1:nrow(ww4)){
  
  if(ww4[i,2] != nweeks){
    
    y.telemetry[ww4[i,1], (ww4[i,2] + 1):nweeks] <- 5
    
  }
  
}

#Also, there are 2 individuals that are resighted in the study area after having moved away. That is not supported in our state-transition matrices. Let's fix it.

y.telemetry[82,c(3,12)] <- 5
y.telemetry[112, 17] <- 5



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
  # 4 dead outside study area (shot after after migration)
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
      
      logit(mig[i,t]) <- lm.mean +      ### intercept for mean survival
        b.mig.week*(week[t])     ### migration probability dependent on week
      
      logit(p.obs.in[i,t]) <- lpin.mean +      ### intercept for mean survival
        b.obs.effort*(effort[i,t]) +    ### observation dependent on effort in that week and year
        b.obs.tag*(tag[i])    ### observation dependent on whether animal was tagged
      
      logit(p.obs.out[i,t]) <- lpout.mean +      ### intercept for mean survival
        b.obs.tag*(tag[i])    ### observation dependent on whether animal was tagged
      
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
  mean.mig ~ dunif(0,1)   # fairly uninformative prior for weekly survival probabilities
  lm.mean <- log(mean.mig/(1 - mean.mig))    # logit transformed migration intercept
  
  #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
  mean.p.in ~ dunif(0,1)   # fairly uninformative prior for weekly detection probabilities
  lpin.mean <- log(mean.p.in/(1 - mean.p.in))    # logit transformed detection intercept
  
  mean.p.out ~ dunif(0,0.3)   # fairly uninformative prior for weekly detection probabilities
  lpout.mean <- log(mean.p.out/(1 - mean.p.out))    # logit transformed detection intercept
  
  #### SLOPE PARAMETERS FOR PROBABILITIES ON LOGIT SCALE
  b.mig.week ~ dunif(0,2)         # Prior for week effect on migration probability on logit scale - must be positive
  b.obs.effort ~ dunif(0, 2)         # Prior for effort effect on observation probability on logit scale  - must be positive
  b.obs.tag ~ dunif(0, 2)         # Prior for tag effect on observation probability on logit scale  - must be positive 
  
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
    z[i,f[i]] <- 1 ## alive when first marked
    for (t in (f[i]+1):nweeks){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:5])     ### check different CJS distribution in nimbleEcology dcjs
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])
      rep.states[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])    ## include this as goodness of fit measure
    } #t
  } #i
}) ## end of nimble code chunk


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE NIMBLE RUN WITH DATA INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters monitored
parameters.telemetry <- c("mean.mig","b.obs.effort","b.obs.tag","b.mig.week")

# Initial values for some parameters
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

####I see the z inits also have some issues here. Let's rebuild them

z.telemetry2 <- matrix(nrow = nrow(y.telemetry), ncol = ncol(y.telemetry), data = NA)

for(i in 1:nrow(z.telemetry2)){
  
  z.telemetry2[i,f.telemetry[i]] <- 1
  
}

#State 2: dead inside study area. That means, on all occasions before the individual had to be alive and in the study area. Caution: it matches with event 3
w2 <- which(y.telemetry == 3, arr.ind = T) 

for(i in 1:nrow(w2)){
  
  z.telemetry2[w2[i,1], (f.telemetry[w2[i,1]]+1):w2[i,2]] <- 1
  z.telemetry2[w2[i,1], w2[i,2]] <- 2
  
  if(w2[i,2] != nweeks){
    
    z.telemetry2[w2[i,1], (w2[i,2]+1):nweeks] <- 5
    
  }
  
  
}

#State 3:  Alive outside study area. Matches with event 2. Once the individual is outside, it cannot go back

w3 <- which(y.telemetry == 2, arr.ind = T)

for(i in 1:nrow(w3)){
  
  z.telemetry2[w3[i,1], w3[i,2]:nweeks] <- 3
  
}

#State 4: Dead outside study area. Matches with event 4 this time.

w4 <- which(y.telemetry == 4, arr.ind = T)

for(i in 1:nrow(w4)){
  
  #Was the individual resighted before abroad?
  tmp <- which(z.telemetry2[w4[i,1],] == 3)
  
  if(length(tmp)>0){ #If yes,
    
    z.telemetry2[w4[i,1], tmp:w4[i,1]] <- 3 #Fill with 3s from first detection abroad
    
  }
  
  #And add a 4 to the dead encounter
  
  z.telemetry2[w4[i,1], w4[i,2]] <- 4
  
  #And fill next occasions with 5s (absorbing states)
  
  if(w4[i,2] != nweeks){
    
    z.telemetry2[w4[i,1], (w4[i,2]+1):nweeks] <- 5
    
  }
  
}

#The rest of values after first encounter should be 1s, if I am right

for(i in 1:nrow(z.telemetry2))
  for(t in f.telemetry[i]:nweeks)
    if(is.na(z.telemetry2[i,t]))
      z.telemetry2[i,t] <- 1


#Let's do some checks. Without a super global check things look fine.
y.telemetry[1:5,]
z.telemetry2[1:5,]

# Missing values (NAs) or non-finite values were found in model variables: phi, b.phi.age, b.phi.mig, b.phi.feed, b.phi.sex, b.phi.FRA, b.phi.SUI, b.phi.ESP, p.obs, tag.fail, p.found.dead, ps, po, z, y, rep.states.
smartInit1 <- list(z = z.telemetry2,
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
                   mean.mig = 0.5,   # fairly uninformative prior for weekly survival probabilities
                   
                   #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
                   mean.p.in = 0.2,   # fairly uninformative prior for weekly detection probabilities
                   mean.p.out = 0.05,   # fairly uninformative prior for weekly detection probabilities
                   
                   #### SLOPE PARAMETERS FOR PROBABILITIES ON LOGIT SCALE
                   b.mig.week = 0.5,         # Prior for week effect on migration probability on logit scale - must be positive
                   b.obs.effort = 0.5,         # Prior for effort effect on observation probability on logit scale  - must be positive
                   b.obs.tag = 0.5         # Prior for tag effect on observation probability on logit scale  - must be positive 
                   
)

# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.mig.model,
                    constants=telemetry.constants,
                    data = telemetry.data,
                    inits = smartInit1,
                    calculate=TRUE)

test$initializeInfo()

test$calculate("z")

#I did some error tracking in the model to address former fails. I leave it here if that helps in further problems. 

#In the first case, I had NAs in my z latent state matrix, I saw where the problems were this way:

znames <- matrix(nrow = nrow(y.telemetry), ncol = ncol(y.telemetry))
for(i in 1:nrow(y.telemetry))
  for(t in 1:ncol(y.telemetry))
    znames[i,t] <-paste0("z[", i, ",", t,"]")

fails <- matrix(nrow = nrow(znames), ncol = ncol(znames))

for(i in 1:nrow(fails))
  for(t in 1:ncol(fails))
    fails[i,t] <- test$calculate(znames[i,t])

#4 encounter histories with errors
#fails[10,] #Starts at t = 18
#y.telemetry[10,] #The individual is found dead abroad on the next occasion
#z.telemetry2[10,] #I see the problem: y keeps being 4 after a death

#fails[12,]
#Again same issue: after the individuals die, we retain the dead value, we have to send it to unseen
#y.telemetry[12,]
#z.telemetry2[12,]

#fails[14,]
#y.telemetry[14,]
#z.telemetry2[14,]
#fails[16,]

#What about dead and recovered in study area
#y.telemetry[338,] #Same issue
#z.telemetry2[338,] #Inconsistent here, also the inits do not look nice at all
#fails[338,] #But interestingly there is no -Inf in the likelihood calculation. However this needs to be modified

#####IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#####

#More fails right now: some y values are inconsistent with the model's transitions and event matrices. 
#A short summary of some of them: always at t = 2 and for individuals tagged at t = 1
#warning: logProb of data node y[32, 2]: logProb is -Inf.
#warning: logProb of data node y[66, 2]: logProb is -Inf.
#warning: logProb of data node y[67, 2]: logProb is -Inf.
#warning: logProb of data node y[73, 2]: logProb is -Inf.
#warning: logProb of data node y[74, 2]: logProb is -Inf.
#warning: logProb of data node y[77, 2]: logProb is -Inf.

#Check some of them
#Individual 32
y.telemetry[32,]
z.telemetry2[32,]
test$calculate("y[32,]")
test$po[1,32,1,5] #Event code 5 at t = 2 is not allowed (prob = 0) according to our event matrix 

#Individual 66: again same problem
y.telemetry[66,]
test$calculate("y[66,]")
test$po[1,66,1,5]

#There's something weird happening at t=2, and it is always happening from t = 1 to t = 2. It is hard to grasp if one does not know well the model. Event 5 means "unobserved", and the observation model is complex: there is an effort covariate, which takes value 0 at t = 2, but at the same time there are tagged and untagged birds (I guess that refers to the GPS).

#Ind. 35: Another individual that is first encountered at 0, but where the likelihood can be calculated appropriately

y.telemetry[35,]
test$calculate("y[35,]") #This is right
#155555555555555555555552
test$po[1,35,1,5] #P = 1.223608e-06 (the event CAN happen, p > 0)
tag[35] #The individual is tagged

#Ind. 32: again with problems
test$calculate("y[32,]") #This is not
y.telemetry[32,]
#151555155555155555555555 
test$po[1,32,1,5] #P = 0, leads to the Inf in the likelihood calculation
tag[32] #The individual is untagged

test$calculate("y[66,]") #-Inf
tag[66] #Again, the individual is untagged

test$calculate("y[75,]")
y.telemetry[75,]#Not Inf
y.telemetry[73,]#Inf
y.telemetry[74,]#Inf

tag[75] #75 was tagged
tag[c(73,74)] #73 and 74 were untagged

#In our complex observation model, it appears that whenever the individual is untagged, it must be observed at t = 2. instead, if the individual is tagged, it can be missed at t=2. This looks strange. Probably a model misspecification?

#Let's test whether the model works
tic()

ctest <- compileNimble(test)

testConf <- configureMCMC(test, print = TRUE)

testConf$addMonitors(parameters.telemetry)
testMCMC <- buildMCMC(testConf)
cmodel <- compileNimble(testMCMC, project = test) 

cmodel$run(10000) #So far the model throws the errors with the ys
samples <- as.matrix(cmodel$mvSamples)


