###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS --------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds

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
library(basicMCMCplots) # for trace plots called chainsPlot

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")

# Load data ---------------------------------------------------------------
# prepared in script WOCO_survival_data_compilation.R
load("data/woco_mig_input.RData")




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
          b.mig.week*(week[t])     ### survival dependent on season
      
      logit(p.obs.in[i,t]) <- lpin.mean +      ### intercept for mean survival
        b.obs.effort*(effort[year[i],week[t]]) +    ### observation dependent on effort in that week and year
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
      
      ps[2,i,t,1]<-0
      ps[2,i,t,2]<-1
      ps[2,i,t,3]<-0
      ps[2,i,t,4]<-0
      
      ps[3,i,t,1]<-0
      ps[3,i,t,2]<-0
      ps[3,i,t,3]<-phi[i,t]
      ps[3,i,t,4]<-(1-phi[i,t])
      
      ps[4,i,t,1]<-0
      ps[4,i,t,2]<-0
      ps[4,i,t,3]<-0
      ps[4,i,t,4]<-1
      
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
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1 ## alive when first marked
    for (t in (f[i]+1):nweeks){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:4])     ### check different CJS distribution in nimbleEcology dcjs
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])
      rep.states[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])    ## include this as goodness of fit measure
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
                            #effort = effort,  # removed because dynamic indexing is not allowed
                            year = year,
                            tag = tag,
                            nind = nind,
                            nweeks=nweeks,
                            nyears = nyears)
telemetry.data <- list(y = y.telemetry,
                       effort = effort)


### SPECIFY DIMENSIONS for arrays used with empty indices
telemetry.dims <- list(ps = c(3, dim(y.telemetry)[1], dim(y.telemetry)[2], 4),
                       po = c(3, dim(y.telemetry)[1], dim(y.telemetry)[2], 5))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters monitored
parameters.telemetry <- c("mean.mig","b.obs.effort","b.obs.tag","b.mig.week")

# Initial values for some parameters
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

# Missing values (NAs) or non-finite values were found in model variables: phi, b.phi.age, b.phi.mig, b.phi.feed, b.phi.sex, b.phi.FRA, b.phi.SUI, b.phi.ESP, p.obs, tag.fail, p.found.dead, ps, po, z, y, rep.states.
smartInit1 <- list(z = ifelse(is.na(z.telemetry),1,z.telemetry),
                   y=ifelse(is.na(y.telemetry),1,y.telemetry),
                   
                   mean.phi = rep(0.995,nyears),
                   mean.p.dead.in<-rep(0.5,nyears),
                   mean.p.dead.out<-rep(0.2,nyears),
                   
                   #### BASELINE FOR MIGRATION PROBABILITY (varies by year)
                   mean.mig <- 0.5,   # fairly uninformative prior for weekly survival probabilities

                   #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
                   mean.p.in <- 0.2,   # fairly uninformative prior for weekly detection probabilities
                   mean.p.out <- 0.05,   # fairly uninformative prior for weekly detection probabilities
                   
                   #### SLOPE PARAMETERS FOR PROBABILITIES ON LOGIT SCALE
                   b.mig.week <- 0.5,         # Prior for week effect on migration probability on logit scale - must be positive
                   b.obs.effort <- 0.5,         # Prior for effort effect on observation probability on logit scale  - must be positive
                   b.obs.tag <- 0.5         # Prior for tag effect on observation probability on logit scale  - must be positive 
                   
)



# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 2000
n.burnin <- 1000
n.chains <- 4


# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.mig.model,
                    constants=telemetry.constants,
                    data = telemetry.data,
                    inits = smartInit1,
                    calculate=TRUE)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
test$initializeInfo()
#help(modelInitialization)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits


# use test output as starting values or check where the NA comes from
test$logProb_b.phi.age
test$logProb_mean.phi
test$b.phi.season
test$mean.phi
test$logProb_y ### there should not be any -Inf in this matrix




### FIT FULL MODEL TO ALL DATA ############################
# RUN NIMBLE MODEL --------------------

tic()
### this takes 40-50 mins for 20000 iterations and converges in that time
woco_surv <- nimbleMCMC(code = surv.model,
                        constants=telemetry.constants,
                        data = telemetry.data,
                        inits = smartInit1,
                        dimensions=telemetry.dims,
                        monitors = parameters.telemetry,
                        thin=3,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains,
                        progressBar = getNimbleOption("MCMCprogressBar"),
                        summary=T)
toc() 


saveRDS(woco_surv,"output/woco_seasonal_survival_output_nimble.rds")
#woco_surv<-readRDS("output/woco_seasonal_survival_output_nimble.rds")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMINE OUTPUT AND DIAGNOSTICS WITH MCMCvis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### MODEL ASSESSMENT ####
out<- as.data.frame(MCMCsummary(woco_surv$samples, params=c("mean.phi","b.phi.age","b.phi.sex","b.phi.mig","b.phi.feed","b.phi.FRA","b.phi.SUI","b.phi.ESP")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_telemetry_seasonal_surv_parm.csv")
#out<-fread("output/woco_telemetry_surv_parm_v1.csv")


MCMCplot(woco_surv$samples, params=c("mean.phi","b.phi.age","b.phi.sex","b.phi.mig","b.phi.feed","b.phi.FRA","b.phi.SUI","b.phi.ESP"))
ggsave("output/woco_seasonal_survival_parameter_estimates.jpg", height=9, width=9)
# ggsave("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/wocopopmod/output/Fig_S1_parameter_estimates.jpg", height=11, width=8)

## look at the chains and whether they mixed well
chainsPlot(woco_surv$samples,
           var=c("base.obs","base.fail","base.recover","b.recover.year","b.obs.FRA","b.obs.SUI"
           ))
chainsPlot(woco_surv$samples,
           var=c("mean.phi","lp.mean",
                 "b.phi.age","b.phi.sex","b.phi.mig","b.phi.feed","b.phi.FRA","b.phi.SUI","b.phi.ESP","b.phi.season"
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
############ PRODUCE ANNUAL SURVIVAL OUTPUT GRAPH SHOWING DIFFERENCES BETWEEN MALES AND FEMALES  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### PREPARE RAW MCMC OUTPUT
parmcols<-dimnames(woco_surv$samples[[1]])[[2]]

# ### COMBINE SAMPLES ACROSS CHAINS
MCMCout<-rbind(woco_surv$samples[[1]],woco_surv$samples[[2]],woco_surv$samples[[3]],woco_surv$samples[[4]])
# str(MCMCout)

### SET UP ANNUAL TABLE FOR PLOTTING THE ANNUAL SURVIVAL GRAPH
## uses mig to also specify time in FRA and ESP

AnnTab<-expand.grid(age=seq(1,6,0.5), sex=c(0,1), feed = c(0,1), mig=c(0,1)) %>%
  #mutate(mig=0.5) %>% ### migration has such a weak effect that we don't include it - set to intermediate probability to migrate
  mutate(scaleage=(age-attr(agescale, 'scaled:center'))/attr(agescale, 'scaled:scale')) %>%
  mutate(season=ifelse(nchar(age)==3,"non-breeding","breeding")) %>%
  mutate(SUI=ifelse(season=="breeding",1,
                    ifelse(mig==0,1,0.2))) %>%
  mutate(FRA=ifelse(season=="breeding",0,
                    ifelse(mig==0,0,0.3))) %>%
  mutate(ESP=ifelse(season=="breeding",0,
                    ifelse(mig==0,0,0.5)))

### CALCULATE PREDICTED VALUE FOR EACH SAMPLE

MCMCpred<-data.frame()
for(s in 1:nrow(MCMCout)) {
  
  X<-  AnnTab %>%
    
    ##CALCULATE MONTHLY SURVIVAL
    mutate(logit.surv=ifelse(season=="breeding",as.numeric(MCMCout[s,match("lp.mean[1]",parmcols)]),as.numeric(MCMCout[s,match("lp.mean[2]",parmcols)]))+
             as.numeric(MCMCout[s,match("b.phi.age",parmcols)])*scaleage +
             as.numeric(MCMCout[s,match("b.phi.mig",parmcols)])*mig +
             as.numeric(MCMCout[s,match("b.phi.sex",parmcols)])*sex +
             as.numeric(MCMCout[s,match("b.phi.feed",parmcols)])*feed*SUI +
             as.numeric(MCMCout[s,match("b.phi.SUI",parmcols)])*SUI +
             as.numeric(MCMCout[s,match("b.phi.FRA",parmcols)])*FRA +
             as.numeric(MCMCout[s,match("b.phi.ESP",parmcols)])*ESP)  %>%
    
    ##BACKTRANSFORM TO NORMAL SCALE
    mutate(surv=plogis(logit.surv)) %>%
    mutate(simul=s)              
  
  MCMCpred<-rbind(MCMCpred,as.data.frame(X)) 
}

### CREATE PLOT

FIGURE<-MCMCpred %>% rename(raw.surv=surv) %>% dplyr::select(-logit.surv,-simul,-scaleage) %>%
  #gather(key="Country", value="time", -age,-sex,-feed,-raw.surv,-mig) %>%
  #filter(time>0) %>%
  mutate(age=as.integer(age)) %>% ## combine the two seasons into the same year
  group_by(sex,feed,age,mig) %>%
  summarise(surv=quantile(raw.surv,0.5),surv.lcl=quantile(raw.surv,0.025),surv.ucl=quantile(raw.surv,0.975)) %>%
  mutate(feed=ifelse(feed==1,"depends on human food","avoids human feeding")) %>%
  mutate(sex=ifelse(sex==1,"male","female")) %>%
  mutate(mig=ifelse(mig==1,"migrant","resident")) %>%
  
  
  ggplot()+
  geom_ribbon(aes(x=age, ymin=surv.lcl, ymax=surv.ucl, fill=sex, linetype=feed), alpha=0.2) +   ##
  geom_line(aes(x=age, y=surv, color=sex, linetype=feed),linewidth=1)+     ##
  facet_wrap(~mig) +
  
  ## format axis ticks
  scale_x_continuous(name="Age in years", limits=c(1,6), breaks=seq(1,6,1), labels=seq(1,6,1)) +
  #scale_y_continuous(name="Monthly survival probability", limits=c(0.8,1), breaks=seq(0.,1,0.05)) +
  labs(y="Annual survival probability",fill="Sex", colour="Sex", type="Feeding") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title=element_text(size=16, color="black"),
        legend.position=c(0.8,0.1), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
FIGURE

ggsave(plot=FIGURE,
       filename="output/woco_seasonal_survival_mig.jpg", 
       device="jpg",width=11, height=8)




