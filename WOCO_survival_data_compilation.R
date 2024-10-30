###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS DATA PREPARATION--------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds

## QUESTIONS TO MURDY:
# - what does 'censor' mean?
# - what temporal resolution do we want (weekly ok)?
# - extent of season (Sept - Nov), what is the intermittent dip in effort?
# do we have a metric of observation effort?
# what are age categories 1-8? EURING Code?


# Clear workspace ---------------------------------------------------------

rm(list=ls()) 
Sys.setenv(LANG = "en") ##change language of error messages to english

# Load libraries, functions etc -------------------------------------------
library(lubridate)
library(tidyverse)
library(readxl)
library(dplyr)
filter<-dplyr::filter
select<-dplyr::select

## set root folder for project
#setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")

# Load data ---------------------------------------------------------------

woco <- read_excel("data/Daten_Waldschnepfe.xlsx", sheet="Data") %>%
  rename(age=`Alter (Euring)`)
woco

### CREATE A TIME SERIES DATA FRAME ###
mindate<-min(woco$Datum)
maxdate<-max(woco$Datum)
timeseries<-data.frame(date=seq(mindate, maxdate, "1 month")) %>%
  mutate(month=month(date),year=year(date)) %>%
  mutate(ymo=as.character(format(date, format="%Y_%m"))) %>%
  mutate(season=ifelse(month %in% c(4,5,6,7,8,9),"summer","winter")) %>%
  mutate(col=seq_along(date))
dim(timeseries)




# Inspect data ---------------------------------------------------------------
table(woco$Ort)
table(woco$age)
hist(month(woco$Datum))
table(woco$Beobachtung)
table(woco$Markierung)
table(woco$Censor)
length(unique(woco$Ring_num))
woco %>% filter(month(Datum) %in% c(9,10,11)) %>%
  mutate(jday=yday(Datum)) %>%
  ggplot(aes(x=jday)) + geom_histogram()

table(woco$Beobachtung, woco$Ort)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PREPARATION APPROACH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## for each individual and year, create weekly occasions for Sept to Nov
## use data from following year to determine whether bird survived or not
# True States (S) - these are often unknown and cannot be observed
# 1 dead
# 2 alive in study area
# 3 alive outside study area (after migration)

  # Observed States (O) - these are based on the actual transmission history
  # 1 Bird (re)captured alive
  # 2 Bird recorded via transmitter in study area
  # 3 Bird recorded via transmitter outside study area
  # 4 Dead bird recovered  in study area
  # 5 Dead bird recovered  outside study area
  # 6 No signal (=not observed)
  # 7 Tag lost and recovered (we will leave this out because there is only a single instance)

### inspect individual that lost the logger (happened 3 days after last live location - does not yield additional info)
woco %>% filter(Beobachtung=="Senderfund")
woco %>% filter(Ring_num==112985) %>% arrange(Datum)

### inspect individual that was caught in France
woco %>% filter(Beobachtung=="Fang") %>% filter(Ort!="UG")
woco %>% filter(Ring_num==115155) %>% arrange(Datum)


# CREATE ANNUAL ENCOUNTER HISTORY  ---------------------------------------------------------------
woco_ch <- woco %>%
  mutate(year=year(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,year,Datum,Beobachtung,Markierung,Ort) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  mutate(obs=1) %>%
  group_by(Ring_num, year) %>%
  summarise(N=sum(obs)) %>%
  spread(key=year, value=N, fill=0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE ANNUAL CAPTURE HISTORIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### weekly encounter history ####
woco_ann_ch<-woco %>%
  mutate(year=year(Datum), week=paste('wk',week(Datum),"")) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  filter(month(Datum) >8) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  mutate(OS=ifelse(Beobachtung=="Fang",1,
                   ifelse((Beobachtung=="Senderlokalisation" & Ort=="UG"),2,
                          ifelse((Beobachtung=="Senderlokalisation" & Ort!="UG"),3,
                                 ifelse((Beobachtung=="Totfund" & Ort=="UG"),4,
                                        ifelse((Beobachtung=="Totfund" & Ort!="UG"),5,6)))))) %>%
  group_by(id, week) %>%
  summarise(N=max(OS)) %>%
  spread(key=week, value=N, fill=6) %>%
  mutate(wk54=6)  ## create blank column for records past the migration season


### CREATE BLANK MATRICES TO HOLD INFORMATION ABOUT TRUE AND OBSERVED STATES ###

woco.obs.matrix<-woco_ann_ch %>% arrange(id)
woco.state.matrix<-woco_ann_ch %>% arrange(id)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### FILL MATRICES WITH STATE INFORMATION ###
for(i in 1:nrow(woco.obs.matrix$id)){
  
  ### extract ind and year
  yr=separate_wider_delim(woco_ann_ch[i,],cols="id",delim="_", names=c('ring','year'))$year
  rn=separate_wider_delim(woco_ann_ch[i,],cols="id",delim="_", names=c('ring','year'))$ring
  
  ### extract start and end dates for selected bird
  daterange<-woco %>%
    filter(Ring_num==rn) %>%
    filter(Beobachtung!="Senderfund") %>%
    filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
    summarise(first=min(Datum),last=max(Datum))
  
  ### specify columns in matrix to populate
  startcol<-min(which(woco.obs.matrix[i,2:20]!=6))
  if(year(daterange$first)<yr) {startcol<-2}
  stopcol<-max(which(woco.obs.matrix[i,2:20]<4))
  if(year(daterange$last)>yr) {stopcol<-20}
  
  ## ASSIGN OBSERVED STATE

  woco.obs.matrix[woco.obs.matrix$individual_id==i,startcol:(stopcol-1)]<-ifelse(CH[CH$individual_id==i,startcol:(stopcol-1)]==1,1,5)  ## if there are no data then label it with 5
  if(startcol==stopcol){woco.obs.matrix[woco.obs.matrix$individual_id==i,2:(stopcol-1)]<-NA} ## for the few cases where stopcol-1 is actually before startcol
  
  ### add the final observation and extend to end of time series
  woco.obs.matrix[woco.obs.matrix$individual_id==i,stopcol:(max(timeseries$col)+1)]<-wocoinds$OS[wocoinds$individual.local.identifier==i]   ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!
  
  ## ASSIGN INITIAL TRUE STATE (to initialise z-matrix of model)
  ### this needs careful manipulation to appropriately configure intermediate 0s
  woco.state.matrix[woco.state.matrix$individual_id==i,(startcol+1):(stopcol-1)]<-2      ## state at first capture is known, hence must be NA in z-matrix
  woco.state.matrix[woco.state.matrix$individual_id==i,stopcol:(max(timeseries$col)+1)]<-wocoinds$TS[wocoinds$individual.local.identifier==i]   ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!
  woco.state.matrix[woco.state.matrix$individual_id==i,2:startcol]<-NA ## error occurs if z at first occ is not NA, so we need to specify that for birds alive for <1 month because stopcol-1 = startcol
  

}

## CHECK FOR NA IN ANY OF THE MATRICES -----------------
## none of these should contain any NA because it will fail the model
any(is.na(woco.dep.matrix))

## these will contain NA before individuals enter the time series
any(is.na(woco.age.matrix))
any(is.na(woco.obs.matrix))
any(is.na(woco.state.matrix))

# which(woco.obs.matrix==5, arr.ind=T)
# woco.obs.matrix[57,20:40]


#### Convert to numeric matrices that NIMBLE can loop over
y.telemetry<-as.matrix(woco.obs.matrix[,2:(max(timeseries$col)+1)])
z.telemetry<-as.matrix(woco.state.matrix[,2:(max(timeseries$col)+1)])
age.mat<-as.matrix(woco.age.matrix[,2:(max(timeseries$col)+1)])
dep.mat<-as.matrix(woco.dep.matrix[,2:(max(timeseries$col)+1)])

#### SCALE THE AGE SO THAT NUMERICAL OVERFLOW DOES NOT OCCUR
max(age.mat, na.rm=T)
age.mat<-ifelse(age.mat>48,48,age.mat)
agescale<-scale(1:max(age.mat, na.rm=T))
#feedscale<-scale(1:max(feed.mat, na.rm=T))
dim(age.mat)
dim(dep.mat)
dim(y.telemetry)
dim(z.telemetry)

#### create vector of first marking and of last alive record
get.first.telemetry<-function(x)min(which(!is.na(x)))
get.last.telemetry<-function(x)max(which(!is.na(x) & x==1))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)
l.telemetry<-apply(y.telemetry,1,get.last.telemetry)

#### CREATE VECTOR OF PFDP
pfdp.vec<-wocoinds %>%
  mutate(PFPD=ifelse(is.na(PFPD),0,PFPD)) %>%
  arrange(individual.local.identifier) %>%
  select(PFPD)
pfdp.scale<-as.numeric(scale(as.numeric(pfdp.vec$PFPD)))


############# SAVE AND RE-LOAD PREPARED DATA ----------
save.image("data/woco_surv_tradeoff_input.RData")
load("data/woco_surv_tradeoff_input.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP NIMBLE CODE FOR SURVIVAL MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify model in NIMBLE format
surv.model<-nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi: seasonal survival probability intercept
  # tag.fail: probability that tag failed or was lost
  
  # probability to be tracked with functioning tag is assumed to be perfect (=1) and not included as parameter
  # p.obs: probability to be observed alive despite the tag being defunct
  # p.found.dead: probability for carcass to be recovered
  
  # -------------------------------------------------
  # States (S):
  # 1 dead
  # 2 alive with functioning tag
  # 3 alive with defunct tag or tag lost
  
  # Observations (O):
  # 1 Tag ok, bird moving
  # 2 Tag ok, bird not moving (dead, or tag lost and no longer on bird)
  # 3 Tag failed, bird observed alive
  # 4 Dead bird recovered
  # 5 No signal (=not seen, status unknown)
  
  # -------------------------------------------------
  
  # Priors and constraints
  
  
  #### MONTHLY SURVIVAL PROBABILITY
  for (i in 1:nind){
    for (t in f[i]:(n.occasions)){
      logit(phi[i,t]) <- lp.mean +      ### intercept for mean survival
        b.phi.season*(season[t]) +     ### survival dependent on season
        b.phi.pfdp*(pfdp[i])*(1-dep[i,t]) +     ### survival dependent on PFDP only for months after dependence
        b.phi.age*(age[i,t])    ### survival dependent on age (in years)

    } #t
  } #i
  
  #### BASELINE FOR SURVIVAL PROBABILITY (varies by season)
  #for(s in 1:2) {
    mean.phi ~ dbeta(15, 2)   # fairly uninformative prior intercept for all monthly survival probabilities
    lp.mean <- log(mean.phi/(1 - mean.phi))    # logit transformed survival intercept
  #}
  
  
  #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
  b.phi.age ~ dnorm(0, tau=0.001)           # Prior for age effect on survival probability on logit scale
  b.phi.season ~ dnorm(0, tau=0.001)           # Prior for season effect on survival probability on logit scale
  b.phi.pfdp ~ dunif(-2,2)         # Prior for migration effect on survival probability on logit scale
  #b.phi.sex ~ dnorm(0, tau=0.001)         # Prior for sex effect on survival probability on logit scale 
 
  
  #### TAG FAILURE AND LOSS PROBABILITY
  for (i in 1:nind){
    p.obs.ind[i] ~ dunif(0,1)
    tag.fail.ind[i] ~ dunif(0,1)
    p.found.dead.ind[i] ~ dunif(0,1)
    
    for (t in f[i]:(n.occasions)){
      # logit(p.obs[i,t]) <- base.obs
      # logit(tag.fail[i,t]) <- base.fail
      # logit(p.found.dead[i,t]) <- base.recover
      p.obs[i,t] <-p.obs.ind[i]
      tag.fail[i,t] <- tag.fail.ind[i]
      p.found.dead[i,t] <- p.found.dead.ind[i]
    } #t
  } #i
  
  
  ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
  # base.obs ~ dnorm(-2, tau=0.001)                # Prior for observation probability of birds with malfunctioning transmitter on logit scale
  # base.fail ~ dnorm(-1, tau=0.001)               # Prior for tag failure probability on logit scale
  # base.recover ~ dnorm(0.5, tau=0.001)            # Prior for probability to recover dead bird on logit scale
  
  
  # -------------------------------------------------
  # Define state-transition and observation matrices 
  # -------------------------------------------------
  
  for (i in 1:nind){
    
    for (t in f[i]:(n.occasions-1)){
      
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      
      ps[1,i,t,1]<-1    ## dead birds stay dead
      ps[1,i,t,2]<-0
      ps[1,i,t,3]<-0
      
      ps[2,i,t,1]<-(1-phi[i,t])
      ps[2,i,t,2]<-phi[i,t] * (1-tag.fail[i,t])
      ps[2,i,t,3]<-phi[i,t] * tag.fail[i,t]
      
      ps[3,i,t,1]<-(1-phi[i,t])
      ps[3,i,t,2]<-0
      ps[3,i,t,3]<-phi[i,t]
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      
      po[1,i,t,1]<-0
      po[1,i,t,2]<-p.obs[i,t] * (1-tag.fail[i,t]) * (1-p.found.dead[i,t])
      po[1,i,t,3]<-0
      po[1,i,t,4]<-p.found.dead[i,t]
      po[1,i,t,5]<-(1-p.obs[i,t]) * tag.fail[i,t] * (1-p.found.dead[i,t])
      
      po[2,i,t,1]<-(1-tag.fail[i,t])
      po[2,i,t,2]<-0
      po[2,i,t,3]<-0
      po[2,i,t,4]<-0
      po[2,i,t,5]<-tag.fail[i,t]
      
      po[3,i,t,1]<-0
      po[3,i,t,2]<-0
      po[3,i,t,3]<-p.obs[i,t]
      po[3,i,t,4]<-0
      po[3,i,t,5]<-(1-p.obs[i,t])
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 2 ## alive when first marked
    for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:3])     ### check different CJS distribution in nimbleEcology dcjs
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
                            season=ifelse(timeseries$season=="summer",1,2),
                            age= ifelse(is.na(age.mat),0,scale(age.mat)),
                            #sex = ifelse(birds$sex_compiled=="f",0,1),
                            dep = dep.mat,  ## range(mig.mat)
                            pfdp = pfdp.scale,  ## single value for each bird
                            nind = dim(y.telemetry)[1],
                            n.occasions = dim(y.telemetry)[2])
telemetry.data <- list(y = y.telemetry)


### SPECIFY DIMENSIONS for arrays used with empty indices
telemetry.dims <- list(ps = c(3, dim(y.telemetry)[1], dim(y.telemetry)[2], 3),
                       po = c(3, dim(y.telemetry)[1], dim(y.telemetry)[2], 5))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY AND SET UP MODEL RUNS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters monitored
parameters.telemetry <- c("base.obs","base.fail","base.recover",
                          "mean.phi","lp.mean",
                          "b.phi.age","b.phi.season","b.phi.pfdp")

# Initial values for some parameters
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

# Missing values (NAs) or non-finite values were found in model variables: phi, b.phi.age, b.phi.mig, b.phi.feed, b.phi.sex, b.phi.FRA, b.phi.SUI, b.phi.ESP, p.obs, tag.fail, p.found.dead, ps, po, z, y, rep.states.
smartInit1 <- list(z = ifelse(is.na(z.telemetry),2,z.telemetry),
                   y=ifelse(is.na(y.telemetry),1,y.telemetry),
                   mean.phi = 0.95, ### mean survival over 6 months
                   lp.mean = logit(0.95),
                   #phi = as.matrix(CH[,2:(dim(y.telemetry)[2]+1)]),
                   # base.obs = rnorm(1,-2, 0.001),                # Start value for intercept of observation probability on logit scale
                   # base.fail = rnorm(1,-1, 0.001),               # Start value for intercept of tag failure probability on logit scale
                   # base.recover = rnorm(1,0.5, 0.001),
                   p.obs.ind = runif(dim(y.telemetry)[1],0, 1),                # Start value for intercept of observation probability on logit scale
                   tag.fail.ind = runif(dim(y.telemetry)[1],0, 1),               # Start value for intercept of tag failure probability on logit scale
                   p.found.dead.ind = runif(dim(y.telemetry)[1],0, 1),
                   b.phi.age = 0,           # Start value for age effect on survival probability on logit scale
                   b.phi.pfdp = 0,         # Start value for migration effect on survival probability on logit scale
                   b.phi.season = 0         # Start value for sex effect on survival probability on logit scale 
)



# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 2000
n.burnin <- 1000
n.chains <- 4


# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = surv.model,
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




