###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS DATA PREPARATION--------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds

## QUESTIONS TO MURDY:
# what does 'censor' mean? check with Michi if censoring necessary for migration
# what temporal resolution do we want (weekly ok)? yes
# extent of season (Sept - Nov), what is the intermittent dip in effort? monthly capture events
# do we have a metric of observation effort?
# what are age categories 1-8? EURING Code? Only needed if age affects departure or survival
# annual variation in survival, detection, migration?


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
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")

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
woco %>% filter(month(Datum) %in% c(8,9,10,11)) %>%
  mutate(jday=yday(Datum)) %>%
  ggplot(aes(x=jday)) + geom_histogram()

table(woco$Beobachtung, woco$Ort)


## check last observation per bird in study area
woco %>% filter(month(Datum) %in% c(8,9,10,11)) %>%
  filter(Ort=="UG") %>%
  mutate(year=year(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  mutate(jday=yday(Datum)) %>%
  mutate(date=ymd("2023-12-31")+days(jday)) %>%
  group_by(id) %>%
  summarise(last=max(date)) %>%
  ggplot(aes(x=last)) +
  geom_histogram(aes(x=last,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=7)+                               
  scale_x_date(name="last obs in study area",date_breaks="1 week", date_labels="%d-%b")+
  labs(y="proportion of woodcocks") +
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=12, color="black",angle=45,hjust = 1),
        axis.text.y=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
ggsave("output/prop_woodcock_per_week_rawdat.jpg", width=7, height=6)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PREPARATION APPROACH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## for each individual and year, create weekly occasions for Sept to Nov
## use data from following year to determine whether bird survived or not
# True States (S) - these are often unknown and cannot be observed
# 1 alive in study area
# 2 dead in study area
# 3 alive outside study area (after migration)
# 4 dead outside study area (shot after after migration)

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
woco_ann_ch_true<-woco %>%
  mutate(year=year(Datum), week=paste('wk',week(Datum),"")) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  filter(month(Datum) >7) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  mutate(TS=ifelse(Beobachtung=="Fang",1,
                   ifelse((Beobachtung=="Senderlokalisation" & Ort=="UG"),1,
                          ifelse((Beobachtung=="Senderlokalisation" & Ort!="UG"),3,
                                 ifelse((Beobachtung=="Totfund" & Ort=="UG"),2,
                                        ifelse((Beobachtung=="Totfund" & Ort!="UG"),4,99)))))) %>% ## there should be no 99filter(TS==99)
  group_by(id, week) %>%
  summarise(N=max(TS)) %>%
  spread(key=week, value=N, fill=NA) %>%
  mutate(wk54=NA)  ## create blank column for records past the migration season

woco_ann_ch_obs<-woco %>%
  mutate(year=year(Datum), week=paste('wk',week(Datum),"")) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  filter(month(Datum) >7) %>%
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

woco.obs.matrix<-woco_ann_ch_obs %>% arrange(id)
woco.state.matrix<-woco_ann_ch_true %>% arrange(id)
dim(woco.obs.matrix)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### FILL MATRICES WITH STATE INFORMATION ###
for(i in 1:nrow(woco.obs.matrix)){
  
  ### extract ind and year
  yr<-separate_wider_delim(woco_ann_ch_obs[i,],cols="id",delim="_", names=c('ring','year'))$year
  rn<-separate_wider_delim(woco_ann_ch_obs[i,],cols="id",delim="_", names=c('ring','year'))$ring
  
  ### extract start and end dates for selected bird
  daterange<-woco %>%
    filter(Ring_num==rn) %>%
    filter(Beobachtung!="Senderfund") %>%
    filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
    summarise(first=min(Datum),last=max(Datum))
  
  ### ASSIGN OBSERVED STATES
  startcol<-min(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]!=6))
  if(startcol>1){
    if(year(daterange$first)<yr) {woco.obs.matrix[i,2:(startcol)]<-6} else {
      woco.obs.matrix[i,2:(startcol)]<-NA}
  }

  stopcol<-max(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]<6))
  if((stopcol+1)<dim(woco.obs.matrix)[2]){
    if(year(daterange$last)>yr) {woco.obs.matrix[i,dim(woco.obs.matrix)[2]]<-3} else {  ## birds that survived until next year are labelled to have been recorded outside study area
      woco.obs.matrix[i,(stopcol+2):dim(woco.obs.matrix)[2]]<-6}
  }

  
  ### ENSURE THAT DEAD BIRDS STAY DEAD
  if(TRUE %in% (c(4,5) %in% woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]])){
    deadcol<-min(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]] %in% c(4,5)))
    woco.obs.matrix[i,(deadcol+2):dim(woco.obs.matrix)[2]]<-woco.obs.matrix[i,(deadcol+1)]  ## assign the same state as last observed for rest of time series
  }

  ## ASSIGN INITIAL TRUE STATE (to initialise z-matrix of model)
  ### this needs careful manipulation to appropriately configure intermediate 0s
  startcol<-min(which(!is.na(woco.state.matrix[i,2:dim(woco.state.matrix)[2]])))
  stopcol<-max(which(!is.na(woco.state.matrix[i,2:dim(woco.state.matrix)[2]])))  
  startstate<-woco.state.matrix[i,startcol+1]
  endstate<-woco.state.matrix[i,stopcol+1]
  if(startcol>1){
    if(year(daterange$first)<yr) {woco.state.matrix[i,2:startcol]<-startstate} ## anything before first obs gets same state as first obs
  }
  if((stopcol+1)<dim(woco.state.matrix)[2]){
    woco.state.matrix[i,(stopcol+2):dim(woco.state.matrix)[2]]<-endstate ## anything after last obs gets same state as last obs
  }
  if(year(daterange$last)>yr) {
    woco.state.matrix[i,dim(woco.state.matrix)[2]]<-3  ## last column is always alive outside study area if bird also recorded next year
  }
}



#### Convert to numeric matrices that NIMBLE can loop over
y.telemetry<-as.matrix(woco.obs.matrix[,2:(dim(woco.state.matrix)[2])])
z.telemetry<-as.matrix(woco.state.matrix[,2:(dim(woco.state.matrix)[2])])


#### PREPARE A MATRIX OF WEEKS
nyears<-dim(woco_ch)[2]-1
nweeks<-dim(y.telemetry)[2]
nind<-dim(y.telemetry)[1]
week<-seq(1:nweeks)
year<-as.numeric(as.factor(separate_wider_delim(woco_ann_ch_obs,cols="id",delim="_", names=c('ring','year'))$year))
tag<-
effort<-

#### create vector of first marking and of last alive record
get.first.telemetry<-function(x)min(which(!is.na(x)))
get.last.telemetry<-function(x)max(which(!is.na(x) & x<6))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)
l.telemetry<-apply(y.telemetry,1,get.last.telemetry)


############# SAVE AND RE-LOAD PREPARED DATA ----------
save.image("data/woco_mig_input.RData")
load("data/woco_mig_input.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP NIMBLE CODE FOR SURVIVAL MODEL
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
  # 1 Bird (re)captured alive
  # 2 Bird recorded via transmitter in study area
  # 3 Bird recorded via transmitter outside study area
  # 4 Dead bird recovered  in study area
  # 5 Dead bird recovered  outside study area
  # 6 No signal (=not observed)
  
  # -------------------------------------------------
  
  # Priors and constraints
  
  
  #### WEEKLY SURVIVAL, MIGRATION and DETECTION PROBABILITIES
  for (i in 1:nind){
    for (t in f[i]:(n.occasions)){
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
    mean.phi[s] ~ runif(0.9,1)   # fairly uninformative prior for weekly survival probabilities
    #lp.mean[s] <- log(mean.phi[s]/(1 - mean.phi[s]))    # logit transformed survival intercept
    
    #### DEAD RECOVERY PROBABILITY VARIES BY YEAR
      mean.p.dead.in[s] ~ dunif(0,1)
      mean.p.dead.out[s] ~ dunif(0,0.3)
  }
  
  #### BASELINE FOR MIGRATION PROBABILITY (varies by year)
  mean.mig ~ runif(0,1)   # fairly uninformative prior for weekly survival probabilities
  lm.mean <- log(mean.mig/(1 - mean.mig))    # logit transformed migration intercept
  
  #### BASELINE FOR OBSERVATION PROBABILITY (varies by year)
  mean.p.in ~ runif(0,1)   # fairly uninformative prior for weekly detection probabilities
  lpin.mean <- log(mean.p.in/(1 - mean.p.in))    # logit transformed detection intercept
  
  mean.p.out ~ runif(0,0.3)   # fairly uninformative prior for weekly detection probabilities
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
    
    for (t in f[i]:(n.occasions-1)){
      
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      
      ps[1,i,t,1]<-1    ## dead birds stay dead
      ps[1,i,t,2]<-0
      ps[1,i,t,3]<-0
      ps[1,i,t,4]<-0
      
      ps[2,i,t,1]<-(1-phi[i,t])
      ps[2,i,t,2]<-phi[i,t] * (1-tag.fail[i,t])
      ps[2,i,t,3]<-phi[i,t] * tag.fail[i,t]
      ps[2,i,t,4]<-phi[i,t] * tag.fail[i,t]
      
      ps[3,i,t,1]<-(1-phi[i,t])
      ps[3,i,t,2]<-0
      ps[3,i,t,3]<-phi[i,t]
      ps[3,i,t,4]<-phi[i,t]
      
      ps[4,i,t,1]<-(1-phi[i,t])
      ps[4,i,t,2]<-0
      ps[4,i,t,3]<-phi[i,t]
      ps[4,i,t,4]<-phi[i,t]
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      
      po[1,i,t,1]<-0
      po[1,i,t,2]<-p.obs[i,t] * (1-tag.fail[i,t]) * (1-p.found.dead[i,t])
      po[1,i,t,3]<-0
      po[1,i,t,4]<-p.found.dead[i,t]
      po[1,i,t,5]<-(1-p.obs[i,t]) * tag.fail[i,t] * (1-p.found.dead[i,t])
      po[1,i,t,6]<-(1-p.obs[i,t]) * tag.fail[i,t] * (1-p.found.dead[i,t])
      
      po[2,i,t,1]<-(1-tag.fail[i,t])
      po[2,i,t,2]<-0
      po[2,i,t,3]<-0
      po[2,i,t,4]<-0
      po[2,i,t,5]<-tag.fail[i,t]
      po[2,i,t,6]<-tag.fail[i,t]
      
      po[3,i,t,1]<-0
      po[3,i,t,2]<-0
      po[3,i,t,3]<-p.obs[i,t]
      po[3,i,t,4]<-0
      po[3,i,t,5]<-(1-p.obs[i,t])
      po[3,i,t,6]<-(1-p.obs[i,t])
      
      po[4,i,t,1]<-0
      po[4,i,t,2]<-0
      po[4,i,t,3]<-p.obs[i,t]
      po[4,i,t,4]<-0
      po[4,i,t,5]<-(1-p.obs[i,t])
      po[4,i,t,6]<-(1-p.obs[i,t])
      
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




