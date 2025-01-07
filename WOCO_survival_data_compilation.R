###########################################################################
##### WOODCOCK DEPARTURE TIME ANALYSIS DATA PREPARATION--------  ############################
###########################################################################

## written by Steffen Oppel on 30 Oct 2024 to initiate analysis
## goal is to estimate WHEN birds depart so that hunters can be confident they don't shoot Swiss birds
## federal regulation permits shooting from 15 Sept to 16 December

## QUESTIONS TO MURDY:
# what does 'censor' mean? check with Michi if censoring necessary for migration
# what temporal resolution do we want (weekly ok)? yes
# extent of season (Sept - Nov), what is the intermittent dip in effort? monthly capture events
# do we have a metric of observation effort?
# what are age categories 1-8? EURING Code? Only needed if age affects departure or survival
# annual variation in survival, detection, migration?

## NEED TO DO:
# weed out non-informative birds
# overwrite short excursions outside of the study area

## UPDATED 2 JAN 2025 to include actual netting effort data (provided by Pierre Mollet)
## used data for those occasions where available, and extrapolated remaining occasions based on linear regression with capture numbers

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
  geom_histogram(aes(x=last,y=(after_stat(count))/tapply(after_stat(count),after_stat(PANEL),sum)[after_stat(PANEL)]),binwidth=7)+                               
  scale_x_date(name="last obs in study area",date_breaks="1 week", date_labels="%d-%b")+
  labs(y="proportion of woodcocks") +
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=12, color="black",angle=45,hjust = 1),
        axis.text.y=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
#ggsave("output/prop_woodcock_per_week_rawdat.jpg", width=7, height=6)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER DATA by removing individuals with no record in Aug - Dec
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
removals<-woco %>%
  mutate(targetObs=ifelse(month(Datum)>7,1,0)) %>%
  group_by(Ring_num) %>%
  summarise(obs=sum(targetObs)) %>%
  filter(obs==0)

woco<-woco %>% filter(!(Ring_num %in% removals$Ring_num))


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
  # 1 Bird (re)captured alive or recorded via transmitter in study area
  # 2 Bird recorded via transmitter outside study area
  # 3 Dead bird recovered  in study area
  # 4 Dead bird recovered  outside study area
  # 5 No signal (=not observed)

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
### we only want weeks from Aug - December
### but some birds are marked earlier in the year and are then shot - cannot have a dead state without a previous marking occasion
### summarised all encounters prior to August and entered it in first August week
### weekly encounter history ####
woco_ann_ch_true<-woco %>%
  mutate(year=year(Datum), week=week(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  mutate(week=ifelse(week<31,31,week)) %>%
  filter(week>30) %>%
  mutate(week=paste0('wk',week)) %>%
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
  mutate(year=year(Datum), week=week(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  mutate(week=ifelse(week<31,31,week)) %>%
  filter(week>30) %>%
  mutate(week=paste0('wk',week)) %>%
  #filter(month(Datum) >7) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  mutate(OS=ifelse(Beobachtung=="Fang",1,
                   ifelse((Beobachtung=="Senderlokalisation" & Ort=="UG"),1,
                          ifelse((Beobachtung=="Senderlokalisation" & Ort!="UG"),2,
                                 ifelse((Beobachtung=="Totfund" & Ort=="UG"),3,
                                        ifelse((Beobachtung=="Totfund" & Ort!="UG"),4,5)))))) %>%
  group_by(id, week) %>%
  summarise(N=max(OS)) %>%
  spread(key=week, value=N, fill=5) %>%
  mutate(wk54=5)  ## create blank column for records past the migration season


woco_tag_mat<- woco %>%
  mutate(year=year(Datum), week=week(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  select(id, Ring_num,week,Datum,Beobachtung,Markierung,Ort) %>%
  mutate(week=ifelse(week<31,31,week)) %>%
  filter(week>30) %>%
  mutate(week=paste0('wk',week)) %>%
  #filter(month(Datum) >7) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  mutate(tag=ifelse(Markierung=="Sender",1,0)) %>%
  group_by(id) %>%
  summarise(tag=max(tag))


### CREATE BLANK MATRICES TO HOLD INFORMATION ABOUT TRUE AND OBSERVED STATES ###

woco.obs.matrix<-woco_ann_ch_obs %>% arrange(id)
woco.state.matrix<-woco_ann_ch_true %>% arrange(id)
tag<-as.data.frame(woco_tag_mat %>% arrange(id) %>% select(tag))[,1]
woco.eff.matrix<-woco_ann_ch_obs %>% arrange(id)
dim(woco.obs.matrix)
length(tag)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE EFFORT MATRIX BASED ON ACTUAL NETTING EFFORT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### full matrix of records to be used to calibrate effort data from netting
effort_mat<- woco %>%
  mutate(year=year(Datum), week=week(Datum)) %>%
  mutate(week=ifelse(week<31,31,week)) %>%
  mutate(obs=1) %>%
  filter(week>30) %>%
  mutate(week=paste0('wk',week)) %>%
  filter(Beobachtung!="Senderfund") %>%
  filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
  group_by(year,week) %>%
  summarise(eff=sum(obs)) %>%
  spread(key=week,value=eff, fill=0) %>%
  mutate(wk54=2)  ## create blank column for records past the migration season

### summary of actual netting effort 
netting <- read_excel("data/Data_Aufwand_Nov2024.xlsx", sheet="Tabelle1") %>%
  mutate(effort=(`Länge (m)`*`Dauer (min)`)/60) %>%
  mutate(year=year(Datum), week=week(Datum)) %>%
  mutate(week=ifelse(week<31,31,week)) %>%
  filter(week>30) %>%
  mutate(week=paste0('wk',week)) %>%
  group_by(year,week) %>%
  summarise(tot_eff=sum(effort))
  
### calibrate relationship between netting effort and number of records
effort_mat %>% filter(year %in% c(2016,2017,2018)) %>%
  gather(key="week", value="exp_eff",-year) %>%
  left_join(netting, by=c("year","week")) %>%
  filter(!is.na(tot_eff)) %>%
  
  ggplot(aes(x=tot_eff, y=exp_eff)) +
  geom_point(aes(col=week),size=2.5) +
  geom_smooth(method="lm") +
  labs(x="Weekly net metre hours", y="Total N of woodcock captures in UG") +
  theme(legend.position="none")

## fit linear regression
calib_dat<-effort_mat %>% filter(year %in% c(2016,2017,2018)) %>%
  gather(key="week", value="exp_eff",-year) %>%
  left_join(netting, by=c("year","week")) %>%
  filter(!is.na(tot_eff))
eff_pred<-lm(tot_eff~0+exp_eff, data=calib_dat) ## forced through origin


## interpolate all occassions without actual netting effort data
effort<-effort_mat %>%
  gather(key="week", value="exp_eff",-year) %>%
  ungroup() %>%
  left_join(netting, by=c("year","week"))
effort$pred_eff<-predict(eff_pred, newdat=effort) ## predict effort based on linear model
effort_mat<-effort %>%
  mutate(eff=ifelse(is.na(tot_eff),pred_eff,tot_eff)) %>% ## use the predicted effort for occasions when no actual effort is available
  mutate(scale_eff=scale(eff)) %>%
  select(year,week,scale_eff) %>%
  spread(key=week,value=scale_eff)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### FILL MATRICES WITH STATE INFORMATION ###
for(i in 1:nrow(woco.obs.matrix)){
  
  ### extract ind and year
  yr<-separate_wider_delim(woco_ann_ch_obs[i,],cols="id",delim="_", names=c('ring','year'))$year
  rn<-separate_wider_delim(woco_ann_ch_obs[i,],cols="id",delim="_", names=c('ring','year'))$ring
  
  ### add effort to individual matrix (because nimble does not allow dynamic indexing)
  woco.eff.matrix[i,2:dim(woco.obs.matrix)[2]]<-effort_mat[effort_mat$year==yr,2:dim(woco.obs.matrix)[2]]
  woco.eff.matrix[i,dim(woco.obs.matrix)[2]]<-sum(effort_mat[effort_mat$year==yr,2:20])  ### effort for following year is the sum across many months
  
  ### extract start and end dates for selected bird
  daterange<-woco %>%
    filter(Ring_num==rn) %>%
    filter(Beobachtung!="Senderfund") %>%
    filter(!(Beobachtung=="Fang" & Ort!="UG")) %>%
    summarise(first=min(Datum),last=max(Datum))
  
  ### ASSIGN OBSERVED STATES
  startcol<-min(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]!=5))
  if(startcol>1){
    if(year(daterange$first)<yr) {woco.obs.matrix[i,2:(startcol)]<-5} else {
      woco.obs.matrix[i,2:(startcol)]<-NA}
  }

  stopcol<-max(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]<5))
  if((stopcol+1)<dim(woco.obs.matrix)[2]){
    if(year(daterange$last)>yr) {woco.obs.matrix[i,dim(woco.obs.matrix)[2]]<-2} else {  ## birds that survived until next year are labelled to have been recorded outside study area
      woco.obs.matrix[i,min((stopcol+2),dim(woco.obs.matrix)[2]):dim(woco.obs.matrix)[2]]<-5}
  }

  
  ### ENSURE THAT DEAD BIRDS STAY DEAD
  if(TRUE %in% (c(3,4) %in% woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]])){
    deadcol<-min(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]] %in% c(3,4)))
    woco.obs.matrix[i,min((deadcol+2),dim(woco.obs.matrix)[2]):dim(woco.obs.matrix)[2]]<-woco.obs.matrix[i,(deadcol+1)]  ## assign the same state as last observed for rest of time series
  }
  
  ### ENSURE THAT LOCAL BIRDS STAY LOCAL
  ## there are very few 'excursions' out of the study area, but we do not want to model those, so we eliminate any state other than 1 prior to or between states of 1
  startcolIN<-min(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]<5))
  stopcolIN<-max(which(woco.obs.matrix[i,2:dim(woco.obs.matrix)[2]]==1))
  if(stopcolIN>startcolIN){
    for(col in (startcolIN+1):(stopcolIN+1)) {
      if(!is.na(woco.obs.matrix[i,col])){woco.obs.matrix[i,col]<-ifelse(woco.obs.matrix[i,col]==5,5,1)}  ## assign everything to 1 except the not observed state
    }
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
    woco.state.matrix[i,min((stopcol+2),dim(woco.state.matrix)[2]):dim(woco.state.matrix)[2]]<-endstate ## anything after last obs gets same state as last obs
  }
  
  if(stopcol>(startcol+1)){ ## only needed if there is a chance of intermediate gaps between start and stopcol
    for(col in (startcol+1):stopcol) {
    woco.state.matrix[i,col]<-ifelse(is.na(woco.state.matrix[i,col]),woco.state.matrix[i,col-1],woco.state.matrix[i,col]) ## anything inbetween the start and stop column gets the previous state unless it is known
  }}
  
  if(year(daterange$last)>yr) {
    woco.state.matrix[i,dim(woco.state.matrix)[2]]<-3  ## last column is always alive outside study area if bird also recorded next year
  }
  if(stopcolIN>startcolIN){
    for(col in (startcolIN+1):(stopcolIN+1)) {
      woco.state.matrix[i,col]<-1  ## assign everything to 1 for true state
    }
  }
}



#### Convert to numeric matrices that NIMBLE can loop over
y.telemetry<-as.matrix(woco.obs.matrix[,2:(dim(woco.state.matrix)[2])])
z.telemetry<-as.matrix(woco.state.matrix[,2:(dim(woco.state.matrix)[2])])


#### REMOVE individuals with no information (i.e. those that left or were shot before August in a given year)
noninfobirds<-which(apply(y.telemetry, 1, function(x) length(unique(x)) == 1) == TRUE)
y.telemetry<-y.telemetry[-noninfobirds,]
z.telemetry<-z.telemetry[-noninfobirds,]
woco_ann_ch_obs<-woco_ann_ch_obs[-noninfobirds,]
woco.eff.matrix<-woco.eff.matrix[-noninfobirds,]
tag<-tag[-noninfobirds]

#### RETAIN ONLY individuals that were once seen alive in study area (all others have no value for estimating WHEN live birds leave study area)
UKbirds<-which(apply(y.telemetry, 1, function(x) 1 %in% unique(x)) == TRUE)
y.telemetry<-y.telemetry[UKbirds,]
z.telemetry<-z.telemetry[UKbirds,]
woco_ann_ch_obs<-woco_ann_ch_obs[UKbirds,]
woco.eff.matrix<-woco.eff.matrix[UKbirds,]
tag<-tag[UKbirds]

#### PREPARE A MATRIX OF WEEKS
nyears<-dim(woco_ch)[2]-1
nweeks<-dim(y.telemetry)[2]
nind<-dim(y.telemetry)[1]
week<-seq(1:nweeks)
year<-as.numeric(as.factor(separate_wider_delim(woco_ann_ch_obs,cols="id",delim="_", names=c('ring','year'))$year))
effort<-woco.eff.matrix[,-1]

#### create vector of first marking and of last alive record
get.first.telemetry<-function(x)min(which(!is.na(x)))
get.last.telemetry<-function(x)max(which(!is.na(x) & x<6))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)
l.telemetry<-apply(y.telemetry,1,get.last.telemetry)


############# SAVE  PREPARED DATA ----------
save.image("data/woco_mig_input.RData")
