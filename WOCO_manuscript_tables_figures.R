###########################################################################
##### WOODCOCK AUTUMN MIGRATION MANUSCRIPT TABLES AND FIGURES--------  ############################
###########################################################################

## written by Steffen Oppel on 25 Sept 2025
## creating tables and figures for manuscript by loading in output from models and raw data


# Clear workspace ---------------------------------------------------------

rm(list=ls()) 
Sys.setenv(LANG = "en") ##change language of error messages to english

# Load libraries, functions etc -------------------------------------------
library(lubridate)
library(tidyverse)
library(readxl)
library(dplyr)
library(janitor)
filter<-dplyr::filter
select<-dplyr::select
require(png)
library(grid)
library(gtable)
require(jpeg)
library(magick)



## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")
imgGun<-readPNG("manuscript/rifleicon.png")
gunicon <- rasterGrob(imgGun, interpolate=TRUE)
imgWOCO<-image_read("manuscript/woodcock.jpg") %>% image_transparent("white", fuzz=5)
wocoicon <- rasterGrob(imgWOCO, interpolate=TRUE)


# Load data ---------------------------------------------------------------

woco <- read_excel("data/Daten_Waldschnepfe.xlsx", sheet="Data") %>%
  mutate(Beobachtung=ifelse(is.na(Beobachtung),"Fang",Beobachtung)) %>%  ## fill in missing values for type of transmitter
  mutate(Sendertyp=ifelse(is.na(Sendertyp),Markierung,Sendertyp)) %>%  ## fill in missing values for type of transmitter
  rename(age=`Alter (Euring)`)
woco


# 1. TABLE S1 ---------------
unique(woco$Beobachtung)
# weirdobs<-woco %>% filter(is.na(Beobachtung))
# woco %>% filter(Ring %in% weirdobs$Ring)


all_birds<-woco %>%
  dplyr::filter(Beobachtung=="Fang") %>%
  group_by(Ring) %>%
  summarise(first=min(Datum),total=length(Datum)) %>%
  mutate(season=if_else(month(first) %in% c(5,6,7,8),"breeding","migration")) %>%
  ungroup() %>%
  group_by(season) %>%
  summarise(N=length(unique(Ring))) %>%
  adorn_totals()

dead_birds<-woco %>%
  dplyr::filter(Beobachtung=="Totfund") %>%
  group_by(Ort, Todesursache) %>%
  summarise(N=length(unique(Ring))) %>%
  adorn_totals()

tagged_birds<-woco %>%
  dplyr::filter(Beobachtung=="Fang") %>%
  group_by(Ring, Sendertyp) %>%
  summarise(first=min(Datum),total=length(Datum))

argos_birds<-tagged_birds %>%
  dplyr::filter(Sendertyp=="ARGOS")
argos_birds %>%
  ungroup() %>%
  summarise(N=length(unique(Ring)))

vhf_birds<-tagged_birds %>%
  dplyr::filter(Sendertyp=="VHF") %>%
  dplyr::filter(!(Ring %in% argos_birds$Ring))
vhf_birds %>%
  ungroup() %>%
  summarise(N=length(unique(Ring)))

ringonly_birds<-tagged_birds %>%
  dplyr::filter(!(Ring %in% vhf_birds$Ring)) %>%
  dplyr::filter(!(Ring %in% argos_birds$Ring))
ringonly_birds %>%
  ungroup() %>%
  summarise(N=length(unique(Ring)))


# doubletagged_birds<-woco %>%
#   dplyr::filter(Beobachtung=="Fang") %>%
#   dplyr::filter((Ring %in% argos_birds$Ring)) %>%
#   group_by(Ring) %>%
#   summarise(N=length(unique(Sendertyp))) %>%
#   dplyr::filter(N>1)
# 
# woco %>%
#   dplyr::filter(Ring %in% doubletagged_birds$Ring) %>%
#   print(n=50)

singlecaps <- woco %>%
  group_by(Ring) %>%
  summarise(N=length(Datum)) %>%
  filter(N==1) %>%
  ungroup() %>%
  summarise(N=length(unique(Ring)))






# 2. FIGURE S1 ---------------
medlast<-woco %>% filter(month(Datum) %in% c(8,9,10,11)) %>%
  filter(Ort=="UG") %>%
  mutate(year=year(Datum)) %>%
  mutate(id=paste(Ring_num,year, sep="_")) %>%
  mutate(jday=yday(Datum)) %>%
  mutate(date=ymd("2023-12-31")+days(jday)) %>%
  group_by(id) %>%
  summarise(last=max(date)) %>%
  ungroup() %>%
  summarise(day=median(last))

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
  ### add vertical lines to specify median date of last observation
  geom_vline(aes(xintercept=medlast$day), linetype="dashed", col="firebrick", linewidth=1.5) +
  scale_x_date(name="last observation in study area",date_breaks="1 week", date_labels="%d-%b")+
  labs(y="proportion of woodcocks") +
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=12, color="black",angle=45,hjust = 1),
        axis.text.y=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
#ggsave("manuscript/Figure_S1.jpg", width=8, height=6)






# 3. FIGURE 1 ---------------

woco_mig<-readRDS("output/woco_mig_depart_simulation.rds")

FIGURE1<-woco_mig %>% 
  group_by(week) %>%
  summarise(mig=quantile(prop_mig,0.5),mig.lcl=quantile(prop_mig,0.025),mig.ucl=quantile(prop_mig,0.975)) %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>%
  
  ggplot()+
  geom_ribbon(aes(x=Date, ymin=mig.lcl, ymax=mig.ucl), alpha=0.2, fill="firebrick") +   ##
  geom_line(aes(x=Date, y=mig),linewidth=1, col="firebrick")+     ##

  ### add vertical lines to specify key dates OF OLD HUNTING TIMES
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
  scale_y_continuous(name="Cumulative proportion that have left study area", limits=c(0,1), breaks=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE1

ggsave(plot=FIGURE1,
       filename="manuscript/Figure_1.jpg", 
       device="jpg",width=11, height=8)








# 4. FIGURE 2 ---------------

mean.p.nonlocal.migprior<- fread("output/WOCO_nonlocal_probs_mig_prior.csv")
mean.p.nonlocal<- fread("output/WOCO_nonlocal_probs_comb_prior.csv")

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
  facet_wrap(~prior, ncol = 2) +
  
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
        legend.position.inside=c(0.90,0.85),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE2


ggsave(plot=FIGURE2,
       filename="output/woco_iso_time_origin_probability_estimates_comb_prior.jpg", 
       device="jpg",width=9, height=12)



# 5. FIGURE S4 -----------------------

woco_mig<-readRDS("output/woco_mig_depart_simulation.rds") %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>%
  mutate(prop.remain=1-prop_mig) %>% 
  group_by(week, Date) %>%
  summarise(mig=quantile(prop.remain,0.5),mig.lcl=quantile(prop.remain,0.025),mig.ucl=quantile(prop.remain,0.975))

woco_abundance<-fread("data/AnnualPhenology_abundance_index.csv") %>%
  mutate(Date=lubridate::ymd("2023-12-27") + lubridate::days(Pentade*5)) %>%
  mutate(abund=SOPM/max(SOPM)) %>%
  dplyr::filter(Date>=min(woco_mig$Date))

shot_dates<-hist(lubridate::yday(UNK_WC$DATE), breaks=seq(250,365,7),plot=F)
woco_shot<-tibble(yday=shot_dates$mids, N=shot_dates$counts) %>%
  mutate(Date=parse_date_time(as.integer(yday), orders="j")) %>%
  mutate(abund=N/max(N)) %>%
  mutate(Date=as.Date(Date-years(1)))

colors <- c("All birds" = "darkolivegreen", "Local birds" = "firebrick", "Shot birds" = "gray23")

FIG_s4<-ggplot()+
  geom_line(data=woco_mig, aes(x=Date, y=mig, color="Local birds"),linewidth=1) +
  geom_line(data=woco_abundance, aes(x=Date, y=abund, color="All birds"),linewidth=1) +
  geom_col(data=woco_shot, aes(x=Date, y=abund, color="Shot birds"),width = 6, alpha=0.5) +
  labs(color = "Origin") +
  scale_color_manual(values = colors) +
  
  ## format axis ticks
  scale_x_date(name="Date", date_labels = "%d %b") +
  scale_y_continuous(name="Relative abundance of woodcocks", limits=c(0,1), breaks=seq(0,1,0.2)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title=element_text(size=16, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.position="inside",
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.position.inside=c(0.1,0.5),
        strip.background=element_rect(fill="white", colour="black"))
FIG_s4

ggsave(plot=FIG_s4,
       filename="manuscript/Figure_S4.jpg",
       device="jpg",width=11, height=8)






# 7. FIGURE S2 ----------------------------------
load("data/woco.input.data.RData")
woco$ORIGINE<-ifelse(woco$ORIGINE=="SCHWEIZ","Switzerland","unknown")

binwidth <- 5
breaks <- seq(min(woco$dH,na.rm=T), max(woco$dH,na.rm=T) + binwidth, by = binwidth)
FIGURES2<-woco %>%
  mutate(bin = cut(dH, breaks = breaks, right = FALSE)) %>%
  group_by(bin) %>%
  mutate(bin_median = median(dH)) %>%
  ungroup() %>%
  group_by(ORIGINE, bin, bin_median) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ORIGINE) %>%
  mutate(percentage = count / sum(count)) %>%
  
  ggplot(aes(x = round(bin_median, 1), y = percentage, fill = ORIGINE)) +
  geom_bar(stat = "identity", position = position_dodge(width=2.5,preserve = "single"), alpha = 0.5, width = binwidth*0.8) +

  #geom_histogram(alpha=0.5,position = position_dodge(width=1)) +
  #geom_histogram(aes(y = ..count.. / tapply(..count.., ..PANEL.., sum)[..PANEL..]),binwidth = 5,  alpha = 0.5,position = "dodge") +
  #geom_histogram(aes(y=stat(density)),binwidth=,alpha=0.5,position = position_dodge(width=2.5))+  
  
  ## format axis ticks
  labs(y="Proportion of woodcock feathers",
       x=expression(paste(delta^{2}, "H (\u2030)")),
       col="Origin", fill="Origin") +
  scale_y_continuous(labels = scales::percent_format(scale=100)) +
  
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
        legend.position.inside=c(0.13,0.75),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURES2
#ggsave(FIGURES2, filename="manuscript/FIGURE_S2.jpg",device="jpg",width=12, height=9)


## report numbers in manuscript
table(ORIG_WC$AGE)
length(UNK_WC$dH)/9260
table(UNK_WC$AGE)
summary(ORIG_WC$dH)
summary(UNK_WC$dH)














# 7. FIGURE S3 ----------------------------------
rm(list=ls())
load("data/woco.input.data.RData")
try(rm(isoscape, globcover), silent=T)

woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  #filter(KANTON !="VS") %>% ## remove Valais because only 6 birds from 1 age class, causes imbalance in data
  filter(!is.na(d2h_GS))



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

grid.arrange(grobs=plot_list,ncol=3)
ggsave(filename="output/woco_iso_origin_local_maps.jpg", 
       device="jpg",width=8, height=12)
ggsave(filename="output/woco_iso_origin_local_maps.jpg", 
       device="jpg",width=8, height=12)

























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
# REPORT QUANTITIES IN MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inspect data ---------------------------------------------------------------
table(woco$Ort)
table(woco$age)
hist(month(woco$Datum))
table(woco$Beobachtung)
table(woco$Markierung)
table(woco$Markierung, woco$Sendertyp)
table(woco$Censor)
length(unique(woco$Ring_num))
woco %>% filter(month(Datum) %in% c(8,9,10,11)) %>%
  mutate(jday=yday(Datum)) %>%
  ggplot(aes(x=jday)) + geom_histogram()

table(woco$Beobachtung, woco$Ort)


woco %>% group_by(Ring, Markierung, Sendertyp) %>%
  summarise(n_encounter=length(Datum)) %>%
  ungroup() %>%
  mutate(N=1) %>%
  group_by(Markierung, Sendertyp) %>%
  summarise(N_birds=sum(N))





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
save.image("data/woco_mig_input_no_argos.RData")
