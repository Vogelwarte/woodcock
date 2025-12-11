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
require(png)
library(grid)
library(gtable)
require(jpeg)
library(magick)
library(sf)
library(terra)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)
library(tidyterra)
library(cowplot)
filter<-dplyr::filter
select<-dplyr::select


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock"), silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock"), silent=T)

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


## 1.1 TABLE S2 ---------------

TableS2<-woco %>%
  dplyr::filter(Beobachtung=="Fang") %>%
  group_by(Ring) %>%
  summarise(first=min(Datum), age=min(age,na.rm=T)) %>%
  mutate(Year=year(first)) %>%
  mutate(age=ifelse(is.na(age),7,age)) %>%
  mutate(age=ifelse(age<4,"juvenile","adult")) %>%
  ungroup() %>%
  group_by(Year,age) %>%
  summarise(n=length(unique(Ring))) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from="age",values_from="n", values_fill=0) %>%
  adorn_totals()

fwrite(TableS2,"manuscript/Table_S2.csv")







# 2. FIGURE 1 ---------------

woco_mig<-readRDS("output/woco_mig_depart_simulation.rds")

FIGURE1<-woco_mig %>% 
  group_by(week) %>%
  summarise(mig=quantile(prop_mig,0.5),mig.lcl=quantile(prop_mig,0.025),mig.ucl=quantile(prop_mig,0.975)) %>%
  mutate(Date=lubridate::ymd("2024-07-26") + lubridate::weeks(week - 1)) %>% #print(n=35)
  
  ggplot()+
  geom_ribbon(aes(x=Date, ymin=mig.lcl, ymax=mig.ucl), alpha=0.2, fill="firebrick") +   ##
  geom_line(aes(x=Date, y=mig),linewidth=1, col="firebrick")+     ##

  ### add vertical lines to specify key dates OF OLD HUNTING TIMES
  geom_vline(aes(xintercept=min(Date[mig>0.95])), linetype="dashed", col="forestgreen", linewidth=1.5) +
  geom_segment(x=lubridate::ymd("2024-09-15"),y=0,yend=0.6, linetype="dashed", col="grey27", linewidth=2) +
  geom_segment(x=lubridate::ymd("2024-10-01"),y=0,yend=0.6, linetype="dashed", col="grey27", linewidth=2) +
  geom_segment(x=lubridate::ymd("2024-10-15"),y=0,yend=0.6, linetype="dashed", col="grey27", linewidth=2) +
  geom_segment(x=lubridate::ymd("2024-12-15"),y=0,yend=0.6, col="grey27", linewidth=1) +
  # geom_text(x=lubridate::ymd("2024-09-15"),y=0.65,label = "JU\nNE", size=6,col="grey36", vjust = 'bottom')+
  # geom_text(x=lubridate::ymd("2024-10-01"),y=0.65,label = "BE\nVS", size=6,col="grey27", vjust = 'bottom')+
  # geom_text(x=lubridate::ymd("2024-10-15"),y=0.65,label = "FR\nTI\nVD", size=6,col="grey18", vjust = 'bottom')+
  # 
  
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




## 2.1. CALCULATE ANNUAL SURVIVAL ---------
out<-fread("output/woco_telemetry_seasonal_surv_parm.csv")
out %>% filter(startsWith(parameter,"mean.phi")) %>%
  mutate(ann.surv=mean^52,lcl.ann.surv=lcl^52,ucl.ann.surv=ucl^52) %>%
  select(median, lcl, ucl, ann.surv,lcl.ann.surv,ucl.ann.surv) %>%
  slice_min(median) 
out %>% filter(startsWith(parameter,"mean.phi")) %>%
  mutate(ann.surv=mean^52,lcl.ann.surv=lcl^52,ucl.ann.surv=ucl^52) %>%
  select(median, lcl, ucl, ann.surv,lcl.ann.surv,ucl.ann.surv) %>%
  slice_max(median) 



# 3. FIGURE 2 ---------------

mean.p.nonlocal.migprior<- fread("output/WOCO_nonlocal_probs_mig_prior.csv")
mean.p.nonlocal<- fread("output/WOCO_nonlocal_probs_comb_prior.csv")
mean.p.nonlocal.null<- fread("output/WOCO_nonlocal_probs_no_prior.csv")



# mean.p.nonlocal.migprior<- fread("output/WOCO_nonlocal_probs_mig_prior.csv")
# mean.p.nonlocal<- fread("output/WOCO_nonlocal_probs_comb_prior.csv")
# mean.p.nonlocal.null<- fread("output/WOCO_nonlocal_probs_no_prior.csv")



FIGURE2<- bind_rows(mean.p.nonlocal,mean.p.nonlocal.migprior) %>%
  group_by(age,ctn,ind, prior) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(age,ctn, prior) %>%
  summarise(for.med=median(p.nonlocal.mean),for.ucl=quantile(p.nonlocal.mean,0.025), for.lcl=quantile(p.nonlocal.mean,0.975)) %>%
  bind_rows(mean.p.nonlocal.null) %>%
  mutate(prior=factor(prior, levels=c("only migration","combined abundance and migration","uninformative prior"))) %>%
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  #mutate(Kanton=levels(as.factor(woco.unk.sf$KANTON))[ctn]) %>%
  
  ggplot(aes(x=ctn, y=for.med))+
  geom_point(aes(col=Age, shape=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=for.lcl, ymax=for.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  facet_wrap(~prior, ncol = 1) +
  
  # annotation_custom(grob=gunicon, xmin=0.5, xmax=1.5, ymin=0.05, ymax=0.18) +
  # annotation_custom(wocoicon, xmin=0.5, xmax=2.9, ymin=0.10, ymax=0.35) +
  
  ## format axis ticks
  labs(y="Proportion of shot woodcocks of non-local origin",x="Swiss Canton",col="", shape="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  
  # viridis discrete color scale (cividis is very color-blind friendly)
  scale_color_viridis_d(option = "cividis", end = 0.9) +
  # complementary shapes for Age (helps in grayscale/print)
  scale_shape_manual(values = c("Adult" = 16, "Juvenile" = 17)) + # 16 = solid circle, 17 = solid triangle
  
  
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
        legend.position.inside=c(0.87,0.93),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE2

ggsave(plot=FIGURE2,
       filename="manuscript/Figure_2.jpg", 
       device="jpg",width=9, height=12, dpi=600)





### extrapolating annual harvest total of Swiss population

npairs<-c(1000,4000) ## from Knaus 2018
sexratio<-0.5
productivity<-1.6 ## from Kramer et al. 2019: https://www.sciencedirect.com/science/article/pii/S0006320718314149
shot<-1820
tot.birds<-(npairs/sexratio)+npairs*productivity

bind_rows(mean.p.nonlocal,mean.p.nonlocal.migprior) %>%
  group_by(age,ctn,ind, prior) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(prior) %>%
  summarise(for.lcl=quantile(p.nonlocal.mean,0.025), for.ucl=quantile(p.nonlocal.mean,0.975)) %>%
  mutate(min=(1-for.ucl)*shot/tot.birds[2], max=(1-for.lcl)*shot/tot.birds[1])






# 4. FIGURE S1 ---------------
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







# 5. FIGURE S2 ----------------------------------
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
  # viridis discrete color scale (cividis is very color-blind friendly)
  scale_fill_viridis_d(option = "cividis", end = 0.9) +
  
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
length(UNK_WC$dH)/9681
table(UNK_WC$AGE)
summary(ORIG_WC$dH)
summary(UNK_WC$dH)





# 6. FIGURE S3 ----------------------------------

ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
  dplyr::filter(!is.na(dH)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE)

## 6.1. add data of known origin provided by andrew Hoodless
ORIG_WC<-fread("data/WOCO_known_origin_feathers_Hoodless.csv") %>%
  mutate(ADULTE=ifelse(Age=="Adult",1,0),PROVENANCE_voigt=as.character(LocationCode),
         AGE=ifelse(Age=="Adult", "Adulte","Juvénile")) %>%
  mutate(dH=10.774 + (0.852*DHF)) %>%   ## new correction from Soto et al 2017 inserted after discussion with David Soto on 20 Oct 2025
  rename(ID=RefID, KANTON=Region,
         DATE=Date,LATITUDE=Latitude,LONGITUDE=Longitude) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
  bind_rows(ORIG_WC)
dim(ORIG_WC)

## 6.2. create spatial feature for plotting ------
woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)



## 6.3. load spatial background data and plot ------
bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)

WOCO.isoscape17<-readRDS("./data/isoscape17.rds") %>%
  mutate(mean=ifelse(mean==0,NA,mean)) %>%
  terra::crop(.,bbox)
plot(WOCO.isoscape17$mean)

FIGURES3<-ggplot() + 
  tidyterra::geom_spatraster(data=WOCO.isoscape17, aes(fill=mean))+
  geom_sf(data=woco.sf,color="red") +
  geom_sf(data=EUR, colour="grey12", fill=NA) +
  #ggtitle(unique(woco.unk.sf$KANTON)[ct]) +
  scale_fill_gradient(low = 'goldenrod', high = 'darkblue', na.value="white") +
  labs(fill = expression("Annual mean precipitation " * delta^2 * H)) +
  
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12, color="black"),
        axis.title=element_blank(),
        legend.position="inside",
        legend.position.inside=c(0.5,0.05),
        legend.direction = "horizontal",
        legend.background=element_blank(),
        plot.margin= unit(rep(.5, 4), "lines"),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURES3

ggsave(FIGURES3, filename="manuscript/FIGURE_S3.jpg",device="jpg",width=9, height=9)









# 7. FIGURE S4 ----------------------------------

#### completely revised figure on 26 Nov 2025
# prepared in WOCO_rainfall_isoscape_creation.r
load("data/woco.input.data.annual.RData")
woco.sf<-bind_rows(woco.unk.sf,woco.sf) %>%
  dplyr::filter(KANTON %in% unique(woco.unk.sf$KANTON))

bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)

woco.countries.orig <- EUR %>%
  dplyr::filter(admin %in% c("Ukraine","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
                             "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))

WOCO.isoscape02<-readRDS("./data/isoscape02.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape03<-readRDS("./data/isoscape03.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape05<-readRDS("./data/isoscape05.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape07<-readRDS("./data/isoscape07.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape08<-readRDS("./data/isoscape08.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape09<-readRDS("./data/isoscape09.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape10<-readRDS("./data/isoscape10.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape13<-readRDS("./data/isoscape13.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape14<-readRDS("./data/isoscape14.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape15<-readRDS("./data/isoscape15.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape16<-readRDS("./data/isoscape16.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape17<-readRDS("./data/isoscape17.rds") %>%
  terra::mask(woco.countries.orig)
WOCO.isoscape18<-readRDS("./data/isoscape18.rds") %>%
  terra::mask(woco.countries.orig)


isoscapes<-list(WOCO.isoscape02,WOCO.isoscape03,WOCO.isoscape05,WOCO.isoscape07,WOCO.isoscape08,WOCO.isoscape09,WOCO.isoscape10,
                WOCO.isoscape13,WOCO.isoscape14,WOCO.isoscape15,WOCO.isoscape16,WOCO.isoscape17,WOCO.isoscape18)




refscape<-WOCO.isoscape13  ## make sure that all isoscapes have the same extent by setting a reference isoscape
plot_list <- list()
for (ct in 1:length(unique(woco.unk.sf$KANTON))){
  
  ## get canton-wise distribution
  woco.cnt <- woco.sf %>% dplyr::filter(KANTON==unique(woco.unk.sf$KANTON)[ct]) 
  cnt.iso <-  woco.cnt %>% st_drop_geometry() %>%
    dplyr::select(ID,KANTON,d2h_MA,d2h_se_MA) 
  
  ## EXTRACT ISOSCAPE CELLS FALLING INTO THIS RANGE
  # Your observations (numeric vector)
  obs <- cnt.iso$d2h_MA
  probscape<-list()
  for (l in 1:length(isoscapes)){
    
    # Select the layers that hold the mean and (response) variance
    iso <- isoscapes[[l]][[c("mean", "mean_predVar")]]
    iso <-  terra::resample(iso,refscape, method="near")
    
    # mean probability per cell across all years
    probscape[[l]] <- terra::app(
      iso, fun = function(v, obs) {
        mu  <- v[1]
        var <- v[2]
        if (is.na(mu) || is.na(var) || var <= 0) return(NA_real_)
        sigma <- sqrt(var)
        mean(pnorm(obs, mean = mu, sd = sigma, lower.tail = FALSE))
      },
      obs = obs
    )
    
    names(probscape[[l]])<- "loglik"
    
  }
  
  CNT.isoscape <-  terra::app(sds(probscape),mean) %>%
    terra::crop(EUR)
  
  ## create plot
  plot_list[[ct]] <- ggplot(EUR) +
    geom_sf() +
    tidyterra::geom_spatraster(data=CNT.isoscape, aes(fill=loglik))+
    geom_sf(data=woco.cnt,color="red") +
    geom_sf(data=EUR, colour="grey12", fill=NA) +
    #ggtitle(unique(woco.unk.sf$KANTON)[ct]) +
    
    annotation_custom(wocoicon, xmin=-10, xmax=5, ymin=65, ymax=78) +
    annotate("text", x = 35, y = 75, label = unique(woco.unk.sf$KANTON)[ct], size=8, colour="darkolivegreen") +
    #scale_fill_gradient(limits=c(0,1), low = 'white', high = 'darkred', na.value="white") +
    #scale_fill_gradient(limits=c(0,1), low = 'white', high = 'darkred', na.value="white") +
    scale_fill_viridis_b(begin=0,end=1, alpha=0.8, option="E", direction=-1, na.value="white", n.breaks=7) +
    
    labs(fill = expression("Probability of indistinguishable rainfall " * delta^2 * H)) +
    
    
    # Important: remove expansion and keep tight
    coord_sf(expand = FALSE) +
    # Streamlined theme to reduce space
    theme(
      panel.background = element_blank(),   # remove box fill
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_blank(),
      legend.position = "none",             # we'll collect legend outside
      plot.margin = unit(rep(1, 4), "pt"),  # very tight outer margins
      strip.text = element_text(size = 18, color = "black"),
      strip.background = element_rect(fill = "white", colour = "black")
    )
  
  
} #ct




#FIG_S4<-grid.arrange(grobs=plot_list,ncol=2)



# Extract the legend from one plot
legend <- get_legend(plot_list[[1]] + theme(legend.position = "bottom"))

# # Remove legends from all plots
# plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

# Make margins small and panel spacing minimal for all plots
plots_tight <- lapply(plot_list, function(p) {
  p +
    theme(
      legend.position = "none",
      plot.margin   = margin(0, 0, 2, 2),         # tight outer margins (in pts)
      panel.spacing = unit(0, "lines")            # no internal panel spacing
    )
})



# Helper themes to quickly drop axes
drop_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

drop_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# Apply axis suppression by position in the 2-column grid
plots_axes <- lapply(seq_along(plots_tight), function(i) {
  p <- plots_tight[[i]]
  
  # Determine row (1..3) and column (1..2)
  row <- ceiling(i / 2)
  col <- ifelse(i %% 2 == 0, 2, 1)
  
  # Drop x-axis for rows 1 and 2 (keep only bottom row)
  if (row < 3) p <- p + drop_x
  
  # Drop y-axis for right column (keep only left column)
  if (col == 2) p <- p + drop_y
  
  p
})


# Combine plots and legend


# Combine 6 panels into a 2x3 grid
grid_panels <- plot_grid(
  plotlist = plots_axes,
  ncol = 2,
  align = "hv",       # helps align plotting areas
  labels = NULL
)


FIG_S4 <- plot_grid(
  grid_panels,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)  # Adjust legend height
)


# FIG_S4 <- plot_grid(plotlist = plots_no_legend, ncol = 2, align = "none",labels = NULL,label_size = 0 )
# FIG_S4 <- plot_grid(FIG_S4, legend, ncol = 1, rel_heights = c(1, 0.1),rel_widths = c(1,0.5))
FIG_S4



ggsave(filename="manuscript/Figure_S4.jpg", plot=FIG_S4,
       device="jpg",width=10, height=14)











# 8. FIGURE S5 -----------------------
UNK_WC<-woco %>% dplyr::filter(ORIGINE=="unknown") %>%
  dplyr::filter(JAGDT==1) %>%
  dplyr::filter(!is.na(DATE)) %>%
  dplyr::filter(!is.na(dH)) %>%
  mutate(Date=lubridate::parse_date_time(x=DATE,orders="dmy", tz = "UTC", drop=T)) %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,Date,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
  rename(DATE=Date)
UNK_WC$ID<-str_replace_all(UNK_WC$ID, "[^[:alnum:]]", " ")

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

#colors <- c("All birds" = "darkolivegreen", "Local birds" = "firebrick", "Shot birds" = "gray23")
colors <- c("All birds" = "#0C7BDC", "Local birds" = "#FFC20A", "Shot birds" = "gray23")

FIG_s5<-ggplot()+
  geom_line(data=woco_mig, aes(x=Date, y=mig, color="Local birds"),linewidth=2) +
  geom_line(data=woco_abundance, aes(x=Date, y=abund, color="All birds"),linewidth=2) +
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
FIG_s5

ggsave(plot=FIG_s5,
       filename="manuscript/Figure_S5.jpg",
       device="jpg",width=12, height=8)








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9. REPORT QUANTITIES IN MANUSCRIPT -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FILTER DATA by removing individuals with no record in Aug - Dec

removals<-woco %>%
  mutate(targetObs=ifelse(month(Datum)>7,1,0)) %>%
  group_by(Ring_num) %>%
  summarise(obs=sum(targetObs)) %>%
  filter(obs==0)

woco<-woco %>% filter(!(Ring_num %in% removals$Ring_num))

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







# 10. ABANDONED CODE FOR ALTERNATIVE FIGURES ------------------------------------------

mean.p.nonlocal.null<-fread("output/WOCO_nonlocal_probs_no_prior_featherscape.csv")
mean.p.nonlocal.migprior<-fwrite("output/WOCO_nonlocal_probs_mig_prior_featherscape.csv")
mean.p.nonlocal<-fwrite("output/WOCO_nonlocal_probs_featherscape.csv")


FIGURES2<- bind_rows(mean.p.nonlocal,mean.p.nonlocal.migprior) %>%
  group_by(age,ctn,ind, prior) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(age,ctn, prior) %>%
  summarise(for.med=median(p.nonlocal.mean),for.ucl=quantile(p.nonlocal.mean,0.025), for.lcl=quantile(p.nonlocal.mean,0.975)) %>%
  bind_rows(mean.p.nonlocal.null) %>%
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
        legend.position.inside=c(0.90,0.6),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURES2



## ABANDONED CODE FOR PREVIOUS FIG S3 ISOTOPE ORIGIN

# woco<-fread("data/WOCO_isotopes.csv")
# 
# woco$dH<-10.774 + (0.852*woco$dH_scaled)
# woco %>% dplyr::filter(is.na(dH))
# ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
#   dplyr::filter(!is.na(dH)) %>%
#   dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE)
# UNK_WC<-woco %>% dplyr::filter(ORIGINE=="UNBEKANNT") %>%
#   dplyr::filter(JAGDT==1) %>%
#   dplyr::filter(!is.na(DATE)) %>%
#   dplyr::filter(!is.na(dH)) %>%
#   mutate(Date=lubridate::parse_date_time(x=DATE,orders="dmy", tz = "UTC", drop=T)) %>%
#   dplyr::select(ID,ADULTE,AGE,KANTON,Date,dH,PROVENANCE_voigt,LATITUDE,LONGITUDE) %>%
#   rename(DATE=Date)
# UNK_WC$ID<-str_replace_all(UNK_WC$ID, "[^[:alnum:]]", " ")
# rain_orig_wc<-ORIG_WC %>% dplyr::select(-PROVENANCE_voigt,-ADULTE) %>%
#   mutate(sd=parse_date_time(DATE, orders="dmy")) %>%
#   mutate(feather_growth_date=if_else(AGE=="Adulte",
#                                      if_else(month(sd)<7,
#                                              ymd(paste(year(sd)-1,"-06-15")),
#                                              ymd(paste(year(sd),"-06-15"))),
#                                      if_else(month(sd)<9,
#                                              ymd(paste(year(sd)-1,"-05-15")),
#                                              ymd(paste(year(sd),"-08-15"))))) %>%
#   mutate(feather_sampling_date=as.Date(sd)) %>%
#   dplyr::select(-DATE,-sd,-AGE)
# 
# bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
# SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
#   st_transform(st_crs('EPSG:4326')) # Project to WGS84
# EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
#   st_crop(bbox) %>%
#   st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
#   dplyr::select(admin,name,adm0_a3,geometry)
# woco.sf <- ORIG_WC %>%
#   st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
# 
# 
# isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
#   crop(extent(SUI))
# woco.countries <- EUR %>%
#   dplyr::filter(admin %in% c("Ukraine","Switzerland","Sweden","Slovakia","Poland","Norway","Netherlands","Russia","Moldova","Luxembourg","Lithuania","Liechtenstein","Latvia",
#                              "Germany","Finland","Estonia","Denmark","Czechia","Belarus","Austria","Belgium"))
# forest.mat<-matrix(0, nrow=3, ncol=3)
# forest.mat[,1]<-c(11,40,110)  ## from values for conversion matrix
# forest.mat[,2]<-c(30,100,230)  ## to values for conversion matrix
# forest.mat[,3]<-c(0,1,0)  ## replacement values for conversion matrix
# globcover<-terra::rast("data/GLOBCOVER_L4_200901_200912_V2.3.tif") %>%
#   crop(woco.countries) %>%
#   terra::classify(rcl=forest.mat,include.lowest=T,right=NA) %>%
#   terra::project(.,crs(isoscape))
# WOCO.isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
#   crop(woco.countries)
# globcover <- terra::resample(globcover,WOCO.isoscape, method="max")
# WOCO.isoscape <- WOCO.isoscape*globcover
# rain.d2H<-as.numeric(na.omit(terra::values(WOCO.isoscape)[,1]))
# rain.d2H<-rain.d2H[rain.d2H<0]  ## remove the non-forest values (>10,0000 grid cells are removed)
# woco.vect<-terra::vect(woco.sf)
# woco.sf$d2h_MA<-terra::extract(isoscape,woco.vect)$d2h_MA
# woco.sf$d2h_se_MA<-terra::extract(isoscape,woco.vect)$d2h_se_MA
# woco.unk.sf <- UNK_WC %>%
#   st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)
# woco.unk.vect<-terra::vect(woco.unk.sf)
# 
# woco.unk.sf$d2h_MA<-terra::extract(isoscape,woco.unk.vect)$d2h_MA
# woco.unk.sf$d2h_se_MA<-terra::extract(isoscape,woco.unk.vect)$d2h_se_MA
# woco.unk.sf %>% filter(is.na(AGE))
# 
# 
# 
# woco.unk.sf <- woco.unk.sf %>%
#   filter(!is.na(AGE)) %>%
#   filter(!is.na(dH)) %>%
#   #filter(KANTON !="VS") %>% ## remove Valais because only 6 birds from 1 age class, causes imbalance in data
#   filter(!is.na(d2h_MA))
# 
# 
# 
# plot_list <- list()
# for (ct in 1:length(unique(woco.unk.sf$KANTON))){
#   
#   ## get canton-wise distribution
#   woco.cnt <- woco.unk.sf %>% dplyr::filter(KANTON==unique(woco.unk.sf$KANTON)[ct]) 
#   cnt.iso <-  woco.cnt %>% st_drop_geometry() %>%
#     dplyr::select(ID,KANTON,d2h_MA,d2h_se_MA) %>%
#     rowwise() %>%
#     mutate(pot.orig.d2H=rnorm(1,d2h_MA,d2h_se_MA)) %>%
#     ungroup() %>%
#     group_by(KANTON) %>%
#     summarise(min=min(pot.orig.d2H), max=max(pot.orig.d2H))
#   
#   ## EXTRACT ISOSCAPE CELLS FALLING INTO THIS RANGE
#   CNT.isoscape <- ifel(WOCO.isoscape >= cnt.iso$min & WOCO.isoscape <= cnt.iso$max, 1, 0)
#   
#   ## create plot
#   plot_list[[ct]] <- ggplot(EUR) +
#     geom_sf() +
#     tidyterra::geom_spatraster(data=CNT.isoscape, aes(fill=d2h_MA))+
#     geom_sf(data=woco.cnt,color="red") +
#     geom_sf(data=EUR, colour="grey12", fill=NA) +
#     #ggtitle(unique(woco.unk.sf$KANTON)[ct]) +
#     
#     annotation_custom(wocoicon, xmin=-10, xmax=5, ymin=65, ymax=78) +
#     annotate("text", x = 35, y = 75, label = unique(woco.unk.sf$KANTON)[ct], size=8, colour="darkolivegreen") +
#     scale_fill_gradient(low = 'lightgray', high = 'lightgreen', na.value=NA) +
#     
#     ## beautification of the axes
#     theme(panel.background=element_rect(fill="white", colour="black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.text=element_text(size=12, color="black"),
#           axis.title=element_blank(),
#           legend.position="none",
#           plot.margin= unit(rep(.5, 4), "lines"),
#           strip.text=element_text(size=18, color="black"),
#           strip.background=element_rect(fill="white", colour="black"))
#   
# } #ct









