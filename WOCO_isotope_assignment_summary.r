### SUMMARY OF ISOTOPE ASSIGNMENTS ----
## Table 14 in Bohnenstengel et al. Report
## digitized by Steffen Oppel on 8 May 2025

library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
library(readxl)

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")

woco_n<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table11_SampleSize")
woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")
woco_juv<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table15_IsotopeAssignments_JUV")
woco_ad<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table17_IsotopeAssignments_AD")
sum(woco$N)/9260


## total number shot and % of sample (reconstructed from report)
## impossible to know what came from "NorthernAlps" because that is not a canton
## Between 2013 and 2017, Swiss woodcock hunters provided 850 samples of woodcock feathers for isotope analysis (650 feathers from juveniles, 190 from adults and 10 from birds of unknown age)out of 9260 woodcocks taken in the hunt, i.e. 9.2 % of the total harvest for these years
## (Percentages of harvest according to canton: Fribourg: 21.5 %; Jura: 7.0 %; Neuchâtel: 16.7 %; Ticino: 7.1 %; Vaud: 14.9 %; Valais: 14.3 %).

woco_n<-woco_n %>%
	mutate(RegionShot=if_else(Origin %in% c("Tessin","Wallis"), "CentSouthAlps",
					if_else(Origin %in% c("Jura","Waadt","Doubs (F)","Neuchâtel"), "Jura",  # 
						if_else(Origin %in% c("Freiburg"), "Plateau","NorthernAlps",))))

woco_n %>% group_by(RegionShot) %>%
	summarise(N=sum(N_sample))


### the total should add up to 9260 but it doesn't!!
woco %>% group_by(RegionShot) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=c((0.071+0.143)/2,(0.07+0.167+0.149)/3,0.13,0.215)) %>% ## GUESSWORK ASSIGNMENT which canton represents plateau (Fribourg) and Jura (JU,NE,VD) - unknown what is "Northern Alps"??
  mutate(Nshot=as.integer(Ntot/prop)) %>%
  adorn_totals()





## summarise proportion of birds originating from outside Switzerland

ALL_OUT<-tibble()

### CREATE DIFFERENT SUBSETS AND OUTPUT METRICS

# ALL BIRDS ALL CANTONS

ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="all", origin="min") %>%
  bind_rows(ALL_OUT)
  


ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="all", origin="max") %>%
  bind_rows(ALL_OUT)


# ALL BIRDS N Cantons

ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="N_SUI", origin="min") %>%
  bind_rows(ALL_OUT)



ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="N_SUI", origin="max") %>%
  bind_rows(ALL_OUT)


# ALL BIRDS only JU, NE, VD

ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  filter(RegionShot!="NorthernAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="JU_PL", origin="min") %>%
  bind_rows(ALL_OUT)

ALL_OUT<-woco %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  filter(RegionShot!="NorthernAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="all", shot_in="JU_PL", origin="max") %>%
  bind_rows(ALL_OUT)






# JUV BIRDS ALL CANTONS

ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="all", origin="min") %>%
  bind_rows(ALL_OUT)



ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="all", origin="max") %>%
  bind_rows(ALL_OUT)


# JUV BIRDS N Cantons

ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="N_SUI", origin="min") %>%
  bind_rows(ALL_OUT)



ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="N_SUI", origin="max") %>%
  bind_rows(ALL_OUT)


# JUV BIRDS only JU, NE, VD

ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  filter(RegionShot!="NorthernAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="JU_PL", origin="min") %>%
  bind_rows(ALL_OUT)

ALL_OUT<-woco_juv %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  filter(RegionShot!="CentSouthAlps") %>%
  filter(RegionShot!="NorthernAlps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="juv", shot_in="JU_PL", origin="max") %>%
  bind_rows(ALL_OUT)








# ADULT BIRDS ALL CANTONS

ALL_OUT<-woco_ad %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="adult", shot_in="all", origin="min") %>%
  bind_rows(ALL_OUT)



ALL_OUT<-woco_ad %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="adult", shot_in="all", origin="max") %>%
  bind_rows(ALL_OUT)


# ADULT BIRDS only JU, NE, VD

ALL_OUT<-woco_ad %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  filter(RegionShot!="Alps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="adult", shot_in="JU_PL", origin="min") %>%
  bind_rows(ALL_OUT)



ALL_OUT<-woco_ad %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  filter(RegionShot!="Alps") %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot)) %>%
  mutate(age="adult", shot_in="JU_PL", origin="max") %>%
  bind_rows(ALL_OUT)



## CREATE SIMPLE GRAPH

ALL_OUT %>% filter(nonSwiss==1) %>% select(-Ntot) %>%
  spread(key=origin,value=prop) %>%
  ggplot(aes(x=shot_in,colour=age)) +
  geom_errorbar(aes(ymin=min, ymax=max), width=0.1) +
  facet_wrap(~age, ncol=1) +
  scale_y_continuous(name="prop shot woodcock from abroad",limits=c(0,1)) +
  theme_classic() +
  theme(legend.position="none")

ggsave("output/WOCO_prop_shot_origin.jpg", width=6, height=9)



## EXPORT TABLE

TableS1<-ALL_OUT %>% filter(nonSwiss==0) %>% select(-Ntot) %>%
  spread(key=origin,value=prop) %>%
  mutate(propSUI=paste0(round(max*100,1)," - ", round(min*100,1),"%")) %>%
  select(-min, -max, -nonSwiss) %>%
  pivot_wider(names_from=shot_in, values_from=propSUI)


fwrite(TableS1,"output/SUI_prop_shot_woodcock.csv")

