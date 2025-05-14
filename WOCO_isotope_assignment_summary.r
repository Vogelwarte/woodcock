### SUMMARY OF ISOTOPE ASSIGNMENTS ----
## Table 14 in Bohnenstengel et al. Report
## digitized by Steffen Oppel on 8 May 2025

library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
library(readxl)

woco_n<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table11_SampleSize")
woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")
woco_juv<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table15_IsotopeAssignments_JUV")
woco_ad<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table17_IsotopeAssignments_AD")
sum(woco$N)/9260


## total number shot and % of sample (reconstructed from report)

woco %>% group_by(RegionShot) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=c((0.07+0.167+0.149)/3,0.13,0.215)) %>% ## GUESSWORK ASSIGNMENT which canton represents plateau (Fribourg?) and Jura (JU,NE,VD)
  mutate(Nshot=as.integer(Ntot/prop)) %>%
  adorn_totals()
332/2224




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
  geom_errorbar(aes(ymin=min, ymax=max)) +
  facet_wrap(~age, ncol=1) +
  scale_y_continuous(limits=c(0,1)) +
  
  theme_classic()



