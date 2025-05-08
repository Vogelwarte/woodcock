### SUMMARY OF ISOTOPE ASSIGNMENTS ----
## Table 14 in Bohnenstengel et al. Report
## digitized by Steffen Oppel on 8 May 2025

library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)

woco<-fread("output/Table14_IsotopeAssignments_HuntingBag.csv")
sum(woco$N)/9260


## total number shot and % of sample (reconstructed from report)

woco %>% group_by(RegionShot) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=c((0.07+0.167+0.149)/3,0.13,0.215)) %>% ## GUESSWORK ASSIGNMENT which canton represents plateau (Fribourg?) and Jura (JU,NE,VD)
  mutate(Nshot=as.integer(Ntot/prop)) %>%
  adorn_totals()
332/2224

## summarise proportion of birds originating from outside Switzerland

woco %>% mutate(nonSwiss=if_else(OriginRegion>6,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot))


woco %>% mutate(nonSwiss=if_else(OriginRegion>5,1,0)) %>%
  group_by(nonSwiss) %>%
  summarise(Ntot=sum(N)) %>%
  mutate(prop=Ntot/sum(Ntot))