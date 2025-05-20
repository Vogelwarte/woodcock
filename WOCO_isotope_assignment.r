### ISOTOPic ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025

## goal is to estimate proportion of shot birds coming from Switzerland

## NEED A TABLE WITH RAW DATA

library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
#library(readxl)
library(assignR) ## https://cran.r-project.org/web/packages/assignR/vignettes/assignR.html
# library(isocat) ## https://cran.r-project.org/web/packages/isocat/vignettes/isocat.html
# library(isoAssign) ## https://rdrr.io/github/SMBC-NZP/MigConnectivity/man/isoAssign.html
library(terra)



## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")



# 1. READ IN PROCESSED ISOTOPE DATA ----

woco<-fread("data/WOCO_isotopes.csv")
#woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")


## 1.1. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH_reg,dH_scaled,dH_correct,PROVENANCE_voigt)
UNK_WC<-woco %>% filter(ORIGINE=="UNBEKANNT") %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH_reg,dH_scaled,dH_correct,PROVENANCE_voigt)



## 1.2. PLOT THE DIFFERENT d2H isotope values ----

ggplot(woco, aes(x=dH_reg,y=dH_scaled, col=AGE)) +
  geom_point() +
  geom_smooth()


ggplot(woco, aes(x=dH_reg,y=dH_correct, col=AGE)) +
  geom_point() +
  geom_smooth()



## 1.3. PLOT HISTOGRAMS FOR SUISSE AND OTHER BIRDS ----

ggplot(woco, aes(x=dH_correct, col=ORIGINE, fill=ORIGINE)) +
  geom_histogram(alpha=0.5,position = position_dodge(width=1)) 






# 2. Feather fractionation from rainwater hydrogen isotope ----

# converting isoscape rainfall values into feather d2H values
# based on Powell 2022, Table 3.1, page 133 with a sample size of 135 feathers
n <- 135
intcept<- rnorm(100,4.50799, (6.34005*sqrt(n)))
beta.rain<- rnorm(100,0.79760, (0.08049*sqrt(n)))
beta.age<- rnorm(100,-28.73077, (2.14900*sqrt(n)))


d2Hfeather <- intcept + beta.rain*d2Hrain + beta.age * juvenile



## 2. Geographic assignment using assignR ----

plot(naMap)
plot(d2h_lrNA)
names(knownOrig$sites)
names(knownOrig$samples)
plot(wrld_simpl)
points(knownOrig$sites, col = "red")

Ll_d = subOrigData(taxon = "Lanius ludovicianus", mask = naMap, ref_scale = NULL)  ## package data for shrikes
d2h_Ll = calRaster(known = Ll_d, isoscape = d2h_lrNA, mask = naMap)


## simulation of fake data
id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
d2H.sd = runif(5, 1.5, 2.5)
d2H_cal = rep("UT_H_1", 5)
Ll_un = data.frame(id, d2H, d2H.sd, d2H_cal)
print(Ll_un)

## convert if reference scales differ
Ll_un = refTrans(Ll_un, ref_scale = "OldEC.1_H_1")
print(Ll_un)

## posterior probability density map
Ll_prob = pdRaster(d2h_Ll, Ll_un)

## should sum to 1
global(Ll_prob[[1]], 'sum', na.rm = TRUE)


## odds ratio
s1 = states[states$STATE_ABBR == "UT",]
s2 = states[states$STATE_ABBR == "NM",]
plot(naMap)
plot(s1, col = c("red"), add = TRUE)
plot(s2, col = c("blue"), add = TRUE)

## assignment
qtlRaster(Ll_prob, threshold = 0.1)
qtlRaster(Ll_prob, threshold = 0.8, thresholdType = "prob")
jointP(Ll_prob)
Ll_up = unionP(Ll_prob)
qtlRaster(Ll_up, threshold = 0.1)
