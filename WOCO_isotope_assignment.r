### ISOTOPic ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025

## goal is to estimate proportion of shot birds coming from Switzerland



## USEFUL RESOURCES:
# workshop: https://github.com/JacksonKusack/Isotope-Assignment-Workshop
# paper: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/04-0175
# isocat package: https://github.com/cjcampbell/isocat



library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
#library(readxl)
library(assignR) ## https://cran.r-project.org/web/packages/assignR/vignettes/assignR.html
# library(isocat) ## https://cran.r-project.org/web/packages/isocat/vignettes/isocat.html
# library(isoAssign) ## https://rdrr.io/github/SMBC-NZP/MigConnectivity/man/isoAssign.html
library(terra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)



## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock")



# 1. READ IN PROCESSED ISOTOPE DATA ----

### dH_reg stands for rainfall d2H
### dH_correct is the corrected d2H of feathers to 


woco<-fread("data/WOCO_isotopes.csv")
#woco<-read_excel("output/IsotopeAssignment_Tables.xlsx", sheet="Table14_IsotopeAssignments_Hunt")


## 1.1. SPLIT INTO KNOWN ORIGIN AND UNKNOWN ----

ORIG_WC<-woco %>% filter(ORIGINE!="UNBEKANNT") %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,dH_reg,dH_scaled,dH_correct,PROVENANCE_voigt,LATITUDE,LONGITUDE)
UNK_WC<-woco %>% filter(ORIGINE=="UNBEKANNT") %>%
  dplyr::select(ID,ADULTE,AGE,KANTON,DATE,dH,dH_reg,dH_scaled,dH_correct,PROVENANCE_voigt)



## 1.2. PLOT THE DIFFERENT d2H isotope values ----

ggplot(woco, aes(x=dH,y=dH_reg, col=AGE)) +
  geom_point() +
  geom_smooth()


ggplot(woco, aes(x=dH_scaled,y=dH_correct, col=AGE)) +
  geom_point() +
  geom_smooth()



## 1.3. PLOT HISTOGRAMS FOR SUISSE AND OTHER BIRDS ----

ggplot(woco, aes(x=dH_correct, col=ORIGINE, fill=ORIGINE)) +
  geom_histogram(alpha=0.5,position = position_dodge(width=1)) 


distd2H_knownSUI<-hist(ORIG_WC$dH_correct)
distd2H_knownSUI$density   ## this can be used as prior information for the isotope assignment





# 2. Feather fractionation from rainwater hydrogen isotope ----

# converting isoscape rainfall values into feather d2H values
# based on Powell 2022, Table 3.1, page 133 with a sample size of 135 feathers
n <- 135
intcept<- rnorm(100,4.50799, (6.34005*sqrt(n)))
beta.rain<- rnorm(100,0.79760, (0.08049*sqrt(n)))
beta.age<- rnorm(100,-28.73077, (2.14900*sqrt(n)))


#woco$dH <- intcept + beta.rain*d2Hrain + beta.age * juvenile  ## the transformation equation of Powell based on mean annual deuterium (MAD)
woco$dH_Powell <- 4.50799 + 0.79760*woco$dH_reg - 28.73077 * woco$JUVENIL


ggplot(woco, aes(x=dH_reg,y=dH_Powell, col=AGE)) +
  geom_point() +
  geom_smooth()

ggplot(woco, aes(x=dH_correct,y=dH_Powell, col=AGE)) +
  geom_point() +
  geom_smooth()


ggplot(woco, aes(x=dH_scaled,y=dH_Powell, col=AGE)) +
  geom_point() +
  geom_smooth()


## try and figure out how values were created
#summary(lm(dH_correct~dH+AGE, data=woco))
summary(lm(dH_correct~dH_reg+AGE, data=woco))  ### this seems to be the appropriate equation
#summary(lm(dH_correct~dH_scaled+AGE, data=woco))
summary(lm(dH_Powell~dH_reg+AGE, data=woco))




## 2.1. building our own equation using the known origin feather samples ----
## this looks weird like somebody had calculated the dH_reg
summary(lm(dH~dH_reg+AGE, data=ORIG_WC))






# 3. Geographic assignment using assignR ----

## 3.1. load shapefile of Switzerland and EUROPE ----

SUI <- ne_countries(country = "Switzerland", scale=10, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84
plot(SUI)

bbox <- st_sfc(st_point(c(-12, 35)), st_point(c(45, 75)), crs = 4326) %>% st_bbox()
EUR <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(bbox) %>%
  st_transform(st_crs('EPSG:4326')) %>% # Project to WGS84
  dplyr::select(admin,name,adm0_a3,geometry)
plot(EUR)




## 3.2. plot the known origin points of WOCO that were sampled ----

woco.sf <- ORIG_WC %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs=4326)


ggplot(SUI) +
  geom_sf() +
  geom_sf(data=woco.sf,color="red")





## 3.3. extract isoscape ----

### this function throws a weird error - downloaded and ran function independently without problem
# isoscape <- getIsoscapes(isoType = "USTap", timeout = 1200) %>%   ## we use MA because that is what Powell's conversion is based on
#   projectRaster(crs = CRS(SRS_string = 'EPSG:4326'))
# 
# isoscape <- isoscape [[1:2]]  ## retain only dH
# saveRDS(isoscape,"data/global_d2H_MA_isoscape.rds")
isoscape <- readRDS("data/global_d2H_MA_isoscape.rds") %>%
  crop(extent(EUR))

plot(isoscape, xlab="Longitude", ylab="Latitude")





## 3.4. get canned calibration data and calibrate rainwater against feather d2H  ----

# data("knownOrig")
# knownOrig$sources[1:2]
# knownOrig$samples %>% filter(Group=="Ground bird")
# 
# cal.grd.bird <- subOrigData(marker = 'd2H', dataset = c(5), ref_scale = NULL) # Extract data
# str(cal.grd.bird$data)
# 
# d2Hcalib <- calRaster(known = cal.grd.bird, isoscape = isoscape, interpMethod = 1, verboseLM = F, genplot = F) # Calibrate the mean and sd values from our isoscape
# summary(d2Hcalib$lm.model) # Extract model results
#  
# ggplot(data = d2Hcalib$lm.data, aes(y = tissue.iso, x = isoscape.iso)) + 
#   geom_point()  + 
#   stat_smooth(method = "lm", formula = 'y ~ x') + 
#   theme_classic()

## manually rescale the isoscape with the Powell equation

isoscape.AD<- 4.50799 + 0.79760*isoscape
isoscape.JUV<- 4.50799 + 0.79760*isoscape - 28.73077

plot(isoscape.AD, xlab="Longitude", ylab="Latitude")






## 3.5. ASSIGNMENT TO ORIGIN  ----

### because the isotope data are already converted into rainfall isotopes, we do not need the calibrated isoscape
origins <- pdRaster(isoscape, UNK_WC[,c(1,9)], mask = as(EUR, 'Spatial'), genplot = F)















# 4. Geographic assignment using assignR ----

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




# 5. SUGGESTIONS BY COPILOT WITHOUT CONSIDERING ANY ISOSCAPE ----



To estimate the probability that birds come from Switzerland based on isotope ratios in their feathers, you can use a statistical approach known as **Bayesian inference**. Here's a step-by-step outline of how you might proceed:

1. **Collect Data**: 
   - Obtain the isotope ratios (e.g., hydrogen, strontium) in feathers from birds of known Swiss origin.
   - Collect isotope ratios from feathers of birds with unknown origins.

2. **Model the Data**:
   - Assume that the isotope ratios from Swiss birds follow a certain probability distribution (e.g., normal distribution). This distribution is your **prior** knowledge.
   - The isotope ratios from the unknown birds are your **observed data**.

3. **Calculate Likelihood**:
   - For each bird of unknown origin, calculate the likelihood of observing its isotope ratio given the distribution of isotope ratios from Swiss birds. This involves using the probability density function of the assumed distribution.

4. **Apply Bayes' Theorem**:
  - Use Bayes' Theorem to update your prior knowledge with the observed data to get the **posterior probability**. Bayes' Theorem is given by:
  $$
  P(\text{Swiss} \mid \text{Data}) = \frac{P(\text{Data} \mid \text{Swiss}) \cdot P(\text{Swiss})}{P(\text{Data})}
$$
  where:
  - \( P(\text{Swiss} \mid \text{Data}) \) is the posterior probability that the bird is from Switzerland given the observed data.
- \( P(\text{Data} \mid \text{Swiss}) \) is the likelihood of observing the data given that the bird is from Switzerland.
- \( P(\text{Swiss}) \) is the prior probability that a bird is from Switzerland.
- \( P(\text{Data}) \) is the marginal likelihood of the data.

5. **Interpret Results**:
  - The posterior probability gives you an updated estimate of the likelihood that the bird is from Switzerland based on its isotope ratios.

### Example Calculation
Let's say you have a normal distribution for the isotope ratios of Swiss birds with mean \( \mu \) and standard deviation \( \sigma \). For an unknown bird with isotope ratio \( x \), the likelihood \( P(x \mid \text{Swiss}) \) can be calculated using the normal distribution's probability density function:
  $$
  P(x \mid \text{Swiss}) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)
$$
  
  You can then plug this into Bayes' Theorem to get the posterior probability.

### Tools and Software
- **R** or **Python**: Both have libraries for Bayesian analysis (e.g., `rjags` in R, `PyMC3` in Python).
- **Stable Isotope Mixing Models**: Tools like IsoSource or MixSIAR can help in analyzing isotope data [1](https://www.epa.gov/eco-research/stable-isotope-mixing-models-estimating-source-proportions).

Would you like more detailed guidance on any of these steps?


To incorporate an informative prior based on the isotope ratios from Swiss birds into your Nimble model, you can use the data from Swiss birds to define the prior distributions for the parameters (mean `mu` and standard deviation `sigma`). Here's how you can modify the previous code to include this informative prior:

1. **Calculate the Prior Parameters**: Use the isotope ratios from Swiss birds to estimate the mean and standard deviation. These estimates will be used to define the prior distributions.

2. **Modify the Model**: Update the prior distributions for `mu` and `sigma` based on the calculated parameters.

Here's an example:
  
  ```r
library(nimble)

# Simulated data: isotope ratios for birds of Swiss origin
isotope_ratio_swiss <- c(1.1, 1.3, 1.4, 1.6, 1.7)

# Calculate the mean and standard deviation from Swiss data
mu_prior <- mean(isotope_ratio_swiss)
sigma_prior <- sd(isotope_ratio_swiss)

# Define the model with informative priors
code <- nimbleCode({
  for (i in 1:N) {
    isotope_ratio[i] ~ dnorm(mu, sd = sigma)
  }
  mu ~ dnorm(mu_prior, sd = sigma_prior)
  sigma ~ dunif(0, 50)
})

# Simulated data: isotope ratios for birds of unknown origin
isotope_ratio_unknown <- c(1.2, 1.5, 1.8, 2.0, 2.3)

# Number of observations
N <- length(isotope_ratio_unknown)

# Define constants and data
constants <- list(N = N, mu_prior = mu_prior, sigma_prior = sigma_prior)
data <- list(isotope_ratio = isotope_ratio_unknown)

# Initial values for the parameters
inits <- list(mu = mu_prior, sigma = sigma_prior)

# Create the model
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configure and build MCMC
mcmc_conf <- configureMCMC(model)
mcmc <- buildMCMC(mcmc_conf)

# Compile the model and MCMC
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model)

# Run the MCMC
samples <- runMCMC(Cmcmc, niter = 10000)

# Summary of the posterior samples
summary(samples)
```

### Explanation:
1. **Prior Parameters**: The mean (`mu_prior`) and standard deviation (`sigma_prior`) are calculated from the isotope ratios of Swiss birds.
2. **Model Definition**: The prior for `mu` is now a normal distribution with mean `mu_prior` and standard deviation `sigma_prior`. The prior for `sigma` remains a uniform distribution.
3. **Constants**: The calculated prior parameters are included in the constants list.
4. **Initial Values**: The initial values for `mu` and `sigma` are set to the calculated prior parameters.

This approach ensures that the prior information from Swiss birds is incorporated into the model, providing a more informed estimate of the probability that the birds of unknown origin are from Switzerland.

If you have any further questions or need additional adjustments, feel free to ask!