### ISOTOPic AND SHOOTING TIME ASSIGNMENT OF WOODCOCK FEATHERS ----
## to properly analyse the origins of birds in Bohnenstengel et al. Report
## initiated by Steffen Oppel on 16 May 2025


rm(list=ls())
library(data.table)
library(dplyr)
library(tidyverse)
library(janitor)
library(terra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(coda)
library(MCMCvis)
library(nimble)
library(tictoc)
library(basicMCMCplots) # for trace plots called chainsPlot
filter<-dplyr::filter
select<-dplyr::select
library(gridExtra)
library(tidyterra)


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Woodcock"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/Woodcock"),silent=T)



# 1. READ IN PROCESSED ISOTOPE DATA ----------------
load("data/woco.input.data.RData")




# 2. COMBINED PROBABILITY MODEL TO ESTIMATE PROB OF LOCAL ORIGIN IN NIMBLE ----

woco.orig.model<-nimbleCode({
  
  # Parameters:
  # p.nonlocal: probability of local origin (that shot woodcock was a local bird) - varies by age and canton
  
  # b.age: regression parameter for age that translates rainfall d2H to feather d2H 
  # b.rain: regression parameter for rainfall d2H that translates rainfall d2H to feather d2H
  # int.rain: regression intercept that translates rainfall d2H to feather d2H 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CALIBRATION REGRESSION FOR KNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Likelihood:  
  for (i in 1:nind.known){
        d2H_feather.known[i] ~ dnorm(mu.known[i], sd=sigma.calib)
        mu.known[i] <- int.rain + b.age*age.known[i] + b.rain*d2H_rain.known[i]
      }
    
  # Priors informed by known origin woodcock feathers in UK (Powell thesis 2012)
  int.rain ~ dnorm(4.5, sd=5) # informative prior based on Powell
  b.age ~ dnorm(-28.7, sd=5) # informative prior based on Powell
  b.rain ~ dnorm(0.8, sd=1) # informative prior based on Powell
  sigma.calib ~ dunif(0,20) # standard deviation
  dispersion ~ dnorm(15,tau=100) # dispersion parameter to convert prior probability into beta distribution - almost fixed quantity
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # TIME DISTRIBUTION FOR OVERALL ABUNDANCE
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # negative binomial distribution:
  for (i in 1:n.count.weeks){
    abund.known[i] ~ dnegbin(mu.abd.known[i], r.abd)
    mu.abd.known[i] <- int.abd + b.countday*count.day[i] + b2.countday*count.day.sq[i]
  }

  # Priors not informed by anything
  int.abd ~ dnorm(0.2, sd=5) # uninformative prior
  b.countday ~ dnorm(1, sd=5) # uninformative prior
  b2.countday ~ dnorm(1, sd=5) # uninformative prior for quadratic effect
  r.abd ~ dgamma(0.01,0.01) # dispersion parameter for neg bin distribution

  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PROBABILITY ESTIMATION FOR SHOT UNKNOWN ORIGIN BIRDS 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # Likelihood 
  sd.unknown[2]<-sigma.calib+sd.rain.d2H  ### overall distribution across Europe
  sd.unknown[1]<-sigma.calib
  
  for (i in 1:nind.unkn){
    
    # Potential local distribution based on timing used as prior
    logit.mig.unknown[i] <- lm.mean.mig +      ### intercept from migration model
      b.mig.week*(unk.week[i])     ### week effect from migration model
    p.nonlocal.prior1[i] <-ilogit(logit.mig.unknown[i]) ### predicted probability from departure model
    
    
    # Potential local distribution based on timing
    abd.unknown[i] <- int.abd + b.countday*unk.day[i] + b2.countday*unk.day.sq[i]    ### overall abundance distribution from abundance data
    p.nonlocal.prior2[i] <-ilogit(abd.unknown[i])   ### predicted probability from abundance time model
    
    
    # combine the probabilities of migration and abundance into a single probability using log odds
    # Convert each probability to odds:
    odds_mig[i] <- p.nonlocal.prior1[i] / (1 - p.nonlocal.prior1[i])
    odds_abd[i] <- p.nonlocal.prior2[i]  / (1 - p.nonlocal.prior2[i])
    
    # Combine odds
    combined_odds[i] <- odds_mig[i] * odds_abd[i]
    
    # Convert back to probability
    p.nonlocal.prior[i] <- combined_odds[i] / (1 + combined_odds[i])


    # Latent indicator for isotope distribution based on beta distribution of shooting time prior
    alpha[i] <- p.nonlocal.prior[i] * dispersion
    beta[i] <- (1 - p.nonlocal.prior[i]) * dispersion
    p.nonlocal[i] ~ dbeta(alpha[i], beta[i])
    z[i] ~ dbern(p.nonlocal[i])  ## 
    
    # Potential local distribution based on isotope ratio
    mu.unknown[i,2] <- int.rain + b.age*age.unknown[i] + b.rain*mean.rain.d2H    ### overall distribution across Europe
    mu.unknown[i,1] <- int.rain + b.age*age.unknown[i] + b.rain*d2H_rain.unknown[i]  ### local expected feather isotope distribution

    # evaluate origin feather isotope value from plausible target distributions
    d2H_feather.unknown[i] ~ dnorm(mu.unknown[i,z[i] + 1], sd=sd.unknown[z[i]+1])
    
    
  } #i
  

}) ## end of nimble code chunk






# 3. prepare the data needed for NIMBLE input --------
woco.unk.sf <- woco.unk.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_GS))

woco.sf <- woco.sf %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(dH)) %>%
  filter(!is.na(d2h_GS))

table(woco.unk.sf$AGE,woco.unk.sf$KANTON)

## 3.1. DISTINGUISH CONSTANTS AND DATA----
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
iso.constants <- list(nind.unkn = dim(woco.unk.sf)[1],
                      nind.known = dim(woco.sf)[1],
                      mean.rain.d2H = mean.rain.d2H,
                      sd.rain.d2H = sd.rain.d2H,
                      d2H_rain.known=woco.sf$d2h_GS,
                      d2H_rain.unknown=woco.unk.sf$d2h_GS,
                      age.unknown = ifelse(woco.unk.sf$AGE=="Adulte",0,1),
                      age.known = ifelse(woco.sf$AGE=="Adulte",0,1),
                      lm.mean.mig=logit(readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[2,1]),
                      b.mig.week=readRDS("output/woco_mig_depart_output_nimble.rds")$summary$all.chains[1,1],
                      unk.day=yday(UNK_WC$DATE),
                      unk.day.sq=(yday(UNK_WC$DATE))^2,
                      unk.week=week(UNK_WC$DATE)-week(ymd("2024-07-26")),   ## migration weeks start only in August
                      count.day=yday(woco_abundance$Date),
                      count.day.sq=yday(woco_abundance$Date)^2,
                      n.count.weeks=length(woco_abundance$SOPM)
                      )

iso.data <- list(d2H_feather.known = woco.unk.sf$dH,
                 abund.known=woco_abundance$SOPM,
                 d2H_feather.unknown = woco.unk.sf$dH)




## 3.2. specify NIMBLE run settings ----

# Parameters monitored
parameters.iso <- c("b.rain","b.age","int.rain","sigma.calib","dispersion","p.nonlocal","p.nonlocal.prior","p.nonlocal.prior1","p.nonlocal.prior2")

# Initial values  FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

iso.inits <- list(z = ifelse(woco.unk.sf$dH < 4.5+0.8*woco.unk.sf$d2h_GS-28*ifelse(woco.unk.sf$AGE=="Adulte",0,1),0,1),
                  int.rain = rnorm(1,4.5, sd=5), # informative prior based on Powell
                  b.age = rnorm(1,-28.7, sd=5), # informative prior based on Powell
                  b.rain = rnorm(1,0.8, sd=1), # informative prior based on Powell
                  dispersion = 15, # dispersion parameter for beta distribution set to sensible value
                  
                  ## INITS FOR ABUNDANCE REGRESSION
                  int.abd = rnorm(1,0.2, sd=5), # uninformative prior
                  b.countday = rnorm(1,1, sd=5), # uninformative prior
                  b2.countday = rnorm(1,1, sd=5), # uninformative prior for quadratic effect
                  r.abd  = rgamma(1,0.01,0.01), # dispersion parameter for neg bin distribution
                  
                  sigma.calib = runif(1,0,20) # standard deviation

)


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 50000
n.burnin <- 25000
n.chains <- 4



# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = woco.orig.model,
                    constants=iso.constants,
                    data = iso.data,
                    inits = iso.inits,
                    calculate=TRUE)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
test$initializeInfo()
#help(modelInitialization)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits



## 3.3. run model in NIMBLE ------------------------------------------
### this takes 400 sec for 50000 iterations and converges in that time
tic()

woco.iso <- nimbleMCMC(code = woco.orig.model,
                        constants=iso.constants,
                        data = iso.data,
                        inits = iso.inits,,
                        monitors = parameters.iso,
                        thin=4,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains,
                        progressBar = getNimbleOption("MCMCprogressBar"),
                        summary=T)
toc() 



saveRDS(woco.iso,"output/woco_iso_time_origin_model.rds")
#woco_surv<-readRDS("output/woco_iso_origin_model.rds")





# 4. assess model and examine convergence ------------------------------------------

out<- as.data.frame(MCMCsummary(woco.iso$samples, params=c("b.rain","b.age","int.rain","sigma.calib","dispersion"))) #"int.abd","b.countday","b2.countday","r.abd")))
out$parameter<-row.names(out)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
#out<-out %>%  select(parameter,Mean, median, lcl, ucl,SSeff,psrf)
out
fwrite(out,"output/woco_iso_time_origin_parm_estimates.csv")
#out<-fread("output/woco_iso_origin_parm_estimates.csv")


#MCMCplot(woco.iso$samples, params=c("mean.p.nonlocal"))


## look at the chains and whether they mixed well
chainsPlot(woco.iso$samples,
           var=c("b.rain","b.age","int.rain","sigma.calib","dispersion"))



## 4.1 summarise probabilities across individuals ------------------------------------------

# compile all the samples
samples <- rbind(woco.iso$samples$chain1,woco.iso$samples$chain2,woco.iso$samples$chain3,woco.iso$samples$chain4)
mean.p.nonlocal <- as_tibble(samples[,grep("p.nonlocal", colnames(samples))]) %>%
  gather(key="parameter",value="p.nonlocal") %>%
  mutate(ind=as.numeric(str_extract(parameter,pattern="\\d+"))) %>%
  mutate(age=iso.constants$age.unknown[ind]) %>%
  mutate(ctn=woco.unk.sf$KANTON[ind]) 

hist(mean.p.nonlocal$p.nonlocal)






## 4.2. summarise output in graphical form ------------------------------------------

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


FIGURE<-mean.p.nonlocal %>%
  group_by(age,ctn,ind) %>%
  summarise(p.nonlocal.mean=mean(p.nonlocal)) %>%
  ungroup() %>%
  group_by(age,ctn) %>%
  summarise(for.med=median(p.nonlocal.mean),for.ucl=quantile(p.nonlocal.mean,0.025), for.lcl=quantile(p.nonlocal.mean,0.975)) %>%
  
  mutate(Age=ifelse(age==1,"Adult","Juvenile")) %>%
  #mutate(Kanton=levels(as.factor(woco.unk.sf$KANTON))[ctn]) %>%
  
  ggplot(aes(x=ctn, y=for.med))+
  geom_point(aes(col=Age), position=position_dodge(width=0.2), size=2.5) +
  geom_errorbar(aes(ymin=for.lcl, ymax=for.ucl, col=Age), width=0.05, linewidth=1, position=position_dodge(width=0.2)) +
  
  ## format axis ticks
  labs(y="Prop. of shot woodcock that are not local",x="Swiss Canton",col="") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,1,0.2)) +
  # scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  # scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  # 
  ### add the bird icons
  annotation_custom(gunicon, xmin=5.2, xmax=5.8, ymin=0.88, ymax=0.96)+
  annotation_custom(wocoicon, xmin=5.5, xmax=6.9, ymin=0.85, ymax=1.02)+
  
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
        legend.position.inside=c(0.65,0.85),
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
FIGURE

ggsave(plot=FIGURE,
       filename="output/woco_iso_time_origin_probability_estimates.jpg", 
       device="jpg",width=11, height=8)




# 5. CREATE MAPS FOR EACH CANTON WHAT LOCAL RAINFALL ENCOMPASSES -----------
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

grid.arrange(grobs=plot_list,ncol=2)
ggsave(filename="output/woco_iso_origin_local_maps.jpg", 
       device="jpg",width=8, height=12)

