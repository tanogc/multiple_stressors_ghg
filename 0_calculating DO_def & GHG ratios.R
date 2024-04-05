# Setting working folder
setwd("C:/Users/tanog/Dropbox/Trabajos en curso/streameco_ghg/script/def")

# Loading packages
library(streamMetabolizer)

# Loading data
dat<-read.table("dat.txt", sep="\t", h=T)

# Estimating DO deficit 

# Estimating air pressure
air_temp<-dat$air_mean_temp # Air temperature
alt<-dat$alt # altitude
air_pres<-calc_air_pressure(temp.air = air_temp, elevation = alt) # air pressure
wt<-dat$tmean # water temperature

do.sat <- calc_DO_sat(temp.water = wt, pressure.air = air_pres) # saturation concentration
dat$do.def<-do.sat-dat$DOmean

# CH4:CO2 molar ratio
pCH4_mmol<-dat$pCH4/(12+1*4)
pCO2_mmol<-dat$pCO2/(12+16*2)
dat$ch4co2_ratio<-pCH4_mmol/pCO2_mmol

# Calculating contribution of CH4 to total warming potential
dat$pceq<-28*dat$pCH4+dat$pCO2 # Calculating CO2-eq concentration
dat$perCH4<-28*dat$pCH4/dat$pceq

summary(dat$ch4co2_ratio)
hist(dat$ch4co2_ratio)

dat$NEP.mle<-dat$GPP.mle+dat$ER.mle

# Saving dataset
write.table(dat, "dat_completed.txt",sep="\t")
