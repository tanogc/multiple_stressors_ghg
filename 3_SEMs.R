# Setting working folder
setwd("C:/Users/tanog/Dropbox/Trabajos en curso/streameco_ghg/script/def")

# Setting results plots
dir.create("plots")
dir.create("plots/assum")
plot_folder<-paste(getwd(),"/plots/",sep="")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")

# Loading packages
library(usdm)
library(piecewiseSEM)
library(MuMIn)
library(semEff)

# Loading data
dat<-read.table("prepared_dat.txt", sep="\t", h=T)

# Removing NAs
dat.m<-dat[which(is.na(dat$GPP.mle)==F),]

# Evaluating predictor collinearity (SEMs)
vifstep(dat.m[,c("din","discharge","tmean")])

vifstep(dat.m[,c("algal.production","din","discharge","flow_vel")])
vifstep(dat.m[,c("GPP.mle","din","discharge","flow_vel")])
vifstep(dat.m[,c("ER.mle","din","discharge","flow_vel")])
vifstep(dat.m[,c("NEP.mle","din","discharge","flow_vel")])

vifstep(dat.m[,c("do.def","din","discharge","flow_vel")])

################################################################
# Evaluating the role of local and distal processes (SEMs)
################################################################

############
# CO2
############

sem1 <- summary(mod1<-psem(lm(NEP.mle ~ din+discharge+rip_open, dat.m),
                           lm(pCO2 ~ NEP.mle+din+discharge, dat.m)))

# AICc values
AIC_sem<-unlist(as.numeric(c(sem1$AIC[2])))

# Model results
sem_res<-data.frame(mod=1, 
                    rbind(sem1$Cstat), 
                    AICc=round(AIC_sem,2),
                    r2=round(rbind(rsquared(mod1)$R.squared),2))
sem_res

# Calculating total standardized effects

# Direct effects
sem_distal<-data.frame(din=sem1$coef[5,8],
                       tmean=0,
                       rip_open=sem1$coef[3,8],
                       discharge=sem1$coef[6,8])

# localect effects
sem_local<-data.frame(din=sem1$coef[1,8]*sem1$coef[4,8],
                      tmean=0,
                      rip_open=0,
                      discharge=sem1$coef[2,8]*sem1$coef[4,8])

ste_co2<-as.matrix(rbind(sem_distal,sem_local))
rownames(ste_co2)<-c("distal", "local")

############
# CH4
############

sem1 <- summary(mod1<-psem( lm(NEP.mle~din+discharge+tmean, dat.m), 
                            lm(do.def ~ NEP.mle+din+discharge+tmean, dat.m), 
                            lm(pCH4 ~ do.def+din+discharge+tmean, dat.m)))

# AICc values
AIC_sem<-unlist(as.numeric(c(sem1$AIC[2])))

# Model results
sem_res<-data.frame(mod=1, 
                    rbind(sem1$Cstat), 
                    AICc=round(AIC_sem,2),
                    r2=round(rbind(rsquared(mod1)$R.squared),2))
sem_res

# Calculating total standardized effects

sem_distal<-data.frame(din=sem1$coef[9,8]+(sem1$coef[5,8]*sem1$coef[8,8]),
                    tmean=sem1$coef[11,8]+(sem1$coef[7,8]*sem1$coef[8,8]),
                    rip_open=0,
                    discharge=sem1$coef[10,8]+(sem1$coef[6,8]*sem1$coef[8,8]))

sem_local<-data.frame(din=((sem1$coef[1,8]*sem1$coef[4,8]*sem1$coef[8,8])),
                      tmean=((sem1$coef[3,8]*sem1$coef[4,8]*sem1$coef[8,8])),
                      rip_open=0,
                      discharge=((sem1$coef[2,8]*sem1$coef[4,8]*sem1$coef[8,8])))

ste<-sem_distal+sem_local

ste_ch4<-as.matrix(rbind(sem_distal,sem_local))
rownames(ste_ch4)<-c("distal", "local")

ste_ch4<-as.matrix(rbind(sem_distal,sem_local))

# Plotting standardized total effects

pdf(file=paste(plot_folder,"sem_effects.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=8*2,height=4*2)

par(mfrow=c(1,2),cex.axis=1.5*2, cex.lab=1.6*2, mar=c(5,5,7,5))

barplot(abs(ste_co2[,4:1]),col=c("gold","light blue"), las=1, xlim=c(0,1),names.arg=c("discharge","canopy open","temp.","DIN"),border=NA,horiz=T)
abline(v=0)

barplot(abs(ste_ch4[,4:1]),col=c("gold","light blue"), las=1, xlim=c(0,1),names.arg=c("discharge","canopy open","temp.","DIN"),border=NA,horiz=T)
abline(v=0)

dev.off()


