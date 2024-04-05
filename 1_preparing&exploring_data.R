# Setting working folder
setwd("C:/Users/tanog/Dropbox/Trabajos en curso/streameco_ghg/script/def")

# Setting results plots
dir.create("plots")
dir.create("plots/assum")
plot_folder<-paste(getwd(),"/plots/",sep="")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")

# Loading packages
library(usdm)
library(corrplot)
library(ade4)
library(viridis)

# Loading data
dat<-read.table("dat_completed.txt", sep="\t", h=T)
dat_raw<-dat # saving raw values

# Removing sites with bad metabolic estimations
dat[c(19, 28, 38, 41), c("GPP.mle","ER.mle")]<-NA

############################
# Variable transformation
############################

q<-sapply(dat,class)=="numeric" | sapply(dat,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(q==T)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

# Water chemistry - LUI
dat$din<-log(dat$din)
dat$do.def<-log(dat$do.def)
dat$co2per<-log(dat$pCO2/1000/dat$DOmean)

# Riparian canopy openness
dat$rip_open<-car::logit(1-dat$NDVI100m)

# Flow reduction
dat$discharge<-log(dat$discharge)
dat$flow_vel<-log(dat$flow_vel)

# Algal production
dat$algal.production<-log(dat$algal.production)

# Metabolism
dat$GPP.mle<-log(dat$GPP.mle)
dat$ER.mle<-log(-dat$ER.mle)
dat$NEP.mle<-(-log(-dat$NEP.mle))

# GHGs
dat$pCO2<-log(dat$pCO2)
dat$pCH4<-log(dat$pCH4)
dat$FCO2.h<-log(dat$FCO2.h+94)
dat$FCH4.h<-log(dat$FCH4.h)
dat$ch4co2_ratio<-sqrt(dat$ch4co2_ratio)

par(mfrow=c(2,3))
for (i in which(q==T)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

# Scaling variables
dat[,-1]<-scale(dat[,-1])

# Saving prepared data.frame for models
write.table(dat, "prepared_dat.txt",sep="\t")

############################
# Exploratory analysis
############################

# Characterizing stress gradients
summary(dat_raw$din)
summary(dat_raw$do.def)
summary(dat_raw$tmean)
summary(dat_raw$NDVI100m)
summary(dat_raw$discharge)

stressor_ranges1<-apply(dat_raw[,c("din","do.def", "tmean","NDVI100m","discharge")],2,median)
stressor_ranges2<-apply(dat_raw[,c("din","do.def", "tmean","NDVI100m","discharge")],2,range)

stressor_ranges1
stressor_ranges2

# Saving stressor ranges      
write.table(rbind(stressor_ranges1, stressor_ranges2), "stressor_range.txt",sep="\t")

# Relationship between K600 and ER

pdf(file=paste(plot_folder,"rel_k600_ER.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5,height=5.5)

cor1 <- bquote(italic(r[P]) == .(round(cor(log(-dat_raw$ER.mle), dat_raw$K600_hyd), 2)))

par(mfrow=c(1,1),cex=1.6, cex.lab=0.9, cex.axis=0.6)

# GPP vs ER
plot(dat_raw$K600_hyd~abs(dat_raw$ER.mle), xlab=expression(k[600]~(m~ d^{-1})), 
     ylab=expression(ER~(mmol~ C~ m^{-2} ~ d^{-1})),pch=21,bg="red")

mod<-lm(K600_hyd~-ER.mle, data=dat_raw, na.action="na.omit")

mtext(cor1,  line = -1.35, adj = 0.75, cex = 1.3, font = 2)

abline(mod, col="red", lwd=3)

dev.off()

# Histogram of metabolic rates, algal production and GHGs

var_names<-c("GPP", "ER","NEP","Algal production", "pCO2" ,"pCH4")

var_col<-c("red","red","red","orange","blue","green")

pdf(file=paste(plot_folder,"hist_ghg.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=9,height=6)

par(mfrow=c(2,3),cex.lab=1.4,cex.axis=1.2) 

hist(dat_raw$GPP.mle, col="red", main="", xlab=expression(GPP~(mmol~ C~ m^{-2} ~ d^{-1})))

hist(-dat_raw$ER.mle, col="red", main="", ylab="",xlab=expression(ER~(mmol~ C~ m^{-2} ~ d^{-1})))

hist(dat_raw$NEP.mle, col="red", main="", ylab="",xlab=expression(NEP~(mmol~ C~ m^{-2} ~ d^{-1})))

hist(dat_raw$algal.production, col="orange", main="", xlab=expression(Algal~production~(Chl ~italic("a")~ m^{-2} ~ d^{-1})))

hist(dat_raw$pCO2, col="blue", main="", ylab="",xlab=expression(italic("p")~ "CO"[2] ~ (mu ~ atm)))

hist(dat_raw$pCH4, col="green", main="", ylab="",xlab=expression(italic("p")~ "CH"[4] ~ (mu ~ atm)))

dev.off()

# Histogram of K600 values

pdf(file=paste(plot_folder,"hist_k600.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=4,height=4)

par(mfrow=c(1,1),cex.lab=1.4,cex.axis=1.2) 

hist(dat_raw$K600_hyd, col="light blue", main="", xlab=expression(k[600]~(m~ d^{-1})))

dev.off()

# Relationship between metabolic rates, DO deficit and GHGs

pdf(file=paste(plot_folder,"rel_ghg.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=10,height=11)

cor1 <- bquote(italic(r[P]) == .(round(cor(dat.m$GPP.mle, dat.m$ER.mle), 2)))
cor2 <- bquote(italic(r[P]) == .(round(cor(dat.m$NEP.mle, dat.m$do.def), 2)))
cor3 <- bquote(italic(r[P]) == .(round(cor(dat$pCO2, dat$FCO2.h), 2)))
cor4 <- bquote(italic(r[P]) == .(round(cor(dat$pCH4, dat$FCH4.h), 2)))

par(mfrow=c(2,2), cex=1.8, cex.lab=1, cex.axis=0.6, mar=c(4,5,4,1))

# GPP vs ER
plot(-dat_raw$ER.mle~dat_raw$GPP.mle, ylab=expression(ER~(mmol~ C~ m^{-2} ~ d^{-1})), 
     xlab=expression(GPP~(mmol~ C~ m^{-2} ~ d^{-1})),pch=21,bg="red")

mod<-lm(-ER.mle~GPP.mle, data=dat_raw)

mtext(cor1,  line = -1.25, adj = 0.75, cex = 1.75, font = 2)

abline(mod, col="red", lwd=3)
abline(a=0, b=1, col="black", lwd=1, lty=2)


# NEP vs DO def

plot(dat_raw$do.def~dat_raw$NEP.mle, ylab=expression(DO~deficit~(mg~ DO~ L^{-1})), 
     xlab=expression(NEP~(mmol~ C~ m^{-2} ~ d^{-1})),pch=21,bg="red")

mod<-lm(do.def~NEP.mle, data=dat_raw)

mtext(cor2,  line = -1.25, adj = 0.95, cex = 1.75, font = 2)

abline(mod, col="red", lwd=3)

# pCO2 vs FCO2
# Ientifying negative fluxes
co2_col<-rep("blue",nrow(dat))
co2_col[which(dat_raw$FCO2.h<0)]<-"gold"

plot(dat_raw$FCO2.h~dat_raw$pCO2, pch=21, bg=co2_col, main="",ylab=expression(italic("F")~CO[2]~(mmol~ "CO"[2]~ m^{-2} ~ d^{-1})),
     xlab=expression(italic("p")~ "CO"[2] ~ (mu ~ atm)))

mod<-lm(FCO2.h~pCO2, dat_raw)

mtext(cor3,  line = -1.25, adj = 0.95, cex = 1.75, font = 2)

abline(mod, col="blue", lwd=3)
abline(h=0, col="black", lty=2)

# pCH4 vs FCH4

plot(dat_raw$FCH4.h ~ dat_raw$pCH4, pch=21, bg="green", main="",ylab=expression(italic("F")~CH[4]~(mmol~ "CO"[2]~ m^{-2} ~ d^{-1})),
     xlab=expression(italic("p")~ "CH"[4] ~ (mu ~ atm)))

mod<-lm(FCH4.h~pCH4, dat_raw)

mtext(cor4,  line = -1.25, adj = 0.95, cex = 1.75, font = 2)

abline(mod, col="green", lwd=3)

dev.off()




