# Setting working folder
setwd("C:/Users/tanog/Dropbox/Trabajos en curso/streameco_ghg/script/def")

# Setting results plots
dir.create("plots")
dir.create("plots/assum")
plot_folder<-paste(getwd(),"/plots/",sep="")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")

# Loading packages
library(MuMIn)
library(usdm)
library(variancePartition)
library(corrplot)
library(ade4)
library(viridis)
library(lattice)
source("0_part_r2.R")

# Loading data
dat<-read.table("prepared_dat.txt", sep="\t", h=T)

# Removing sites with bad metabolic estimations
dat.m<-dat[-c(19, 28, 38, 41),]

# Evaluating predictor collinearity
vifstep(dat[,c("din","do.def","tmean","discharge","rip_open")])

################################################################################
# Multiple stressor effects on metabolism and GHGs (Multi-model inference)
################################################################################

# Creating a list
# Creating lists to store model results
ci<-avg<-r2_order<-mod_set<-res_d<-res_d_full<-res<-list() # to store data

#Models for multiple stressors analysis
options (na.action="na.fail") # necessary to run dredge()

#######################################################################
# GPP
#######################################################################

i=1

#Global model
mod1<-lm(GPP.mle~din+tmean+discharge+rip_open,data=dat.m)

# Models including interactions
mod2<-lm(GPP.mle~din+tmean+discharge+rip_open+din:rip_open,data=dat.m)
mod3<-lm(GPP.mle~din+tmean+discharge+rip_open+din:tmean,data=dat.m)
mod4<-lm(GPP.mle~din+tmean+discharge+rip_open+din:discharge,data=dat.m)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

mod<-mod1 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_gpp <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_gpp)->n.mod

for (k in 1:n.mod) 
{
  mods_gpp[[k]]->mod
  pdf(file=paste(assum_folder,"gpp_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:5]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_gpp, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# ER
#######################################################################

i=2

#Global model
mod1<-lm(ER.mle~din+tmean+discharge+rip_open,data=dat.m)

# Models including interactions
mod2<-lm(ER.mle~din+tmean+discharge+rip_open+din:rip_open,data=dat.m)
mod3<-lm(ER.mle~din+tmean+discharge+rip_open+din:tmean,data=dat.m)
mod4<-lm(ER.mle~din+tmean+discharge+rip_open+din:discharge,data=dat.m)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

# Testing if mod2 has better explanatory capacity than mod1
anova(mod1, mod2)

mod<-mod1 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_er <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_er)->n.mod

for (k in 1:n.mod) 
{
  mods_er[[k]]->mod
  pdf(file=paste(assum_folder,"er_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:5]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_er, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# NEP
#######################################################################

i=3

#Global model
mod1<-lm(NEP.mle~din+tmean+discharge+rip_open,data=dat.m)

# Models including interactions
mod2<-lm(NEP.mle~din+tmean+discharge+rip_open+din:rip_open,data=dat.m)
mod3<-lm(NEP.mle~din+tmean+discharge+rip_open+din:tmean,data=dat.m)
mod4<-lm(NEP.mle~din+tmean+discharge+rip_open+din:discharge,data=dat.m)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

mod<-mod1 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_nep <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_nep)->n.mod

for (k in 1:n.mod) 
{
  mods_nep[[k]]->mod
  pdf(file=paste(assum_folder,"ner_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:5]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_nep, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# algal production
#######################################################################

i=4

#Global model
mod1<-lm(algal.production~din+tmean+discharge+rip_open,data=dat)

# Models including interactions
mod2<-lm(algal.production~din+tmean+discharge+rip_open+din:rip_open,data=dat)
mod3<-lm(algal.production~din+tmean+discharge+rip_open+din:tmean,data=dat)
mod4<-lm(algal.production~din+tmean+discharge+rip_open+din:discharge,data=dat)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

mod<-mod1 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_algal <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_algal)->n.mod

for (k in 1:n.mod) 
{
  mods_algal[[k]]->mod
  pdf(file=paste(assum_folder,"algal_prod_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:5]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_algal, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# pCO2
#######################################################################

i=5

#Global model
mod1<-lm(pCO2~din+tmean+discharge+rip_open,data=dat)

# Models including interactions
mod2<-lm(pCO2~din+tmean+discharge+rip_open+din:rip_open,data=dat)
mod3<-lm(pCO2~din+tmean+discharge+rip_open+din:tmean,data=dat)
mod4<-lm(pCO2~din+tmean+discharge+rip_open+din:discharge,data=dat)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

mod<-mod1 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_co2 <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_co2)->n.mod

for (k in 1:n.mod) 
{
  mods_co2[[k]]->mod
  pdf(file=paste(assum_folder,"pco2_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:5]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_co2, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# pCH4
#######################################################################

i=6

#Global model
mod1<-lm(pCH4~din+do.def+tmean+discharge+rip_open,data=dat)

# Models including interactions
mod2<-lm(pCH4~din+do.def+tmean+discharge+rip_open+din:rip_open,data=dat)
mod3<-lm(pCH4~din+do.def+tmean+discharge+rip_open+din:tmean,data=dat)
mod4<-lm(pCH4~din+do.def+tmean+discharge+rip_open+din:discharge,data=dat)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

# Comparing the explanatory capacity of mod1 and mod3
anova(mod1, mod3)

mod<-mod3 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_ch4 <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_ch4)->n.mod

for (k in 1:n.mod) 
{
  mods_ch4[[k]]->mod
  pdf(file=paste(assum_folder,"pch4_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:7]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_ch4, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

#######################################################################
# Molar ratio CH4:CO2
#######################################################################

i=7

#Global model
mod1<-lm(ch4co2_ratio~din+do.def+tmean+discharge+rip_open,data=dat)

# Models including interactions
mod2<-lm(ch4co2_ratio~din+do.def+tmean+discharge+rip_open+din:rip_open,data=dat)
mod3<-lm(ch4co2_ratio~din+do.def+tmean+discharge+rip_open+din:tmean,data=dat)
mod4<-lm(ch4co2_ratio~din+do.def+tmean+discharge+rip_open+din:discharge,data=dat)

# Identifying best model structure (with or without interactions)
AICc(mod1, mod2, mod3, mod4)

# Comparing the explanatory capacity of mod1 and mod3
anova(mod1, mod3)

mod<-mod3 # Choosing most supported structure

# Producing all potential predictor combinations
mod_d<- dredge(mod, rank="AICc", extra=c("R^2"))

mod_d[which(mod_d$delta<=7), ]->res_d[[i]]
mod_d[as.character(1:nrow(mod_d)),]->r2_order[[i]]

# Selecting models to do inference
mod_set[[i]] <- subset(mod_d, delta<=7)
mod_set[[i]]

# coefficients average and CI
avg[[i]]<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,1]
se<-summary(model.avg(get.models(mod_d, delta<=7)))$coefmat.subset[-1,2]

ci1<-avg[[i]] - qnorm(0.975)*as.numeric(se)
ci2<-avg[[i]] + qnorm(0.975)*as.numeric(se)

ci[[i]]<-data.frame(ci1,ci2)

# Getting models (we can use these models for other operations)
mods_ch4co2_ratio <- get.models (mod_d, subset=delta<=7) # subset 

# Checking model assumptions (see PDFs in plots/assum/ folder)
length(mods_ch4co2_ratio)->n.mod

for (k in 1:n.mod) 
{
  mods_ch4co2_ratio[[k]]->mod
  pdf(file=paste(assum_folder,"ratio_",k,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model_int",k),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}

# Variance partitioning
sel.pred<-sort(names(mod_d)[2:7]) # Predictor in your model in alphabetical order
round(hier.r2.res(mods_ch4co2_ratio, sel.pred=sel.pred),3)->var.part

# Weighted partitioned variance
res[[i]]<-apply(var.part, 2, function(x) weighted.mean(x, mod_set[[i]]$weight))
res[[i]]

# Saving model results

write.table(mod_set[1], paste(plot_folder,"gpp_mods.txt", sep=""),sep="\t")
write.table(mod_set[2], paste(plot_folder,"er_mods.txt", sep=""),sep="\t")
write.table(mod_set[3], paste(plot_folder,"nep_mods.txt", sep=""),sep="\t")
write.table(mod_set[4], paste(plot_folder,"algal_mods.txt", sep=""),sep="\t")
write.table(mod_set[5], paste(plot_folder,"pco2_mods.txt", sep=""),sep="\t")
write.table(mod_set[6], paste(plot_folder,"pch4_mods.txt", sep=""),sep="\t")
write.table(mod_set[7], paste(plot_folder,"ch4co2_ratio_mods.txt", sep=""),sep="\t")

######################
# Model plots
#######################

var_col<-c("red","red","red","orange","blue","green")

var_names<-c("GPP", "ER","NEP","Algal production", "pCO2" ,"pCH4")

# GPP

ci.d<-ci[[1]][order(avg[[1]]),]
avg.d<-avg[[1]][order(avg[[1]])]

pdf(file=paste(plot_folder,"gpp_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="red")
          panel.xyplot(x, y, pch=15, cex=2,col="red")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# ER

ci.d<-ci[[2]][order(avg[[2]]),]
avg.d<-avg[[2]][order(avg[[2]])]

pdf(file=paste(plot_folder,"er_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="red")
          panel.xyplot(x, y, pch=15, cex=2,col="red")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# NEP

ci.d<-ci[[3]][order(avg[[3]]),]
avg.d<-avg[[3]][order(avg[[3]])]

pdf(file=paste(plot_folder,"nep_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="red")
          panel.xyplot(x, y, pch=15, cex=2,col="red")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# algal

ci.d<-ci[[4]][order(avg[[4]]),]
avg.d<-avg[[4]][order(avg[[4]])]

pdf(file=paste(plot_folder,"algal_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="orange")
          panel.xyplot(x, y, pch=15, cex=2,col="orange")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# pCO2

ci.d<-ci[[5]][order(avg[[5]]),]
avg.d<-avg[[5]][order(avg[[5]])]

pdf(file=paste(plot_folder,"pco2_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="blue")
          panel.xyplot(x, y, pch=15, cex=2,col="blue")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# pCH4

ci.d<-ci[[6]][order(avg[[6]]),]
avg.d<-avg[[6]][order(avg[[6]])]

pdf(file=paste(plot_folder,"pch4_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="green")
          panel.xyplot(x, y, pch=15, cex=2,col="green")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# CH4:CO2 ratio

ci.d<-ci[[7]][order(avg[[7]]),]
avg.d<-avg[[7]][order(avg[[7]])]

pdf(file=paste(plot_folder,"ch4co2_ratio_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(avg.d, xlab="Standardized coefficient", xlim=c(-1.05, 1.05),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci.d[,1], as.numeric(y), ci.d[,2], as.numeric(y), lty=1, col="violet")
          panel.xyplot(x, y, pch=15, cex=2,col="violet")
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

########################################
# Variance partition plots
########################################

# Metabolism

pdf(file=paste(plot_folder,"var_part_metab.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=22)

par(mfrow=c(4,1), cex=2, cex.lab=1, cex.axis=1, mar=c(4,4,3,4)) 

barplot(100*sort(res[[1]]),col=var_col[1],xlab="",
        horiz=T,las=1, xlim=c(0,40),main=var_names[1]) 

barplot(100*sort(res[[2]]),col=var_col[2],xlab="",
        horiz=T,las=1, xlim=c(0,40),main=var_names[2]) 

barplot(100*sort(res[[3]]),col=var_col[3],xlab="",
        horiz=T,las=1, xlim=c(0,40),main=var_names[3]) 

barplot(100*sort(res[[4]]),col=var_col[4],xlab="Explained variance (%)",
        horiz=T,las=1, xlim=c(0,40),main=var_names[4]) 

dev.off()

# GHGs
pdf(file=paste(plot_folder,"var_part_ghg.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=11)

par(mfrow=c(2,1), cex=2, cex.lab=1, cex.axis=1, mar=c(4,4,3,4)) 

barplot(100*sort(res[[5]]),col=var_col[5],xlab="",
        horiz=T,las=1, xlim=c(0,40),main=var_names[5]) 

barplot(100*sort(res[[6]]),col=var_col[6],xlab="Explained variance (%)",
        horiz=T,las=1, xlim=c(0,40),main=var_names[6]) 

dev.off()

# CH4:CO2 ratio

pdf(file=paste(plot_folder,"var_partitioning_ch4co2_ratio.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

par(mfrow=c(1,1), cex=2, cex.lab=1, cex.axis=1, mar=c(4,4,3,4)) 

barplot(100*sort(res[[7]]),col="violet",xlab="Explained variance (%)",
        horiz=T,las=1, xlim=c(0,40)) 

dev.off()

# Save model variables
save.image("saved_models.RData")