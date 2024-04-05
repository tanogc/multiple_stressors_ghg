# Setting working folder
setwd("C:/Users/tanog/Dropbox/Trabajos en curso/streameco_ghg/script/def")

# Setting results plots
plot_folder<-paste(getwd(),"/plots/",sep="")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")

# Loading packages
library(usdm)
library(corrplot)
library(ade4)
library(viridis)

# Loading data
load("saved_models.RData") # Previously, you need to run "2_multimodel_inference.txt"
dat_raw<-read.table("dat_completed.txt", sep="\t", h=T) # saving raw values

# back-std
sd(log(dat_raw$pCO2))->sd.co2
mean(log(dat_raw$pCO2))->m.co2
sd(log(dat_raw$pCH4))->sd.ch4
mean(log(dat_raw$pCH4))->m.ch4

sd(log(dat_raw$din))->sd.din
mean(log(dat_raw$din))->m.din
sd(dat_raw$tmean)->sd.tmean
mean(dat_raw$tmean)->m.tmean
sd(log(dat_raw$discharge))->sd.disc
mean(log(dat_raw$discharge))->m.disc
sd(log(dat_raw$do.def))->sd.do
mean(log(dat_raw$do.def))->m.do

# Creating sequences
n.seq<-seq(min(dat$din), max(dat$din), length.out=1000)
n.seqr<-seq(min(dat_raw$din), max(dat_raw$din), length.out=1000)

col_ghg=viridis(4)

# Calculating predictions

# Calculating weights for weighted-mean predictions
mod.weights <- Weights(as.numeric(lapply(mods_co2, AICc)))

# Only DIN
mod.pred <- data.frame(lapply(mods_co2, function(x) { 
  
  predict(x,data.frame(din=n.seq,
                       tmean=mean(dat$tmean),
                       discharge=mean(dat$discharge),
                       rip_open=mean(dat$rip_open)))
  
}))

co2_pred <- exp((apply(mod.pred, 1, function(x) sum(x * mod.weights))*sd.co2)+m.co2)

a=0

for (i in c(0.25, 1.5, 7)) {

  a=a+1

# CH4 - Only DIN

# Calculating weights for weighted-mean predictions
mod.weights <- Weights(as.numeric(lapply(mods_ch4, AICc)))
  
mod.pred <- data.frame(lapply(mods_ch4, function(x) { 
  
  predict(x,data.frame(din=n.seq,
                       tmean=mean(dat$tmean),
                       do.def=(log(c(i))-m.do)/sd.do,
                       discharge=mean(dat$discharge),
                       rip_open=mean(dat$rip_open)))
  
}))

ch4_pred<-exp((apply(mod.pred, 1, function(x) sum(x * mod.weights))*sd.ch4)+m.ch4)

# Summing CO2-equivalent
ghg<-co2_pred+(28*ch4_pred)

# Plotting line
lines(n.seqr, ghg, col=col_ghg[a], lwd=3)

}

# Calculating weights for weighted-mean predictions
mod.weights <- Weights(as.numeric(lapply(mods_co2, AICc)))

# Only DIN
mod.pred <- data.frame(lapply(mods_co2, function(x) { 
  
  predict(x,data.frame(din=(log(c(0.25, 0.5,1,2.5,5))-m.din)/sd.din,
                       tmean=mean(dat$tmean),
                       discharge=mean(dat$discharge),
                       rip_open=mean(dat$rip_open)))
  
}))

co2_pred <- exp((apply(mod.pred, 1, function(x) sum(x * mod.weights))*sd.co2)+m.co2)


ch4_pred<-matrix(NA, nrow=3, ncol=5)

a=0
for (i in c(0.25, 1.5, 7)) {
  
  a=a+1
  
  # CH4 - Only DIN
  # Calculating weights for weighted-mean predictions
  mod.weights <- Weights(as.numeric(lapply(mods_ch4, AICc)))
  
  mod.pred <- data.frame(lapply(mods_ch4, function(x) { 
    
    predict(x,data.frame(din=(log(c(0.25, 0.5,1,2.5,5))-m.din)/sd.din,
                         tmean=mean(dat$tmean),
                         do.def=(log(c(i))-m.do)/sd.do,
                         discharge=mean(dat$discharge),
                         rip_open=mean(dat$rip_open)))
    
  }))
  
  ch4_pred[a,]<-exp((apply(mod.pred, 1, function(x) sum(x * mod.weights))*sd.ch4)+m.ch4)
  
}

# Calculating predicted CO2-eq concentrations

ch4_pred<-ch4_pred*28
value=c(rep(co2_pred,4), ch4_pred[1,],ch4_pred[2,],ch4_pred[3,])

value=c(rep(co2_pred,4), ch4_pred[1,],ch4_pred[2,],ch4_pred[3,])

v1=c(co2_pred,NA,co2_pred,NA,co2_pred)
v2=c(ch4_pred[1,],NA,ch4_pred[2,],NA,ch4_pred[3,])

valmat<-rbind(v1,v2)

#################
# Scenarios plots
#################

stacked_data1 <- data.frame(
  din = rep(c(0.25, 0.5,1,2.5,5), 2), # DIN
  co2_eq = c(co2_pred, ch4_pred[1,]),
  gas = rep(c("CO2","CH4"), each = 5)
)

stacked_data2 <- data.frame(
  din = rep(c(0.25, 0.5,1,2.5,5), 2), # DIN
  co2_eq = c(co2_pred, ch4_pred[2,]),
  gas = rep(c("CO2","CH4"), each = 5)
)

stacked_data3 <- data.frame(
  din = rep(c(0.25, 0.5,1,2.5,5), 2), # DIN
  co2_eq = c(co2_pred, ch4_pred[3,]),
  gas = rep(c("CO2","CH4"), each = 5)
)

stacked_labels <- data.frame(
  labels = c(
    "CO2", 
    "CH4"
  ),
  x = c(2, 2),
  y = c(1000, 1750),
  color = c("white","black")
)

###############################
# Basic stacked line plot
###############################

plt1 <- ggplot(stacked_data1) +
  # color = "white" indicates the color of the lines between the areas
  geom_area(aes(din, co2_eq, fill = gas), color = "white") +
  geom_hline(yintercept = 413, lty=2, lwd=1.25, col="grey") +
  scale_fill_manual(values = c("green","blue")) +
  theme(legend.position = "None") # no legend 

plt1

plt2 <- ggplot(stacked_data2) +
  # color = "white" indicates the color of the lines between the areas
  geom_area(aes(din, co2_eq, fill = gas), color = "white") +
  geom_hline(yintercept = 413, lty=2, lwd=1.25, col="grey") +
  scale_fill_manual(values = c("green","blue")) +
  theme(legend.position = "None") # no legend

plt2

plt3 <- ggplot(stacked_data3) +
  # color = "white" indicates the color of the lines between the areas
  geom_area(aes(din, co2_eq, fill = gas), color = "white") +
  geom_hline(yintercept = 413, lty=2, lwd=1.25, col="grey") +
  scale_fill_manual(values = c("green","blue")) +
  theme(legend.position = "None") # no legend

plt3

###############################
# Basic stacked line plot
###############################

plt1 <- plt1 + 
  scale_x_continuous(
    # Note: Data goes from 2008 to 2020. Extra space is added on the right
    # so there's room for the grid line labels ;)
    limits = c(0, 5),
    expand = c(0, 0),
    breaks = seq(1, 5, by = 1), 
    labels = c("1", "2","3","4","5")
  ) +
  scale_y_continuous(
    limits = c(0, 2500),
    breaks = seq(0, 2500, by = 500), 
    expand = c(0, 0)
  ) + 
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Remove all grid lines
    panel.grid = element_blank(),
    # But add grid lines for the vertical axis, customizing color and size 
    panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length.y = unit(0, "mm"), 
    axis.ticks.length.x = unit(2, "mm"),
    # Remove the title for both axes
    axis.title = element_blank(),
    # Only bottom line of the vertical axis is painted in black
    axis.line.x.bottom = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Helvetica", size = 16)
  )

plt1

plt2 <- plt2 + 
  scale_x_continuous(
    # Note: Data goes from 2008 to 2020. Extra space is added on the right
    # so there's room for the grid line labels ;)
    limits = c(0, 5),
    expand = c(0, 0),
    breaks = seq(1, 5, by = 1), 
    labels = c("1", "2","3","4","5")
  ) +
  scale_y_continuous(
    limits = c(0, 2500),
    breaks = seq(0, 2500, by = 500), 
    expand = c(0, 0)
  ) + 
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Remove all grid lines
    panel.grid = element_blank(),
    # But add grid lines for the vertical axis, customizing color and size 
    panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length.y = unit(0, "mm"), 
    axis.ticks.length.x = unit(2, "mm"),
    # Remove the title for both axes
    axis.title = element_blank(),
    # Only bottom line of the vertical axis is painted in black
    axis.line.x.bottom = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Helvetica", size = 16)
  )

plt2


plt3 <- plt3 + 
  scale_x_continuous(
    # Note: Data goes from 2008 to 2020. Extra space is added on the right
    # so there's room for the grid line labels ;)
    limits = c(0, 5),
    expand = c(0, 0),
    breaks = seq(1, 5, by = 1), 
    labels = c("1", "2","3","4","5")
  ) +
  scale_y_continuous(
    limits = c(0, 2500),
    breaks = seq(0, 2500, by = 500), 
    expand = c(0, 0)
  ) + 
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Remove all grid lines
    panel.grid = element_blank(),
    # But add grid lines for the vertical axis, customizing color and size 
    panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length.y = unit(0, "mm"), 
    axis.ticks.length.x = unit(2, "mm"),
    # Remove the title for both axes
    axis.title = element_blank(),
    # Only bottom line of the vertical axis is painted in black
    axis.line.x.bottom = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Helvetica", size = 16)
  )

plt3

#########################
# Complex stacked line plot
#########################

plt3 <- plt3 + 
  geom_text(
    aes(x, y, label = labels, color = color),
    data = stacked_labels,
    hjust = 0,
    vjust = 0.5,
    family = "Helvetica",
    size = 6,
    inherit.aes = FALSE
  ) + 
  scale_color_identity()

plt3 <- plt3 +
  geom_text(
    data = data.frame(x = 0, y = seq(0, 2500, by = 500)),
    aes(x, y, label = y),
    hjust = 1,
    vjust = 0,
    nudge_y = 2000 * 0.01, # Again, the pad is equal to 1% of the vertical range.
    family = "Helvetica",
    size = 6,
    inherit.aes = FALSE
  )

# Note again we use `element_markdown()` to render Markdown content
plt1 <- plt1 + 
  labs(
    title = "DO deficit=0.25 mg/l",
  ) + 
  theme(
    plot.title = element_markdown(
      family = "Helvetica", 
      size = 14
    )
  )

plt2 <- plt2 + 
  labs(
    title = "DO deficit=3 mg/l",
  ) + 
  theme(
    plot.title = element_markdown(
      family = "Helvetica", 
      size = 14
    )
  )

plt3 <- plt3 + 
  labs(
    title = "DO deficit=7 mg/l",
  ) + 
  theme(
    plot.title = element_markdown(
      family = "Helvetica", 
      size = 14
    )
  )

plt3

##################
# Combining plots
##################

plt1 <- plt1 + theme(plot.margin = margin(0, 0.05, 0.05, 0, "npc"))
plt2 <- plt2 + theme(plot.margin = margin(0, 0.05, 0.05, 0, "npc"))
plt3 <- plt3 + theme(plot.margin = margin(0, 0.05, 0.05, 0, "npc"))
plt <- plt1 | plt2| plt3

title_theme <- theme(
  plot.title = element_text(
    family = "Helvetica", 
    face = "bold",
    size = 22,
    margin = margin(0.8, 0, 0, 0, "npc")
  ),
  plot.subtitle = element_text(
    family = "Helvetica",
    size = 20,
    margin = margin(0.4, 0, 0, 0, "npc")
  )
)

plt <- plt + plot_annotation(
  title = labs(title = parse(text = "X[1]")),
  subtitle = "Different*",
  theme = title_theme
) +
  theme(
    plot.margin = margin(0.075, 0, 0.1, 0, "npc"),
  ) 

pdf(file=paste(plot_folder,"scenarios_both2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=10,height=5)

plt

dev.off()

