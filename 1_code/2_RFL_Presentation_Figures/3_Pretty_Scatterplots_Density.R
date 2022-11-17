library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read.csv("2_incremental/20220310_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE)

names(alldata)

n<-7
n2<-49


dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]


Xlabel=expression(bold(paste("Standing Crop"~"(g/m"^"2"~")")))
Ylabel=expression(bold("Zn Burden (mg/g)"))


#windows(width = 15, height = 25)
#dev.off()     
#dev.new()

par(mfrow = c(2, 1)) 

#par(mfrow = c(1, 3))

########### FIG 1

par(fig = c(0.0, 0.9, 0.665, 1), mar=c(5.5,5.5,1,3)+.1)

BW<-10
ALPHA<-0.6
aty <- seq(0.0, max(density(dfEPIP[,7], na.rm = TRUE, bw = BW)$y), length.out=5)


plot(density(dfEPIP[,7], na.rm = TRUE, bw = BW), 
     xlim = range(alldata[,n],na.rm=TRUE), 
     col=NA,
     lwd=3,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     main=NA,
     ylim=c(0.001,max(density(dfEPIP[,7], na.rm = TRUE, bw = BW)$y)),
     
)     

axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=0.70)

polygon(density(dfEPIP[,7], na.rm = TRUE, bw = BW), 
     xlim = range(alldata[,n],na.rm=TRUE), 
     col=adjustcolor("gold",alpha.f=ALPHA),
#     border="gold",
     lwd=3)  
polygon(density(dfEPIL[,7], na.rm = TRUE, bw = BW), col=adjustcolor(col=colors()[89],alpha.f=ALPHA),
#        border=colors()[89],
      lwd=3)          
polygon(density(dfFILA[,7], na.rm = TRUE, bw = BW), col = adjustcolor("chartreuse3",alpha.f=ALPHA),
#        border="chartreuse3",
      lwd=3) 
box(lwd=2)

########## FIG 2

par(fig = c(0.0, 0.9, 0, 0.85), new=T)

plot(alldata[,n], alldata[,n2], 
     ylim = range(alldata[,n2],na.rm=TRUE),
     xlim = range(alldata[,n],na.rm=TRUE),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)


aty <- seq(0, max(alldata[,n2], na.rm=TRUE), length.out=5)

axis(side = 1, at = seq(0, 200, length.out=5), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5, family="Calibri")
title(ylab=Ylabel, line=2.7, cex.lab=1.5, family="Calibri")


points(dfEPIL[,7], dfEPIL[,n2],pch=23, cex=1.5,col="black", bg=colors()[89],lwd=3)

points(dfEPIP[,7], dfEPIP[,n2],pch=23, cex=1.5,col="black", bg="gold",lwd=3)

points(dfFILA[,7], dfFILA[,n2],pch=23, cex=1.5,col="black", bg="chartreuse3",lwd=3)

box(lwd=2)
#title(ylab=Ylabel, line=2.5, cex.lab=1.2, family="Calibri Light")



