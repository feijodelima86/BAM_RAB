library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read.csv("2_incremental/20220313_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE)

#COR.DF<-as.data.frame(na.omit(alldata[,c(13,19,22,25,35,49)]))

names(alldata)

n2<-22
n<-49

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]


Xlabel=expression(bold(paste("Cu Burden (mg/g)")))
Ylabel=expression(bold("Zn Burden (mg/g)"))


dev.new()

par(mar=c(7,7,7,7))

#par(mfrow = c(1, 3))

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

axis(side = 1, at = seq(0, 200, length.out=5), las=1, font.axis=2, cex.axis=1.5)
axis(side = 2, at = aty, labels=format(aty, scientific=TRUE,digits = 1), las=1, font.axis=2, cex.axis=1.5)

title(xlab=Xlabel, line=3, cex.lab=1.75, family="Calibri")
title(ylab=Ylabel, line=5, cex.lab=1.75, family="Calibri")

range(alldata[,n],na.rm=TRUE)



points(dfEPIL[,n], dfEPIL[,n2],pch=23, cex=1.5,col="black", bg=colors()[89],lwd=3)

points(dfEPIP[,n], dfEPIP[,n2],pch=23, cex=1.5,col="black", bg="gold",lwd=3)

points(dfFILA[,n], dfFILA[,n2],pch=23, cex=1.5,col="black", bg="chartreuse3",lwd=3)


#title(ylab=Ylabel, line=2.5, cex.lab=1.2, family="Calibri Light")

box(lwd=2)


