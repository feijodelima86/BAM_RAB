library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read.csv("2_incremental/20220310_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE)

names(alldata)

se <- function(x, ...) sqrt(var(x, ...)/length(x))

n<-7
n2<-13


#Xlabel=expression(bold(paste("Biomass g/m"^"2")))

Ylabel=expression(bold("Burden (mg/g)"))
  

dev.new()
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

axis(side = 1, at = seq(0, 200, length.out=5), las=1, font.axis=2, cex.axis=1.5)
axis(side = 2, at = seq(0, .03, length.out=5), las=1, font.axis=2, cex.axis=1.5)

title(xlab=Xlabel, line=2.5, cex.lab=1.75, family="Calibri Light")

range(alldata[,n],na.rm=TRUE)


dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

points(dfEPIL[,7], dfEPIL[,n2],pch=23, cex=1.5,col="black", bg=colors()[89],lwd=3)

points(dfEPIP[,7], dfEPIP[,n2],pch=23, cex=1.5,col="black", bg="gold",lwd=3)

points(dfFILA[,7], dfFILA[,n2],pch=23, cex=1.5,col="black", bg="chartreuse3",lwd=3)


#title(ylab=Ylabel, line=2.5, cex.lab=1.2, family="Calibri Light")

box(lwd=2)


