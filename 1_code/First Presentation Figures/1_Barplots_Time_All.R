library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read.csv("2_incremental/20220313_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE)

alldata$SAMPLING_DATE<-factor(alldata$SAMPLING_DATE,levels = c("6/22/2021", "7/7/2021", "7/20/2021", "8/3/2021", "8/17/2021","9/9/2021", "9/22/2021", "10/13/2021"), ordered=TRUE)

# Variable date to factor

dfALL <- alldata

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

ssite<-"BN"

se <- function(x, ...) sqrt(var(x, ...)/length(x))

Ylabel=expression(bold("Zn Burden (mg/g)"))

names(alldata)

n1<-49

mult=1

ALL.SUM<-data.frame(aggregate(as.numeric(dfALL[,n1]) ~ SAMPLING_DATE+SAMPLE_DESCRIPTOR, dfALL, mean))
ALL.SUM[,4]<-aggregate(as.numeric(dfALL[,n1]) ~ SAMPLING_DATE+SAMPLE_DESCRIPTOR, dfALL, se)[,3]
#ALL.SUM <- subset(ALL.SUM, SITE == ssite)
names(ALL.SUM)<-c("DATE","SAMPLE_DESCRIPTOR","MEAN.ALL","STDER.ALL")
ALL.SUM[,5]<--2
ALL.SUM[,c(3,4)]<-ALL.SUM[,c(3,4)]*mult 

ALL.SUM[is.na(ALL.SUM)] <- 0

ALL.SUM[ALL.SUM < 0] <- 0


aty <- seq(0, max(ALL.SUM$MEAN.ALL, na.rm=TRUE), length.out=5)

# Uniform color

names(ALL.SUM)

#dev.new()

tabbedMeans <- tapply(ALL.SUM$MEAN.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR,
                                             ALL.SUM$DATE),
                      function(x) c(x = x))
tabbedSE <- tapply(ALL.SUM$STDER.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR,
                                           ALL.SUM$DATE),
                   function(x) c(x = x))

barCenters <- barplot(MEAN.ALL ~ SAMPLE_DESCRIPTOR+DATE, data = ALL.SUM, 
                      beside = TRUE, 
                      ylim = c(0, max(ALL.SUM$MEAN.ALL+ALL.SUM$STDER.ALL, na.rm=T)*1.1),
                      xlab="",
                      ylab=Ylabel,
                      #                     col=colors()[89],
                      col=c(colors()[89],"gold" , "chartreuse"),
                      font=2,
                      lwd=2
)

segments(barCenters, tabbedMeans, barCenters,
         tabbedMeans + tabbedSE, lwd = 2)

arrows(barCenters, tabbedMeans, barCenters,
       tabbedMeans + tabbedSE, lwd = 2, angle = 90,
       code = 2, length = 0.05)

box(lwd=3)

















