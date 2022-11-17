library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read_csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-factor(alldata$SAMPLING_DATE,levels = c("6/22/2021", "7/7/2021", "7/20/2021", "8/3/2021", "8/17/2021","9/9/2021", "9/22/2021", "10/13/2021"), ordered=TRUE)

# Variable date to factor

dfALL <- data.frame(alldata[which(alldata$SAMPLE_DESCRIPTOR == c("EPIL","EPIP","FILA")),])

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

se <- function(x, ...) sqrt(var(x, ...)/length(x))

names(alldata)

Ylabel=expression(bold("Cu Contents (mg/g)"))

n1<-22

mult=1

ALL.SUM<-data.frame(aggregate(as.numeric(dfALL[,n1]) ~ SITE+SAMPLE_DESCRIPTOR, dfALL, mean))
ALL.SUM[,4]<-aggregate(as.numeric(dfALL[,n1]) ~ SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,3]
#ALL.SUM <- subset(ALL.SUM, SITE == ssite)
names(ALL.SUM)<-c("SITE","SAMPLE_DESCRIPTOR","MEAN.ALL","STDER.ALL")
ALL.SUM[,5]<--2
ALL.SUM[,c(3,4)]<-ALL.SUM[,c(3,4)]*mult 

ALL.SUM[is.na(ALL.SUM)] <- 0

ALL.SUM[ALL.SUM < 0] <- 0

ALL.SUM$SITE <- factor(ALL.SUM$SITE, c("WS","DL","GR","GC","BG","BN"))


aty <- seq(0, max(ALL.SUM$MEAN.ALL, na.rm=TRUE), length.out=5)

# Uniform color

names(ALL.SUM)


#dev.new()

tabbedMeans <- tapply(ALL.SUM$MEAN.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR,
                                             ALL.SUM$SITE),
                      function(x) c(x = x))
tabbedSE <- tapply(ALL.SUM$STDER.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR,
                                           ALL.SUM$SITE),
                   function(x) c(x = x))


barCenters <- barplot(MEAN.ALL ~ SAMPLE_DESCRIPTOR+SITE, data = ALL.SUM, 
                      beside = TRUE, 
                      ylim = c(0, max(ALL.SUM$MEAN.ALL+ALL.SUM$STDER.ALL, na.rm=T)*1.1),
                      xlab="",
                      ylab=Ylabel,
                      #                     col=colors()[89],
                      col=c(colors()[89],"gold" , "chartreuse"),
                      font=2,
                      lwd=2,
                      space=c(0, 4)
)

try(segments(barCenters, tabbedMeans, barCenters,tabbedMeans + tabbedSE, lwd = 2), silent=TRUE)

try(arrows(barCenters, tabbedMeans, barCenters, tabbedMeans + tabbedSE, lwd = 2, angle = 90, code = 2, length = 0.05), silent=TRUE)

try(box(lwd=3), silent=TRUE)














