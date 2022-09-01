library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- read.csv("2_incremental/20220313_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-factor(alldata$SAMPLING_DATE,levels = c("6/22/2021", "7/7/2021", "7/20/2021", "8/3/2021", "8/17/2021","9/9/2021", "9/22/2021", "10/13/2021"), ordered=TRUE)

# Variable date to factor

dfALL <- alldata

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

se <- function(x, ...) sqrt(var(x, ...)/length(x))

Ylabel=expression(bold("As Burden (mg/g)"))

names(alldata)

n1<-49

mult=1

ALL.SUM<-data.frame(aggregate(as.numeric(dfALL[,n1]) ~ SAMPLE_DESCRIPTOR*SAMPLING_DATE, dfALL, mean))
ALL.SUM[,4]<-aggregate(as.numeric(dfALL[,n1]) ~ SAMPLE_DESCRIPTOR*SAMPLING_DATE, dfALL, se)[,3]

#ALL.SUM <- subset(ALL.SUM, SITE == ssite)

names(ALL.SUM)<-c("SAMPLE_DESCRIPTOR","SAMPLING_DATE","MEAN.ALL","STDER.ALL")
ALL.SUM[,5]<--2
ALL.SUM[,c(3,4)]<-ALL.SUM[,c(3,4)]*mult 

ALL.SUM[is.na(ALL.SUM)] <- 0

ALL.SUM[ALL.SUM < 0] <- 0

# Uniform color

dev.new()

names(ALL.SUM)

tabbedMeans <- tapply(ALL.SUM$MEAN.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR
),
function(x) c(x = x))
tabbedSE <- tapply(ALL.SUM$STDER.ALL, list(ALL.SUM$SAMPLE_DESCRIPTOR),
                   function(x) c(x = x))

a<-((ALL.SUM[which(ALL.SUM$SAMPLE_DESCRIPTOR == "FILA"),]))
names(a)<-c("SAMPLE_DESCRIPTOR","SAMPLING_DATE","FILA.ALL","STDER.ALL")
a<-((ALL.SUM[which(ALL.SUM$SAMPLE_DESCRIPTOR == "FILA"),]))
names(a)<-c("SAMPLE_DESCRIPTOR","SAMPLING_DATE","FILA.ALL","STDER.ALL")


dev.new()
par(mar=c(9, 9, 9, 9))
barCenters <- barplot(MEAN.ALL ~ SAMPLING_DATE, data = ALL.SUM, 
                      beside = TRUE, 
                      ylim = c(0.001, max(ALL.SUM$MEAN.ALL+ALL.SUM$STDER.ALL, na.rm=T)*1.3),
                      xlab = NA,
                      yaxt = "n",
                      ylab = NA,
                      lwd=3,
                      col=c(colors()[89],"gold" , "chartreuse"),
                      log="y",
                      font=2
)

segments(barCenters, tabbedMeans, barCenters,
         tabbedMeans + tabbedSE, lwd = 2)

arrows(barCenters, tabbedMeans, barCenters,
       tabbedMeans + tabbedSE, lwd = 2, angle = 90,
       code = 2, length = 0.125)

box(lwd=3)

lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

aty <- lseq(0.001, max(ALL.SUM$MEAN.ALL, na.rm=TRUE), length.out=5)

#aty <- seq(0.00, max(ALL.SUM$MEAN.ALL, na.rm=TRUE), length.out=5)


axis(side = 2, at = aty, labels=format(aty, scientific=T,digits = 2), las=1, font.axis=2, cex.axis=1)

#title(xlab=Xlabel, line=3, cex.lab=1.75, family="Calibri")

title(ylab=Ylabel, line=4, cex.lab=1.5, family="Calibri")

one.way <- aov(alldata[,n] ~ SAMPLE_DESCRIPTOR, data = alldata)

summary(one.way)

tukey.two.way<-TukeyHSD(one.way)

tukey.two.way










