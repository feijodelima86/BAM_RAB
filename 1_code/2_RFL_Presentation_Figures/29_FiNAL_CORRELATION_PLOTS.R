library(tidyverse)
library(readr)
library(plyr)
library(dplyr)
library("corrplot")                               # Load corrplot

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

#Selected

EPIL.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

EPIL.COR.DF3<-as.data.frame(na.omit(EPIL.COR.DF3[,c(13,19,22,25,31,35,40,49)]))

names(EPIL.COR.DF3)<- c("As","Cd","Cu","Fe","Mo","Pb","Se","Zn")

EPIP.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

EPIP.COR.DF3<-as.data.frame(na.omit(EPIP.COR.DF3[,c(13,19,22,25,31,35,40,49)]))

names(EPIP.COR.DF3)<- c("As","Cd","Cu","Fe","Mo","Pb","Se","Zn")

FILA.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

FILA.COR.DF3<-as.data.frame(na.omit(FILA.COR.DF3[,c(13,19,22,25,31,35,40,49)]))

names(FILA.COR.DF3)<- c("As","Cd","Cu","Fe","Mo","Pb","Se","Zn")

dev.new()

par(mfrow = c(1,3), mar = c(2, 4, 4, 2))

EPIL.TEST<-cor.mtest(EPIL.COR.DF3)

corrplot(cor(EPIL.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1, cex.main=2.5, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epilithon", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

corrplot(cor(EPIP.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1, cex.main=2.5, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epiphytes", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

corrplot(cor(FILA.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1, cex.main=2.5, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Filamentous", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

cor(EPIL.COR.DF3)
cor(EPIP.COR.DF3)
cor(FILA.COR.DF3)
