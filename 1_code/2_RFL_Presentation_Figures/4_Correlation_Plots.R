library(tidyverse)
library(readr)
library(plyr)
library(dplyr)
library("corrplot")                               # Load corrplot

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

#Big 5



COR.DF<-as.data.frame(na.omit(alldata[,c(13,19,22,25,35,28,49,34,36,32,26,28,43,11,33)+39]))

names(COR.DF)<- c("As","Cd","Cu","Fe","Pb","Zn","P","S","Na","Mg","K","Ca","Si","Al","Ni")

cor(COR.DF)

corrplot(cor(COR.DF), method = "circle", lwd=2) 

#Alldata

EPIL.COR.DF2 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

EPIP.COR.DF2 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

FILA.COR.DF2 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

EPIL.COR.DF2<-as.data.frame(na.omit(EPIL.COR.DF2[,c(13:49)]))

EPIP.COR.DF2<-as.data.frame(na.omit(EPIP.COR.DF2[,c(13:49)]))

FILA.COR.DF2<-as.data.frame(na.omit(FILA.COR.DF2[,c(13:49)]))

#dev.new()

corrplot(cor(FILA.COR.DF2), method = "circle", lwd=2) 

#Selected

EPIL.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

EPIP.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

FILA.COR.DF3 <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]


EPIL.COR.DF3<-as.data.frame(na.omit(EPIL.COR.DF3[,c(13,19,22,25,35,28,49,11,33,34,36,32,26,18,43)]))

EPIP.COR.DF3<-as.data.frame(na.omit(EPIP.COR.DF3[,c(13,19,22,25,35,28,49,11,33,34,36,32,26,18,43)]))

FILA.COR.DF3<-as.data.frame(na.omit(FILA.COR.DF3[,c(13,19,22,25,35,28,49,11,33,34,36,32,26,18,43)]))



dev.new()

par(mfrow = c(1,3), mar = c(2, 4, 4, 2))

EPIL.TEST<-cor.mtest(EPIL.COR.DF3)

corrplot(cor(EPIL.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1.5, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epilithon", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

corrplot(cor(EPIP.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1.5, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epiphytes", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

corrplot(cor(FILA.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=1.5, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Filamentous", diag=FALSE, order = "original", type="upper", col.lim=c(-1,1)) 

cor(EPIL.COR.DF3)
cor(EPIP.COR.DF3)
cor(FILA.COR.DF3)
