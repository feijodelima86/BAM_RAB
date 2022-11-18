library("corrplot")                               # Load corrplot
library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

#Big 5

COR.DF<-as.data.frame(na.omit(alldata[,c(13,19,22,25,35,49)]))

names(COR.DF)<- c("As","Cd","Cu","Fe","Pb","Zn")

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


#EPIL.COR.DF3<-as.data.frame(na.omit(EPIL.COR.DF3[,c(13,16,17,18,19,20,21,22,23,26,31,35,36,40,41,44,46,47,48)]))

#EPIP.COR.DF3<-as.data.frame(na.omit(EPIP.COR.DF3[,c(13,16,17,18,19,20,21,22,23,26,31,35,36,40,41,44,46,47,48)]))

#FILA.COR.DF3<-as.data.frame(na.omit(FILA.COR.DF3[,c(13,16,17,18,19,20,21,22,23,26,31,35,36,40,41,44,46,47,48)]))

EPIL.COR.DF3<-as.data.frame(na.omit(EPIL.COR.DF3[,c(52,58,61,64,74,79,88)]))

EPIP.COR.DF3<-as.data.frame(na.omit(EPIP.COR.DF3[,c(52,58,61,64,74,79,88)]))

FILA.COR.DF3<-as.data.frame(na.omit(FILA.COR.DF3[,c(52,58,61,64,74,79,88)]))



par(mfrow = c(1,3), mar = c(2, 4, 4, 2))

EPIL.TEST<-cor.mtest(EPIL.COR.DF3)

corrplot(cor(EPIL.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=2, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epilithon", diag=FALSE, order = "hclust", type="upper") 

corrplot(cor(EPIP.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=2, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Epiphytes", diag=FALSE, order = "hclust", type="upper") 

corrplot(cor(FILA.COR.DF3), method = "circle", lwd=2, tl.cex = 2, cl.cex=2, cex.main=3, sig.level = 0.05, insig = "blank", 
         mar=c(0,0,0,0), title="\n\n Filamentous", diag=FALSE, order = "hclust", type="upper") 

cor(EPIL.COR.DF3)
cor(EPIP.COR.DF3)
cor(FILA.COR.DF3)

###########GEARING UP FOR MULTIVARIATE ANALYSIS###################

###Testing compartment datasets for variance inflation factor.

library(usdm)

max.VIF<-10

EPIL.RED<-data.frame(EPIL.COR.DF3)
VIF.EPIL<-vif(EPIL.RED)

repeat {   
  if (max(VIF.EPIL[,2])>max.VIF) {
    VIF.EPIL<-vif(EPIL.RED)
    EPIL.RED <- within(EPIL.RED, rm(list=VIF.EPIL[which.max(VIF.EPIL[,2]),1]))
  } else {
    break
  }
}

EPIP.RED<-data.frame(EPIP.COR.DF3)
VIF.EPIP<-vif(EPIP.RED)

repeat {   
  if (max(VIF.EPIP[,2])>max.VIF) {
    VIF.EPIP<-vif(EPIP.RED)
    EPIP.RED <- within(EPIP.RED, rm(list=VIF.EPIP[which.max(VIF.EPIP[,2]),1]))
  } else {
    break
  }
}

FILA.RED<-data.frame(FILA.COR.DF3)
VIF.FILA<-vif(FILA.RED)

repeat {   
  if (max(VIF.FILA[,2])>max.VIF) {
    VIF.FILA<-vif(FILA.RED)
    FILA.RED <- within(FILA.RED, rm(list=VIF.FILA[which.max(VIF.FILA[,2]),1]))
  } else {
    break
  }
}

ALL.COR.DF2<-as.data.frame(na.omit(alldata[,c(13:49)]))

#Selected

ALL.COR.DF3<-as.data.frame(na.omit(alldata[,c(13,16,17,18,19,20,21,22,23,26,31,35,36,40,41,44,46,47,48)]))

ALL.RED<-data.frame(ALL.COR.DF3)
VIF.ALL<-vif(ALL.RED)

repeat {   
  if (max(VIF.ALL[,2])>max.VIF) {
    VIF.ALL<-vif(ALL.RED)
    ALL.RED <- within(ALL.RED, rm(list=VIF.ALL[which.max(VIF.ALL[,2]),1]))
  } else {
    break
  }
}





