library("corrplot")                               # Load corrplot

library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- data.frame(read.csv("2_incremental/20220313_STANDING_CROP.csv"))

names(alldata)

COR.DF<-as.data.frame(na.omit(alldata[,c(13,19,22,25,35,49)]))

names(COR.DF)<- c("As","Cd","Cu","Fe","Pb","Zn")

cor(COR.DF)

corrplot(cor(COR.DF), method = "circle", lwd=2) 


COR.DF2<-as.data.frame(na.omit(alldata[,c(13:49)]))

cor(COR.DF2)

corrplot(cor(COR.DF2), method = "circle", lwd=2) 