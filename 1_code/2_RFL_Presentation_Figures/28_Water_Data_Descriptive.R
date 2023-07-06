library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))
library("scales")
library("performance")
library("splines")
library("effects")
library("mgcv")

alldata <- read_csv("2_incremental/wateralgae13_WORKING_MANUAL_W.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIL_F"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG <-aggregate(x = COMPARTMENTS[,colnames(COMPARTMENTS) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
                             by = list(COMPARTMENTS$SAMPLING_DATE,COMPARTMENTS$SITE,COMPARTMENTS$SAMPLE_DESCRIPTOR),
                             FUN = mean,
                             na.rm = T)

write.csv(COMPARTMENTS_AVG, "2_incremental/COMPARTMENTS_AVG_3.csv")

COMPARTMENTS <- read.csv("2_incremental/COMPARTMENTS_AVG_3.csv")

colnames(COMPARTMENTS)[c(2:4)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

COMPARTMENTS<-COMPARTMENTS[,c(2,3,4,8,10:63)]

names(alldata)

COMPARTMENTS$SITE<-as.factor(COMPARTMENTS$SITE)

COMPARTMENTS$SAMPLING_DATE<-as.Date(COMPARTMENTS$SAMPLING_DATE)

plot(COMPARTMENTS$SITE_DISTANCE, COMPARTMENTS$TD.As, col=as.factor(COMPARTMENTS$SAMPLING_DATE))








