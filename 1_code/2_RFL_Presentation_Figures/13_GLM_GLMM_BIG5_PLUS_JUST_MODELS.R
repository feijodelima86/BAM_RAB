library("readr")
library("corrplot")
library("vegan")
library("lme4")
library("sjPlot")
library("ggthemes") 
library("ggplot2"); theme_set(theme_bw())
library("tidyverse")

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
alldata <- read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

EPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

#COMPARTMENTS <- COMPARTMENTS[-217,]

EPIL$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[EPIL$SITE]

BIG5.LABELS<-as.data.frame(na.omit(EPIL[,c(c(3,4,6,89),c(13,19,22,25,35,40,49)+39)]))
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
colnames(BIG5.LABELS[5:11])<-c("As","Cd","Cu","Fe","Pb","Se","Zn")
BIG5.GLM<-BIG5.LABELS
BIG5.DF[,c(7)]
####As####

mix.int <- glm(log(As.1*1000) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Cd####

mix.int <- glm(log(Cd.1*1000000) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Cu####

mix.int <- glm(log(Cu.1+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Fe####

mix.int <- glm(log(Fer.1.1+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Pb####

mix.int <- glm(log(Pb.1+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Se####

mix.int <- glm(log(Se.1+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

####Zn####

mix.int <- glm(log(Znr.1+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

