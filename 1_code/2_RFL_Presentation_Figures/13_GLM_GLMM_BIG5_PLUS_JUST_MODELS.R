library("readr")
library("corrplot")
library("vegan")
library("lme4")
library("sjPlot")
library("ggthemes") 
library("ggplot2"); theme_set(theme_bw())
library("tidyverse")

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
#alldata <- read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

EPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

EPIL <- alldata[which(EPIL$SITE != "WS"),]

?which

#COMPARTMENTS <- COMPARTMENTS[-217,]

EPIL$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[EPIL$SITE]

#+39

BIG5.LABELS<-as.data.frame(na.omit(EPIL[,c(c(3,4,6,89),c(13,19,22,25,35,40,49+39))]))
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
names(BIG5.LABELS)[c(5:11)]<- c("As","Cd","Cu","Fe","Pb","Se","Zn")
BIG5.GLM<-BIG5.LABELS


####As####

mix.int <- glm(As ~ SAMPLING_DATE + SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Cd####

mix.int <- glm(Cd.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Cu####

mix.int <- glm(Cu.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Fe####

mix.int <- glm(Fer.1.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Pb####

mix.int <- glm(Pb.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Se####

mix.int <- glm(Se.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

####Zn####

mix.int <- glm(Znr.1 ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

