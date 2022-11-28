library("readr")
library("corrplot")
library("vegan")
library("lme4")
library("sjPlot")
library("ggthemes") 
library("ggplot2"); theme_set(theme_bw())
library("tidyverse")

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]
COMPARTMENTS <- COMPARTMENTS[-217,]
COMPARTMENTS <- COMPARTMENTS[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]



nrow(COMPARTMENTS)
names(COMPARTMENTS)

COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]

# make a new character column 'cw$diet_name', similar to 'qw$rownames'

 # set legend order with levels



####Least VIF####

#LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,6),c(17,18,20,31,35,36,40,41,44)+39)]))
#LABELS<-LABELS[order(LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
#PCA.DF<-as.data.frame(LABELS[,c(4:ncol(LABELS))])
#names(LABELS)
#names(PCA.DF)<- c("Be","Ca","Co","Mo","Pb","S","Se","Si","Sn")

#cor(PCA.DF)
#corrplot(cor(PCA.DF), method = "circle", lwd=2) 


####Big 5+Fe+Se####

#+39

names(COMPARTMENTS)

BIG5.LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,6,89),c(13,19,22,25,35,40,49))]))
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
BIG5.LABELS<-BIG5.LABELS[order(BIG5.LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
BIG5.DF<-as.data.frame(BIG5.LABELS[,c(5:ncol(BIG5.LABELS))])
names(BIG5.DF)<- c("As","Cd","Cu","Fe","Pb","Se","Zn")
names(BIG5.LABELS)
BIG5.GLM<-BIG5.LABELS

####As####

mix.int <- glm(log(As+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

pframe <- with(BIG5.GLM,
               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=51),
                           SITE=levels(SITE)))

levels(BIG5.GLM$SITE)

BIG5.GLM$SITE <- factor(BIG5.GLM$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = TRUE)

## add predicted values (on response scale) to prediction frame

pframe$As <- predict(mix.int,newdata=pframe,type="response")

ggplot(BIG5.GLM, aes(SAMPLING_DATE, As, col = SITE))+
  geom_point() +
  scale_y_continuous(trans='log10')+
  geom_smooth(method = "glm", se = FALSE,
              method.args = list(family = "log"), linetype = "dashed")+
  geom_line(data=pframe)  ## use prediction data here


####Cd####

mix.int <- glm(log(Cd+1) ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link="log"))

summary(mix.int)

pframe <- with(BIG5.GLM,
               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=51),
                           SITE=levels(SITE)))

levels(BIG5.GLM$SITE)

BIG5.GLM$SITE <- factor(BIG5.GLM$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = TRUE)

## add predicted values (on response scale) to prediction frame

pframe$Cd <- predict(mix.int,newdata=pframe,type="response")

plot2<- ggplot(BIG5.GLM, aes(SAMPLING_DATE, Cd, col = SITE))+
          geom_point() +
          scale_y_continuous(trans='log10')+
          geom_smooth(method = "glm", se = FALSE,
          method.args = list(family = "log"), linetype = "dashed")+
          geom_line(data=pframe)  ## use prediction data here


