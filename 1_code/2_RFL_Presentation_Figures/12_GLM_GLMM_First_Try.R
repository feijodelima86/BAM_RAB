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

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]
#COMPARTMENTS <- COMPARTMENTS[-217,]

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

#COMPARTMENTS <- COMPARTMENTS[which(COMPARTMENTS$SITE != "WS"),]

COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]

#COMPARTMENTS<-COMPARTMENTS[-1,]#<- removed outlier*




colnames(COMPARTMENTS)<-"ROWNUM"
COMPARTMENTS <- COMPARTMENTS[COMPARTMENTS$NA != "1"]

####Big 5+Fe+Se####

#+39

BIG5.LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,6,89),c(13,19,22,25,35,40,49))]))
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
#BIG5.LABELS<-BIG5.LABELS[order(BIG5.LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
names(BIG5.LABELS)[c(5:11)]<- c("As","Cd","Cu","Fe","Pb","Se","Zn")

BIG5.GLM<-BIG5.LABELS

####Cd####

mix.int <- glm(As ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian(link=log))

summary(mix.int)


pframe <- with(BIG5.GLM,
               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=2),
                           SITE=levels(SITE)))

levels(BIG5.GLM$SITE)


## add predicted values (on response scale) to prediction frame

pframe$As <- predict.glm(mix.int,newdata=pframe,type="response")

BIG5.GLM$SITE <- factor(BIG5.GLM$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = F)

plot2<- ggplot(BIG5.GLM, aes(SAMPLING_DATE, log10(As), col = SITE))+
          geom_point() +
          geom_line(data=pframe)+
#          coord_trans(x="log10") 
#         geom_smooth(method = "glm", se = TRUE,
#          method.args = list(family = "gaussian"), linetype = "dashed")
#ylim(-1.5, max(BIG5.GLM$As))+
#          

## use prediction data here
dev.new()

plot2 + facet_grid( ~ SITE)


