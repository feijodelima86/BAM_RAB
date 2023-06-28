library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("performance")
library("splines")

#alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
alldata <- read_csv("2_incremental/PCA_INPUT_TEST/20230505completedata_SHOTGUN_ALL.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS)<-c("ROWNUM",names(COMPARTMENTS[c(2:ncol(COMPARTMENTS))]))
  
COMPARTMENTS$ROWNUM==1

####Big 5+Fe+Se####

#+39

BIG5.LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,5,6,7),c(13,19,22,25,35,40))]))

BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)

#BIG5.LABELS<-BIG5.LABELS[order(BIG5.LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]

names(BIG5.LABELS)[c(6:ncol(BIG5.LABELS))]<- c("As","Cd","Cu","Fe","Pb","Se","Zn")


# summarising all non-grouping variables

BIG5.LABELS<-distinct(BIG5.LABELS, SAMPLING_DATE, SITE,FIELD.REP, .keep_all= TRUE)

BIG5.LABELS[,c(5:12)] %>% mutate_if(is.character, as.numeric)
BIG5.LABELS[,3]<-as.character(BIG5.LABELS[,3])
BIG5.GLM <- BIG5.LABELS[,-4] %>% group_by(SAMPLING_DATE, SITE, FIELD.REP) %>% summarise_all(sum)

####GLMS####

mix.int <- glm(Fe ~ SAMPLING_DATE * SITE, data = BIG5.GLM, 
               family=gaussian)

summary(mix.int)

coefficients(mix.int)

summary(residuals(mix.int))

pframe <- with(BIG5.GLM,
               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=2),
                           SITE=levels(SITE)))

levels(BIG5.GLM$SITE)

## add predicted values (on response scale) to prediction frame

pframe$Fe <- predict.glm(mix.int,newdata=pframe,type="response")


BIG5.GLM$SITE <- factor(BIG5.GLM$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = F)

plot2<- ggplot(BIG5.GLM, aes(SAMPLING_DATE, Fe, col = SITE))+
  geom_point() +
#  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
  method.args = list(family = "gaussian"), linetype = "dashed")
  ## use prediction data here

#dev.new()

plot2 + facet_grid( ~ SITE)

r2_nagelkerke(mix.int)
