library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("AICcmodavg")
library("performance")
library("splines")

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
#alldata <- read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv")

names(alldata)

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

#COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS)<-c("ROWNUM",names(COMPARTMENTS[c(2:ncol(COMPARTMENTS))]))

COMPARTMENTS$ROWNUM==1

#COMPARTMENTS <- COMPARTMENTS[COMPARTMENTS$ROWNUM != 1, ]

COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]


#COMPARTMENTS <- COMPARTMENTS[which(COMPARTMENTS$SITE != "WS"),]

#COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]

#COMPARTMENTS<-COMPARTMENTS[-1,]#<- removed outlier*

####Big 5+Fe+Se####

#+39
names(BIG5.LABELS)

BIG5.LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,5,6,7,9,10,89),c(13,19,22,25,35,40,49))]))
BIG5.LABELS[,c(5:ncol(BIG5.LABELS))] %>% mutate_if(is.character, as.numeric)
BIG5.LABELS[,3]<-as.character(BIG5.LABELS[,3])
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
names(BIG5.LABELS)[c(9:ncol(BIG5.LABELS))]<- c("As","Cd","Cu","Fe","Pb","Se","Zn")
BIG5.GLM<-BIG5.LABELS

# summarising all non-grouping variables

#BIG5.GLM <- BIG5.LABELS[,-4] %>% group_by(SAMPLING_DATE, SITE, FIELD.REP) %>% summarise_all(sum)
#BIG5.LABELS<-distinct(BIG5.LABELS, SAMPLING_DATE, SITE,FIELD.REP, .keep_all= TRUE)


####GLMS####

mix.int.1 <- glm(Se ~ SAMPLING_DATE * SITE_DISTANCE, data = BIG5.GLM, 
               family=gaussian)

mix.int.2 <- glm(Se ~ ns(SAMPLING_DATE,2) * SITE_DISTANCE, data = BIG5.GLM, 
               family=gaussian)

mix.int.3 <- glm(Se ~ ns(SAMPLING_DATE,3) * SITE_DISTANCE, data = BIG5.GLM, 
                 family=gaussian)

mix.int.4 <- glm(Se ~ ns(SAMPLING_DATE,4) * SITE_DISTANCE, data = BIG5.GLM, 
                 family=gaussian)


mod.names <- c('mix.int.1', 'mix.int.2', 'mix.int.3', 'mix.int.4')

models <- list(mix.int.1, mix.int.2, mix.int.3, mix.int.4)

aictab(cand.set = models, modnames = mod.names)

mix.int<-mix.int.1

summary(mix.int.1)

coefficients(mix.int)

summary(residuals(mix.int))

#pframe <- with(BIG5.GLM,
#               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=2),
#                           SITE=levels(SITE)))


pframe <- with(BIG5.GLM,
               expand.grid(SAMPLING_DATE=SAMPLING_DATE,
                           SITE_DISTANCE=SITE_DISTANCE,
                           SITE=levels(SITE)))


levels(BIG5.GLM$SITE)

## add predicted values (on response scale) to prediction frame

pframe$Se <- predict.glm(mix.int,newdata=pframe,type="response")

BIG5.GLM$SITE <- factor(BIG5.GLM$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = F)

plot2<- ggplot(BIG5.GLM, aes(SAMPLING_DATE, Se, col = SITE_DISTANCE))+
  geom_point() +
  geom_line(data=pframe, size=1.25)+
  scale_y_continuous(limits = c(NA,max(COMPARTMENTS[COMPARTMENTS$ROWNUM != 1, ]$Se, na.rm = T)), 
#  scale_y_continuous(limits = c(NA,max(BIG5.GLM$As, na.rm = T)), 
                     trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")
dev.new()

plot2 + facet_grid( ~ SITE_DISTANCE)

r2_nagelkerke(mix.int)

