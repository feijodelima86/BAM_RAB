library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("performance")
library("stringr")
library("lubridate")

alldata <- data.frame(read_csv("2_incremental/PCA_INPUT_TEST/20230508completedata_SHOTGUN_ALL.csv"))

#alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv"))

STANDING.CROP <- alldata

STANDING.CROP[,c(49:87)]<-STANDING.CROP[,c(10:48)]*STANDING.CROP[,7]

STANDING.CROP$SAMPLING_DATE<-as.Date(STANDING.CROP$SAMPLING_DATE, format = "%m/%d/%Y")

#STANDING.CROP$SAMPLING_DATE <- format(STANDING.CROP$SAMPLING_DATE, "%m/%Y")

COMPARTMENTS <- STANDING.CROP %>% rename_at(vars(names(STANDING.CROP[,c(49:87)])), ~ paste(names(STANDING.CROP[,c(10:48)]),"mg_m2"))

names(COMPARTMENTS)<-c("ROWNUM",names(COMPARTMENTS[c(2:ncol(COMPARTMENTS))]))

#COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP")] <- lapply(COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP")], factor)

####Extra Elements Data Frame####

EXTRA.GLM.DF <- COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP","OM.AREA.g.m2","As","Cd","Cu","Fer","Mo","Pb","Se","Zn")]

EXTRA.GLM.DF.SUM <- data.frame(EXTRA.GLM.DF %>% group_by(SAMPLING_DATE, SITE, FIELD.REP) %>% summarise_all(sum))

EXTRA.GLM.DF.SUM[,c(5:12)]<-EXTRA.GLM.DF.SUM[,c(5:12)]/EXTRA.GLM.DF.SUM[,4]

#EXTRA.GLM.DF.SUM <- EXTRA.GLM.DF.SUM %>% rename_at(vars(names(EXTRA.GLM.DF.SUM[,c(6:ncol(EXTRA.GLM.DF.SUM))])), ~ c("P","S","Na","K","Ca","Si","Ni"))

EXTRA.GLM.DF.SUM$SITE_DISTANCE <- as.numeric(str_replace_all(EXTRA.GLM.DF.SUM$SITE, c("WS"="0","DL"="44.9","GR"="64.82","GC"="89.22","BG"="144.32","BN"="167.82")))

EXTRA.GLM.DF.SUM$SAMPLING_DATE<-as.Date(EXTRA.GLM.DF.SUM$SAMPLING_DATE, format = "%m/%d/%Y")

#EXTRA.GLM.DF.SUM<-EXTRA.GLM.DF.SUM %>% mutate(SAMPLING_DATE = floor_date(ymd(SAMPLING_DATE), "month"))

EXTRA.GLM.DF.SUM$SAMPLING_DATE[EXTRA.GLM.DF.SUM$SAMPLING_DATE=="2021-09-09"]<-"2021-09-22" # Gimmick

EXTRA.GLM.DF.SUM<-EXTRA.GLM.DF.SUM[complete.cases(EXTRA.GLM.DF.SUM), ]

z_scores <- as.data.frame(sapply(EXTRA.GLM.DF.SUM[,c(6:12)], function(mydata) (abs(mydata-mean(mydata))/sd(mydata))))

EXTRA.GLM.DF.SUM <- EXTRA.GLM.DF.SUM[!rowSums(z_scores>2.5), ]  

####GLMS####

####As####

mix.int <- glm(As ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)



plotAs <- ggplot(mix.int, aes(SITE_DISTANCE, As, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotAs + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotAs <- 
  ggplot(mix.int, aes(SAMPLING_DATE, As, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotAs + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Cd####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Cd ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotCd <- ggplot(mix.int, aes(SITE_DISTANCE, Cd, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotCd + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotCd <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Cd, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotCd + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Cu####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Cu ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotCu <- ggplot(mix.int, aes(SITE_DISTANCE, Cu, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotCu + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotCu <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Cu, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotCu + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Fe####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Fer ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotFer <- ggplot(mix.int, aes(SITE_DISTANCE, Fer, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
 # scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotFer + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotFer <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Fer, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
#  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotFer + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Mo####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Mo ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotMo <- ggplot(mix.int, aes(SITE_DISTANCE, Mo, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
   scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotMo + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotMo <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Mo, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  #  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotMo + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Pb####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Pb ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotPb <- ggplot(mix.int, aes(SITE_DISTANCE, Pb, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotPb + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotPb <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Pb, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  #  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotPb + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Se####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Se ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotSe <- ggplot(mix.int, aes(SITE_DISTANCE, Se, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotSe + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotSe <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Se, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  #  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotSe + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

####Zn####"As","Cd","Cu","Fer","Mo","Pb","Se","Zn"

mix.int <- glm(Zn ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)


plotZn <- ggplot(mix.int, aes(SITE_DISTANCE, Zn, col = as.factor(SAMPLING_DATE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotZn + facet_grid(vars(as.factor(SAMPLING_DATE)))+ theme(legend.position = "none") 

plotZn <- 
  ggplot(mix.int, aes(SAMPLING_DATE, Zn, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
  #  geom_line(data=pframe)+
  #  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = T,
              method.args = list(family = "gaussian"), linetype = "dashed")

plotZn + facet_grid(vars(as.factor(SITE_DISTANCE)))+ theme(legend.position = "none") 

