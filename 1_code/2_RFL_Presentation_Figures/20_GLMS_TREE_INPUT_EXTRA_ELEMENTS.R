library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("performance")

#alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata <- data.frame(read_csv("2_incremental/PCA_INPUT_TEST/20230508completedata_SHOTGUN_ALL.csv"))

STANDING.CROP <- alldata

STANDING.CROP[,c(50:88)]<-STANDING.CROP[,c(11:49)]*STANDING.CROP[,7]

STANDING.CROP$SAMPLING_DATE<-as.Date(STANDING.CROP$SAMPLING_DATE, format = "%m/%d/%Y")

COMPARTMENTS <- STANDING.CROP %>% rename_at(vars(names(STANDING.CROP[,c(50:88)])), ~ paste(names(STANDING.CROP[,c(11:49)]),"mg_m2"))

COMPARTMENTS <- COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP" | COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS)<-c("ROWNUM",names(COMPARTMENTS[c(2:ncol(COMPARTMENTS))]))

COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP")] <- lapply(COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP")], factor)

####Extra Elements Data Frame####

names(COMPARTMENTS)

EXTRA.GLM.DF <- na.omit(COMPARTMENTS[,c("SAMPLING_DATE","SITE","FIELD.REP","P mg_m2","Sr mg_m2","Nar mg_m2","Kr mg_m2","Car mg_m2","Sir mg_m2","Ni mg_m2")])

EXTRA.GLM.DF.SUM <- EXTRA.GLM.DF %>% group_by(SAMPLING_DATE, SITE, FIELD.REP) %>% summarise_all(sum)

EXTRA.GLM.DF.SUM <- EXTRA.GLM.DF.SUM %>% rename_at(vars(names(EXTRA.GLM.DF.SUM[,c(4:ncol(EXTRA.GLM.DF.SUM))])), ~ c("P","S","Na","K","Ca","Si","Ni"))

EXTRA.GLM.DF.SUM$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[EXTRA.GLM.DF.SUM$SITE]

EXTRA.GLM.DF.SUM$SAMPLING_DATE<-as.Date(EXTRA.GLM.DF.SUM$SAMPLING_DATE)

####GLMS####

mix.int <- glm(S ~ SAMPLING_DATE * SITE_DISTANCE, data = EXTRA.GLM.DF.SUM, 
               family=gaussian)

summary(mix.int)

coefficients(mix.int)

summary(residuals(mix.int))

pframe <- with(EXTRA.GLM.DF.SUM,
               expand.grid(SAMPLING_DATE=seq(min(SAMPLING_DATE),max(SAMPLING_DATE),length=5),
                           SITE_DISTANCE=seq(min(SITE_DISTANCE),max(SITE_DISTANCE),length=5)))

levels(EXTRA.GLM.DF.SUM$SITE)

## add predicted values (on response scale) to prediction frame

pframe$S <- predict.glm(mix.int,newdata=pframe,type="response")


mix.int$SITE <- factor(mix.int$SITE, levels = c("WS","DL","GR","GC","BG","BN"), ordered = F)

plot2 <- ggplot(mix.int, aes(SAMPLING_DATE, S, col = as.factor(SITE_DISTANCE)))+
  geom_point() +
#  geom_line(data=pframe)+
  scale_y_continuous(trans='pseudo_log', oob=squish)+
  geom_smooth(method="loess" ,se = TRUE,
  method.args = list(family = "gaussian"), linetype = "dashed")
  ## use prediction data here

#dev.new()

#dev.new()

plot2 + facet_grid( ~ as.factor(SITE_DISTANCE))

r2_nagelkerke(mix.int)
