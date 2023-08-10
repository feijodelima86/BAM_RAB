library(readr)
library(gbm)
library(dismo)
library(MASS)
library(usdm)
library(dplyr)
library(tidyr)
require(tidyverse)
library(stats)
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))

algaedata <- data.frame(read.csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_2.csv"))

waterdata <- data.frame(read.csv("2_incremental/water6.csv"))

turnoverdata <- read_csv("2_incremental/TURNOVER_REDUX.csv")

### prepping algae data ###

names(algaedata)

algaedata$SAMPLING_DATE<-as.Date(algaedata$SAMPLING_DATE, format = "%m/%d/%Y")

algaedata<-algaedata[,-c(1,2)]

### prepping water data ###

names(waterdata)

colnames(waterdata)[c(2:3)]<-c("SAMPLING_DATE","SITE")

waterdata$SITE <- factor(waterdata$SITE, labels = c("WS","DL","GR","GC","BG","BN"))

waterdata$SAMPLING_DATE<-as.Date(waterdata$SAMPLING_DATE, format = "%m/%d/%Y")

waterdata<-waterdata[,-c(1)]

### prepping turnover data ###

names(turnoverdata)

turnoverdata$SAMPLING_DATE<-as.Date(turnoverdata$SAMPLING_DATE, format = "%m/%d/%Y")

### full join reduce ###

add.comp<-list(algaedata,waterdata,turnoverdata)

alldata<-add.comp %>% reduce(full_join)

### aggredating dataset ###

COMPARTMENTS_AVG <-aggregate(x = alldata[,colnames(alldata) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
                             by = list(alldata$SAMPLING_DATE,alldata$SITE,alldata$SAMPLE_DESCRIPTOR),
                             FUN = mean,
                             na.rm = T)
COMPARTMENTS_AVG

names(COMPARTMENTS_AVG)

colnames(COMPARTMENTS_AVG)[c(1:3)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

write.csv(COMPARTMENTS_AVG, "2_incremental/TURNOVER_Full_Dataset.csv")

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset2.csv"))

#### BRT ####

names(COMPARTMENTS_AVG)

### As: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,114,13)]),c(3,4,6,7,10,114,35)]))
### Cd: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,117,19)]),c(3,4,6,7,10,117,19)]))
### Cu: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,112,22)]),c(3,4,6,7,10,112,22)]))
### Fe: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,111,25)]),c(3,4,6,7,10,111,25)]))
### Mo: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,116,31)]),c(3,4,6,7,10,116,31)]))
### Pb: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,118,35)]),c(3,4,6,7,10,118,35)]))
### P:  br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,109,34)]),c(3,4,6,7,10,109,34)]))
### Zn: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,113,49)]),c(3,4,6,7,10,113,49)]))


br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,6,7,10,117,19)]),c(3,4,6,7,10,117,19)]))
br.sub$SAMPLE_DESCRIPTOR<-as.factor(br.sub$SAMPLE_DESCRIPTOR)
br.sub$SITE<-as.factor(br.sub$SITE)
#br.sub[,6]<-log10(br.sub[,6]+1)
#br.sub[,7]<-log10(br.sub[,7]+1)

#br.sub$TOTAL.MG.M2.RITCHIE<-log10(br.sub$TOTAL.MG.M2.RITCHIE+1)
#br.sub$OM.AREA.g.m2<-log10(br.sub$OM.AREA.g.m2+1)

UCFR.CD.tc6.lr002 <- gbm.step(data=br.sub, 
                              gbm.x = c(2,3,4,6),
                              gbm.y = ncol(br.sub),
                              family = "gaussian",
                              tree.complexity = 4,
                              learning.rate = 0.0007,
                              bag.fraction = 0.5
)


predfit<-lm(UCFR.CD.tc6.lr002$fitted ~ br.sub[,7])

summary(predfit)

plot(br.sub[,7], UCFR.CD.tc6.lr002$fitted, 
     #xlim=c(0.5,3), 
     #ylim=c(0.5,3), 
     col=br.sub$SAMPLE_DESCRIPTOR, 
     xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

#dev.new()

gbm.plot(UCFR.CD.tc6.lr002, write.title = F, nplots = 5, plot.layout= c(2,2),
         las=1,
         lwd=2,
         cex.lab=1.5,
         cex.axis=1,
         smooth =F,
         rug=F,
         y.label = NA
)

gbm.plot.fits(UCFR.CD.tc6.lr002)

find.int.cd <- gbm.interactions(UCFR.CD.tc6.lr002)

find.int.cd$interactions

find.int.cd$rank.list

gbm.perspec(UCFR.CD.tc6.lr002,4,2, 
            z.range=c(1,1.62), 
            theta = 315,
            phi=45, 
            cex.lab = 1.2, font.lab = 2, cex.axis = 1, font.axis= 1,
            perspective = T,
            smooth=T)

##### GLM #####

require(scatterplot3d)
library(effects)
library(pscl)

names(br.sub)

tdsp<-scatterplot3d(br.sub[,3],br.sub[,6],br.sub[,7], pch = 19, color = factor(br.sub$SAMPLE_DESCRIPTOR, labels = c("darkgreen","gold","chartreuse")), main="3D Scatterplot",
              type="h", angle = 120)

fitm <- lm(br.sub[,7] ~ br.sub[,3]+br.sub[,6]) 

summary(fitm)

tdsp$plane3d(fitm)

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

mix.int_As_TD <- glm(As ~ (TD.As + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_As_TD))

mix.int_Ca_TD <- glm(Car ~ (TD.Ca + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Ca_TD)

with(summary(mix.int_Ca_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Ca_TD))

mix.int_Cd_TD <- glm(Cd ~ (TD.Cd + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cd_TD))

mix.int_Cu_TD <- glm(Cu ~ (TD.Cu + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cu_TD))

mix.int_Fe_TD <- glm(Fe ~ (TD.Fe + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Fe_TD))

mix.int_Mo_TD <- glm(Mo ~ (TD.Mo + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Mo_TD)

plot(allEffects(mix.int_Mo_TD))

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

mix.int_P_TD <- glm(P ~ (TD.P + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_P_TD))

mix.int_Pb_TD <- glm(Pb ~ (TD.Pb + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD))

mix.int_Se_TD <- glm(Se ~ (TD.Se + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Se_TD))

mix.int_Zn_TD <- glm(Zn ~ (TD.Zn + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Zn_TD))

names(COMPARTMENTS_AVG)


