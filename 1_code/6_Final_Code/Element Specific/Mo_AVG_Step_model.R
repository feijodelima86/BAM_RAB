library(readr)
library(effects)
library(AICcmodavg)
library(nlme)
library(car)
library(interactions)
library(ggplot2)
library(lme4)
library(plotly)
library(stringr)
library(reshape2)
library(gridExtra)
library(MASS)
library(relaimpo)
library(usdm)

#### Housekeeping

replace_outliers_with_na <- function(x, threshold = 3) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}

## Data Onborading ###

#### Mo ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Mo_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Mo_Oshift","colloidal_Mo_Oshift","trulydissolved_Mo_Oshift","OM.AREA.g.m2","TURNOVER","Mo")])

names(Mo_DF)<-c("COMP","Mo_SPM","Mo_COL","Mo_TRD","OM","TURNOVER", "Mo")

sum(Mo_DF$Mo_SPM, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_COL, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_TRD, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))

Mo_DF<-Mo_DF[Mo_DF$OM>0,]

#Mo_DF$OM <- log10(Mo_DF$OM+1)

Mo_DF$TURNOVER <- sqrt(Mo_DF$TURNOVER)

Mo_DF <- Mo_DF[complete.cases(Mo_DF),]

max.VIF<-5

Mo.RED <- data.frame(scale(Mo_DF[, c(2:(ncol(Mo_DF)-1))],center=T, scale=T))

Mo.RED <- replace_outliers_with_na(Mo.RED)

VIF.Mo<-vif(Mo.RED)

repeat {   
  if (max(VIF.Mo[,2])>max.VIF) {
    VIF.Mo<-vif(Mo.RED)
    Mo.RED <- within(Mo.RED, rm(list=VIF.Mo[which.max(VIF.Mo[,2]),1]))
  } else {
    break
  }
}

Mo_DF<-data.frame(Mo_DF[,1],Mo.RED,Mo_DF[,ncol(Mo_DF)])

Mo_DF <- Mo_DF[complete.cases(Mo_DF),]

names(Mo_DF)[1]<-"COMP"

names(Mo_DF)[ncol(Mo_DF)]<-"Mo"

Mo.full.model <- glm(Mo ~. ^2, data = Mo_DF, family = gaussian)
coef(Mo.full.model)
Anova(Mo.full.model, test="LR")
with(summary(Mo.full.model), 1 - deviance/null.deviance)

Mo.step.model <- Mo.full.model %>% stepAIC(trace = T)
coef(Mo.step.model)
Anova(Mo.step.model)
summary(Mo.step.model)
with(summary(Mo.step.model), 1 - deviance/null.deviance)

Mo.final.model <- glm(Mo ~ (Mo_TRD+OM+TURNOVER)^2, data = Mo_DF, family = gaussian)
coef(Mo.final.model)
Anova(Mo.final.model, test="LR")
with(summary(Mo.final.model), 1 - deviance/null.deviance)

plot(allEffects(Mo.final.model))

plot(Mo_DF$Mo_TRD, Mo_DF$Mo)

plot(Mo_DF$Mo_COL, Mo_DF$Mo)

plot(Mo_DF$Mo_SPM, Mo_DF$Mo)

plot(Mo_DF$TURNOVER, Mo_DF$Mo)

plot(Mo_DF$OM, Mo_DF$Mo)

calc.relimp(Mo.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Mo_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Mo_Oshift","colloidal_Mo_Oshift","trulydissolved_Mo_Oshift","OM.AREA.g.m2","TURNOVER","Mo")])

names(Mo_DF)<-c("COMP","Mo_SPM","Mo_COL","Mo_TRD","OM","TURNOVER", "Mo")

sum(Mo_DF$Mo_SPM, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_COL, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_TRD, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))

Mo_DF<-Mo_DF[Mo_DF$OM>0,]

#Mo_DF$OM <- log10(Mo_DF$OM+1)

Mo_DF$TURNOVER <- sqrt(Mo_DF$TURNOVER)

Mo_DF <- Mo_DF[complete.cases(Mo_DF),]

max.VIF<-5

Mo.RED <- data.frame(scale(Mo_DF[, c(2:(ncol(Mo_DF)-1))],center=T, scale=T))

Mo.RED <- replace_outliers_with_na(Mo.RED)

VIF.Mo<-vif(Mo.RED)

repeat {   
  if (max(VIF.Mo[,2])>max.VIF) {
    VIF.Mo<-vif(Mo.RED)
    Mo.RED <- within(Mo.RED, rm(list=VIF.Mo[which.max(VIF.Mo[,2]),1]))
  } else {
    break
  }
}

Mo_DF<-data.frame(Mo_DF[,1],Mo.RED,Mo_DF[,ncol(Mo_DF)])

Mo_DF <- Mo_DF[complete.cases(Mo_DF),]

names(Mo_DF)[1]<-"COMP"

names(Mo_DF)[ncol(Mo_DF)]<-"Mo"

Mo_DF <- Mo_DF[,-1]

Mo.full.model <- glm(Mo ~. ^2, data = Mo_DF, family = gaussian)
coef(Mo.full.model)
Anova(Mo.full.model, test="LR")
with(summary(Mo.full.model), 1 - deviance/null.deviance)

Mo.step.model <- Mo.full.model %>% stepAIC(trace = T)
coef(Mo.step.model)
Anova(Mo.step.model)
summary(Mo.step.model)
with(summary(Mo.step.model), 1 - deviance/null.deviance)

Mo.final.model <- glm(Mo ~ (OM), data = Mo_DF, family = gaussian)
coef(Mo.final.model)
Anova(Mo.final.model, test="LR")
summary(Mo.final.model)
with(summary(Mo.final.model), 1 - deviance/null.deviance)

plot(allEffects(Mo.final.model))

plot(Mo_DF$Mo_TRD, Mo_DF$Mo)

plot(Mo_DF$Mo_COL, Mo_DF$Mo)

plot(Mo_DF$Mo_SPM, Mo_DF$Mo)

plot(Mo_DF$TURNOVER, Mo_DF$Mo)

plot(Mo_DF$OM, Mo_DF$Mo)

calc.relimp(Mo.final.model)

