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

#### Cu ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Cu_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Cu_NH3","colloidal_Cu_NH3","trulydissolved_Cu_NH3","OM.AREA.g.m2","TURNOVER","Cu")])

names(Cu_DF)<-c("COMP","Cu_SPM","Cu_COL","Cu_TRD","OM","TURNOVER", "Cu")

sum(Cu_DF$Cu_SPM, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_COL, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_TRD, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))

Cu_DF<-Cu_DF[Cu_DF$OM>0,]

#Cu_DF$OM <- log10(Cu_DF$OM+1)

Cu_DF$TURNOVER <- sqrt(Cu_DF$TURNOVER)

Cu_DF <- Cu_DF[complete.cases(Cu_DF),]

max.VIF<-5

Cu.RED <- data.frame(scale(Cu_DF[, c(2:(ncol(Cu_DF)-1))],center=T, scale=T))

Cu.RED <- replace_outliers_with_na(Cu.RED)

VIF.Cu<-vif(Cu.RED)

repeat {   
  if (max(VIF.Cu[,2])>max.VIF) {
    VIF.Cu<-vif(Cu.RED)
    Cu.RED <- within(Cu.RED, rm(list=VIF.Cu[which.max(VIF.Cu[,2]),1]))
  } else {
    break
  }
}

Cu_DF<-data.frame(Cu_DF[,1],Cu.RED,Cu_DF[,ncol(Cu_DF)])

Cu_DF <- Cu_DF[complete.cases(Cu_DF),]

names(Cu_DF)[1]<-"COMP"

names(Cu_DF)[ncol(Cu_DF)]<-"Cu"

Cu.full.model <- glm(Cu ~. ^2, data = Cu_DF, family = gaussian)
coef(Cu.full.model)
Anova(Cu.full.model, test="LR")
with(summary(Cu.full.model), 1 - deviance/null.deviance)

Cu.step.model <- Cu.full.model %>% stepAIC(trace = T)
coef(Cu.step.model)
Anova(Cu.step.model)
summary(Cu.step.model)
with(summary(Cu.step.model), 1 - deviance/null.deviance)

Cu.final.model <- glm(Cu ~ (Cu_COL + OM +TURNOVER), data = Cu_DF, family = gaussian)
coef(Cu.final.model)
Anova(Cu.final.model, test="LR")
with(summary(Cu.final.model), 1 - deviance/null.deviance)

plot(allEffects(Cu.final.model))

plot(Cu_DF$Cu_TRD, Cu_DF$Cu)

plot(Cu_DF$Cu_COL, Cu_DF$Cu)

plot(Cu_DF$Cu_SPM, Cu_DF$Cu)

plot(Cu_DF$TURNOVER, Cu_DF$Cu)

plot(Cu_DF$OM, Cu_DF$Cu)

calc.relimp(Cu.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Cu_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Cu_NH3","colloidal_Cu_NH3","trulydissolved_Cu_NH3","OM.AREA.g.m2","TURNOVER","Cu")])

names(Cu_DF)<-c("COMP","Cu_SPM","Cu_COL","Cu_TRD","OM","TURNOVER", "Cu")

sum(Cu_DF$Cu_SPM, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_COL, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_TRD, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))

Cu_DF<-Cu_DF[Cu_DF$OM>0,]

#Cu_DF$OM <- log10(Cu_DF$OM+1)

Cu_DF$TURNOVER <- sqrt(Cu_DF$TURNOVER)

Cu_DF <- Cu_DF[complete.cases(Cu_DF),]

max.VIF<-5

Cu.RED <- data.frame(scale(Cu_DF[, c(2:(ncol(Cu_DF)-1))],center=T, scale=T))

Cu.RED <- replace_outliers_with_na(Cu.RED)

VIF.Cu<-vif(Cu.RED)

repeat {   
  if (max(VIF.Cu[,2])>max.VIF) {
    VIF.Cu<-vif(Cu.RED)
    Cu.RED <- within(Cu.RED, rm(list=VIF.Cu[which.max(VIF.Cu[,2]),1]))
  } else {
    break
  }
}

Cu_DF<-data.frame(Cu_DF[,1],Cu.RED,Cu_DF[,ncol(Cu_DF)])

Cu_DF <- Cu_DF[complete.cases(Cu_DF),]

names(Cu_DF)[1]<-"COMP"

names(Cu_DF)[ncol(Cu_DF)]<-"Cu"

Cu_DF <- Cu_DF[,-1]

Cu.full.model <- glm(Cu ~. ^2, data = Cu_DF, family = gaussian)
coef(Cu.full.model)
Anova(Cu.full.model, test="LR")
with(summary(Cu.full.model), 1 - deviance/null.deviance)

Cu.step.model <- Cu.full.model %>% stepAIC(trace = T)
coef(Cu.step.model)
Anova(Cu.step.model)
summary(Cu.step.model)
with(summary(Cu.step.model), 1 - deviance/null.deviance)

Cu.final.model <- glm(Cu ~ (Cu_SPM*Cu_COL*OM), data = Cu_DF, family = gaussian)
coef(Cu.final.model)
Anova(Cu.final.model, test="LR")
summary(Cu.final.model)
with(summary(Cu.final.model), 1 - deviance/null.deviance)

plot(allEffects(Cu.final.model))

plot(Cu_DF$Cu_TRD, Cu_DF$Cu)

plot(Cu_DF$Cu_COL, Cu_DF$Cu)

plot(Cu_DF$Cu_SPM, Cu_DF$Cu)

plot(Cu_DF$TURNOVER, Cu_DF$Cu)

plot(Cu_DF$OM, Cu_DF$Cu)

calc.relimp(Cu.final.model)

