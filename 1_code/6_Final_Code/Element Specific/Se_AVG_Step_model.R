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

#### Se ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Se_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","W_Se82_Oshift","F_Se82_Oshift","trulydissolved_Se82_Oshift","OM.AREA.g.m2","TURNOVER","Se")])

names(Se_DF)<-c("COMP","Se_SPM","Se_COL","Se_TRD","OM","TURNOVER", "Se")

sum(Se_DF$Se_SPM, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_COL, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_TRD, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))

Se_DF<-Se_DF[Se_DF$OM>0,]

#Se_DF$OM <- log10(Se_DF$OM+1)

Se_DF$TURNOVER <- sqrt(Se_DF$TURNOVER)

Se_DF <- Se_DF[complete.cases(Se_DF),]

max.VIF<-5

Se.RED <- data.frame(scale(Se_DF[, c(2:(ncol(Se_DF)-1))],center=T, scale=T))

Se.RED <- replace_outliers_with_na(Se.RED)

VIF.Se<-vif(Se.RED)

repeat {   
  if (max(VIF.Se[,2])>max.VIF) {
    VIF.Se<-vif(Se.RED)
    Se.RED <- within(Se.RED, rm(list=VIF.Se[which.max(VIF.Se[,2]),1]))
  } else {
    break
  }
}

Se_DF<-data.frame(Se_DF[,1],Se.RED,Se_DF[,ncol(Se_DF)])

Se_DF <- Se_DF[complete.cases(Se_DF),]

names(Se_DF)[1]<-"COMP"

names(Se_DF)[ncol(Se_DF)]<-"Se"

Se.full.model <- glm(Se ~. ^2, data = Se_DF, family = gaussian)
coef(Se.full.model)
Anova(Se.full.model, test="LR")
with(summary(Se.full.model), 1 - deviance/null.deviance)

Se.step.model <- Se.full.model %>% stepAIC(trace = T)
coef(Se.step.model)
Anova(Se.step.model)
summary(Se.step.model)
with(summary(Se.step.model), 1 - deviance/null.deviance)

Se.final.model <- glm(Se ~ (COMP*OM*TURNOVER)^2, data = Se_DF, family = gaussian)
coef(Se.final.model)
Anova(Se.final.model, test="LR")
with(summary(Se.final.model), 1 - deviance/null.deviance)

plot(allEffects(Se.final.model))

plot(Se_DF$Se_TRD, Se_DF$Se, col=Se_DF$COMP)

plot(Se_DF$Se_COL, Se_DF$Se, col=Se_DF$COMP)

plot(Se_DF$Se_SPM, Se_DF$Se, col=Se_DF$COMP)

plot(Se_DF$TURNOVER, Se_DF$Se, col=Se_DF$COMP)

plot(Se_DF$OM, Se_DF$Se, col=Se_DF$COMP)

calc.relimp(Se.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Se_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Se82_Oshift","colloidal_Se82_Oshift","trulydissolved_Se82_Oshift","OM.AREA.g.m2","TURNOVER","Se")])

names(Se_DF)<-c("COMP","Se_SPM","Se_COL","Se_TRD","OM","TURNOVER", "Se")

sum(Se_DF$Se_SPM, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_COL, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_TRD, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))

Se_DF<-Se_DF[Se_DF$OM>0,]

#Se_DF$OM <- log10(Se_DF$OM+1)

Se_DF$TURNOVER <- sqrt(Se_DF$TURNOVER)

Se_DF <- Se_DF[complete.cases(Se_DF),]

max.VIF<-5

Se.RED <- data.frame(scale(Se_DF[, c(2:(ncol(Se_DF)-1))],center=T, scale=T))

Se.RED <- replace_outliers_with_na(Se.RED)

VIF.Se<-vif(Se.RED)

repeat {   
  if (max(VIF.Se[,2])>max.VIF) {
    VIF.Se<-vif(Se.RED)
    Se.RED <- within(Se.RED, rm(list=VIF.Se[which.max(VIF.Se[,2]),1]))
  } else {
    break
  }
}

Se_DF<-data.frame(Se_DF[,1],Se.RED,Se_DF[,ncol(Se_DF)])

Se_DF <- Se_DF[complete.cases(Se_DF),]

names(Se_DF)[1]<-"COMP"

names(Se_DF)[ncol(Se_DF)]<-"Se"

Se_DF <- Se_DF[,-1]

Se_DF <- Se_DF[,-2]

Se.full.model <- glm(Se ~. ^2, data = Se_DF, family = gaussian)
coef(Se.full.model)
Anova(Se.full.model, test="LR")
with(summary(Se.full.model), 1 - deviance/null.deviance)

Se.step.model <- Se.full.model %>% stepAIC(trace = T)
coef(Se.step.model)
Anova(Se.step.model)
summary(Se.step.model)
with(summary(Se.step.model), 1 - deviance/null.deviance)

Se.final.model <- glm(Se ~ (Se_TRD + OM + TURNOVER)^2, data = Se_DF, family = gaussian)
coef(Se.final.model)
Anova(Se.final.model, test="LR")
summary(Se.final.model)
with(summary(Se.final.model), 1 - deviance/null.deviance)

plot(allEffects(Se.final.model))

plot(Se_DF$Se_TRD, Se_DF$Se)

plot(Se_DF$Se_COL, Se_DF$Se)

plot(Se_DF$Se_SPM, Se_DF$Se)

plot(Se_DF$TURNOVER, Se_DF$Se)

plot(Se_DF$OM, Se_DF$Se)

calc.relimp(Se.final.model)

