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

#### Fe ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Fe_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Fe_NH3","colloidal_Fe_NH3","trulydissolved_Fe_NH3","OM.AREA.g.m2","TURNOVER","Fer")])

names(Fe_DF)<-c("COMP","Fe_SPM","Fe_COL","Fe_TRD","OM","TURNOVER", "Fe")

sum(Fe_DF$Fe_SPM, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_COL, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_TRD, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))

Fe_DF<-Fe_DF[Fe_DF$OM>0,]

#Fe_DF$OM <- log10(Fe_DF$OM+1)

Fe_DF$TURNOVER <- sqrt(Fe_DF$TURNOVER)

Fe_DF <- Fe_DF[complete.cases(Fe_DF),]

max.VIF<-5

Fe.RED <- data.frame(scale(Fe_DF[, c(2:(ncol(Fe_DF)-1))],center=T, scale=T))

Fe.RED <- replace_outliers_with_na(Fe.RED)

VIF.Fe<-vif(Fe.RED)

repeat {   
  if (max(VIF.Fe[,2])>max.VIF) {
    VIF.Fe<-vif(Fe.RED)
    Fe.RED <- within(Fe.RED, rm(list=VIF.Fe[which.max(VIF.Fe[,2]),1]))
  } else {
    break
  }
}

Fe_DF<-data.frame(Fe_DF[,1],Fe.RED,Fe_DF[,ncol(Fe_DF)])

Fe_DF <- Fe_DF[complete.cases(Fe_DF),]

names(Fe_DF)[1]<-"COMP"

names(Fe_DF)[ncol(Fe_DF)]<-"Fe"

Fe.full.model <- glm(Fe ~. ^2, data = Fe_DF, family = gaussian)
coef(Fe.full.model)
Anova(Fe.full.model, test="LR")
with(summary(Fe.full.model), 1 - deviance/null.deviance)

Fe.step.model <- Fe.full.model %>% stepAIC(trace = T)
coef(Fe.step.model)
Anova(Fe.step.model)
summary(Fe.step.model)
with(summary(Fe.step.model), 1 - deviance/null.deviance)

Fe.final.model <- glm(Fe ~ (TURNOVER+COMP)^2, data = Fe_DF, family = gaussian)
coef(Fe.final.model)
Anova(Fe.final.model, test="LR")
with(summary(Fe.final.model), 1 - deviance/null.deviance)

plot(allEffects(Fe.final.model))

plot(Fe_DF$Fe_TRD, Fe_DF$Fe)

plot(Fe_DF$Fe_COL, Fe_DF$Fe)

plot(Fe_DF$Fe_SPM, Fe_DF$Fe)

plot(Fe_DF$TURNOVER, Fe_DF$Fe)

plot(Fe_DF$OM, Fe_DF$Fe)

calc.relimp(Fe.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Fe_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Fe_NH3","colloidal_Fe_NH3","trulydissolved_Fe_NH3","OM.AREA.g.m2","TURNOVER","Fer")])

names(Fe_DF)<-c("COMP","Fe_SPM","Fe_COL","Fe_TRD","OM","TURNOVER", "Fe")

sum(Fe_DF$Fe_SPM, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_COL, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_TRD, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))

Fe_DF<-Fe_DF[Fe_DF$OM>0,]

#Fe_DF$OM <- log10(Fe_DF$OM+1)

Fe_DF$TURNOVER <- sqrt(Fe_DF$TURNOVER)

Fe_DF <- Fe_DF[complete.cases(Fe_DF),]

max.VIF<-5

Fe.RED <- data.frame(scale(Fe_DF[, c(2:(ncol(Fe_DF)-1))],center=T, scale=T))

Fe.RED <- replace_outliers_with_na(Fe.RED)

VIF.Fe<-vif(Fe.RED)

repeat {   
  if (max(VIF.Fe[,2])>max.VIF) {
    VIF.Fe<-vif(Fe.RED)
    Fe.RED <- within(Fe.RED, rm(list=VIF.Fe[which.max(VIF.Fe[,2]),1]))
  } else {
    break
  }
}

Fe_DF<-data.frame(Fe_DF[,1],Fe.RED,Fe_DF[,ncol(Fe_DF)])

Fe_DF <- Fe_DF[complete.cases(Fe_DF),]

names(Fe_DF)[1]<-"COMP"

names(Fe_DF)[ncol(Fe_DF)]<-"Fe"

Fe_DF <- Fe_DF[,-1]

Fe_DF <- Fe_DF[,-3]

Fe.full.model <- glm(Fe ~. ^2, data = Fe_DF, family = gaussian)
coef(Fe.full.model)
Anova(Fe.full.model, test="LR")
with(summary(Fe.full.model), 1 - deviance/null.deviance)

Fe.step.model <- Fe.full.model %>% stepAIC(trace = T)
coef(Fe.step.model)
Anova(Fe.step.model)
summary(Fe.step.model)
with(summary(Fe.step.model), 1 - deviance/null.deviance)

Fe.final.model <- glm(Fe ~ (OM)^2, data = Fe_DF, family = gaussian)
coef(Fe.final.model)
Anova(Fe.final.model, test="LR")
summary(Fe.final.model)
with(summary(Fe.final.model), 1 - deviance/null.deviance)

plot(allEffects(Fe.final.model))

plot(Fe_DF$Fe_TRD, Fe_DF$Fe)

plot(Fe_DF$Fe_COL, Fe_DF$Fe)

plot(Fe_DF$Fe_SPM, Fe_DF$Fe)

plot(Fe_DF$TURNOVER, Fe_DF$Fe)

plot(Fe_DF$OM, Fe_DF$Fe)

calc.relimp(Fe.final.model)

