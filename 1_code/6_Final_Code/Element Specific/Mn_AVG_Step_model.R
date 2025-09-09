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

#### Mn ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Mn_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Mn","colloidal_Mn","trulydissolved_Mn","OM.AREA.g.m2","TURNOVER","Mn")])

names(Mn_DF)<-c("COMP","Mn_SPM","Mn_COL","Mn_TRD","OM","TURNOVER", "Mn")

sum(Mn_DF$Mn_SPM, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_COL, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_TRD, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))

Mn_DF<-Mn_DF[Mn_DF$OM>0,]

#Mn_DF$OM <- log10(Mn_DF$OM+1)

Mn_DF$TURNOVER <- sqrt(Mn_DF$TURNOVER)

Mn_DF <- Mn_DF[complete.cases(Mn_DF),]

max.VIF<-5

Mn.RED <- data.frame(scale(Mn_DF[, c(2:(ncol(Mn_DF)-1))],center=T, scale=T))

Mn.RED <- replace_outliers_with_na(Mn.RED)

VIF.Mn<-vif(Mn.RED)

repeat {   
  if (max(VIF.Mn[,2])>max.VIF) {
    VIF.Mn<-vif(Mn.RED)
    Mn.RED <- within(Mn.RED, rm(list=VIF.Mn[which.max(VIF.Mn[,2]),1]))
  } else {
    break
  }
}

Mn_DF<-data.frame(Mn_DF[,1],Mn.RED,Mn_DF[,ncol(Mn_DF)])

Mn_DF <- Mn_DF[complete.cases(Mn_DF),]

names(Mn_DF)[1]<-"COMP"

names(Mn_DF)[ncol(Mn_DF)]<-"Mn"

Mn.full.model <- glm(Mn ~. ^2, data = Mn_DF, family = gaussian)
coef(Mn.full.model)
Anova(Mn.full.model, test="LR")
with(summary(Mn.full.model), 1 - deviance/null.deviance)

Mn.step.model <- Mn.full.model %>% stepAIC(trace = T)
coef(Mn.step.model)
Anova(Mn.step.model)
summary(Mn.step.model)
with(summary(Mn.step.model), 1 - deviance/null.deviance)

Mn.final.model <- glm(Mn ~ (Mn_TRD+ OM)^2, data = Mn_DF, family = gaussian)
coef(Mn.final.model)
Anova(Mn.final.model, test="LR")
with(summary(Mn.final.model), 1 - deviance/null.deviance)

plot(allEffects(Mn.final.model))

plot(Mn_DF$Mn_TRD, Mn_DF$Mn)

plot(Mn_DF$Mn_COL, Mn_DF$Mn)

plot(Mn_DF$Mn_SPM, Mn_DF$Mn)

plot(Mn_DF$TURNOVER, Mn_DF$Mn)

plot(Mn_DF$OM, Mn_DF$Mn)

calc.relimp(Mn.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Mn_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Mn","colloidal_Mn","trulydissolved_Mn","OM.AREA.g.m2","TURNOVER","Mn")])

names(Mn_DF)<-c("COMP","Mn_SPM","Mn_COL","Mn_TRD","OM","TURNOVER", "Mn")

sum(Mn_DF$Mn_SPM, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_COL, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_TRD, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))

Mn_DF<-Mn_DF[Mn_DF$OM>0,]

#Mn_DF$OM <- log10(Mn_DF$OM+1)

Mn_DF$TURNOVER <- sqrt(Mn_DF$TURNOVER)

Mn_DF <- Mn_DF[complete.cases(Mn_DF),]

max.VIF<-5

Mn.RED <- data.frame(scale(Mn_DF[, c(2:(ncol(Mn_DF)-1))],center=T, scale=T))

Mn.RED <- replace_outliers_with_na(Mn.RED)

VIF.Mn<-vif(Mn.RED)

repeat {   
  if (max(VIF.Mn[,2])>max.VIF) {
    VIF.Mn<-vif(Mn.RED)
    Mn.RED <- within(Mn.RED, rm(list=VIF.Mn[which.max(VIF.Mn[,2]),1]))
  } else {
    break
  }
}

Mn_DF<-data.frame(Mn_DF[,1],Mn.RED,Mn_DF[,ncol(Mn_DF)])

Mn_DF <- Mn_DF[complete.cases(Mn_DF),]

names(Mn_DF)[1]<-"COMP"

names(Mn_DF)[ncol(Mn_DF)]<-"Mn"

Mn_DF <- Mn_DF[,-1]

Mn.full.model <- glm(Mn ~. ^2, data = Mn_DF, family = gaussian)
coef(Mn.full.model)
Anova(Mn.full.model, test="LR")
with(summary(Mn.full.model), 1 - deviance/null.deviance)

Mn.step.model <- Mn.full.model %>% stepAIC(trace = T)
coef(Mn.step.model)
Anova(Mn.step.model)
summary(Mn.step.model)
with(summary(Mn.step.model), 1 - deviance/null.deviance)

Mn.final.model <- glm(Mn ~ (Mn_SPM + OM), data = Mn_DF, family = gaussian)
coef(Mn.final.model)
Anova(Mn.final.model, test="LR")
summary(Mn.final.model)
with(summary(Mn.final.model), 1 - deviance/null.deviance)

plot(allEffects(Mn.final.model))

plot(Mn_DF$Mn_TRD, Mn_DF$Mn)

plot(Mn_DF$Mn_COL, Mn_DF$Mn)

plot(Mn_DF$Mn_SPM, Mn_DF$Mn)

plot(Mn_DF$TURNOVER, Mn_DF$Mn)

plot(Mn_DF$OM, Mn_DF$Mn)

calc.relimp(Mn.final.model)

