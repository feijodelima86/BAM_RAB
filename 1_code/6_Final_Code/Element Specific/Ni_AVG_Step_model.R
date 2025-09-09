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

#### Ni ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Ni_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Ni","colloidal_Ni","trulydissolved_Ni","OM.AREA.g.m2","TURNOVER","Ni")])

names(Ni_DF)<-c("COMP","Ni_SPM","Ni_COL","Ni_TRD","OM","TURNOVER", "Ni")

sum(Ni_DF$Ni_SPM, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_COL, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_TRD, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))

Ni_DF<-Ni_DF[Ni_DF$OM>0,]

#Ni_DF$OM <- log10(Ni_DF$OM+1)

Ni_DF$TURNOVER <- sqrt(Ni_DF$TURNOVER)

Ni_DF <- Ni_DF[complete.cases(Ni_DF),]

max.VIF<-5

Ni.RED <- data.frame(scale(Ni_DF[, c(2:(ncol(Ni_DF)-1))],center=T, scale=T))

Ni.RED <- replace_outliers_with_na(Ni.RED)

VIF.Ni<-vif(Ni.RED)

repeat {   
  if (max(VIF.Ni[,2])>max.VIF) {
    VIF.Ni<-vif(Ni.RED)
    Ni.RED <- within(Ni.RED, rm(list=VIF.Ni[which.max(VIF.Ni[,2]),1]))
  } else {
    break
  }
}

Ni_DF<-data.frame(Ni_DF[,1],Ni.RED,Ni_DF[,ncol(Ni_DF)])

Ni_DF <- Ni_DF[complete.cases(Ni_DF),]

names(Ni_DF)[1]<-"COMP"

names(Ni_DF)[ncol(Ni_DF)]<-"Ni"

Ni.full.model <- glm(Ni ~. ^2, data = Ni_DF, family = gaussian)
coef(Ni.full.model)
Anova(Ni.full.model, test="LR")
with(summary(Ni.full.model), 1 - deviance/null.deviance)

Ni.step.model <- Ni.full.model %>% stepAIC(trace = T)
coef(Ni.step.model)
Anova(Ni.step.model)
summary(Ni.step.model)
with(summary(Ni.step.model), 1 - deviance/null.deviance)

Ni.final.model <- glm(Ni ~ (Ni_TRD+ OM + TURNOVER)^2, data = Ni_DF, family = gaussian)
coef(Ni.final.model)
Anova(Ni.final.model, test="LR")
with(summary(Ni.final.model), 1 - deviance/null.deviance)

plot(allEffects(Ni.final.model))

plot(Ni_DF$Ni_TRD, Ni_DF$Ni)

plot(Ni_DF$Ni_COL, Ni_DF$Ni)

plot(Ni_DF$Ni_SPM, Ni_DF$Ni)

plot(Ni_DF$TURNOVER, Ni_DF$Ni)

plot(Ni_DF$OM, Ni_DF$Ni)

calc.relimp(Ni.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Ni_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Ni","colloidal_Ni","trulydissolved_Ni","OM.AREA.g.m2","TURNOVER","Ni")])

names(Ni_DF)<-c("COMP","Ni_SPM","Ni_COL","Ni_TRD","OM","TURNOVER", "Ni")

sum(Ni_DF$Ni_SPM, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_COL, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_TRD, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))

Ni_DF<-Ni_DF[Ni_DF$OM>0,]

#Ni_DF$OM <- log10(Ni_DF$OM+1)

Ni_DF$TURNOVER <- sqrt(Ni_DF$TURNOVER)

Ni_DF <- Ni_DF[complete.cases(Ni_DF),]

max.VIF<-5

Ni.RED <- data.frame(scale(Ni_DF[, c(2:(ncol(Ni_DF)-1))],center=T, scale=T))

Ni.RED <- replace_outliers_with_na(Ni.RED)

VIF.Ni<-vif(Ni.RED)

repeat {   
  if (max(VIF.Ni[,2])>max.VIF) {
    VIF.Ni<-vif(Ni.RED)
    Ni.RED <- within(Ni.RED, rm(list=VIF.Ni[which.max(VIF.Ni[,2]),1]))
  } else {
    break
  }
}

Ni_DF<-data.frame(Ni_DF[,1],Ni.RED,Ni_DF[,ncol(Ni_DF)])

Ni_DF <- Ni_DF[complete.cases(Ni_DF),]

names(Ni_DF)[1]<-"COMP"

names(Ni_DF)[ncol(Ni_DF)]<-"Ni"

Ni_DF <- Ni_DF[,-1]

Ni.full.model <- glm(Ni ~. ^2, data = Ni_DF, family = gaussian)
coef(Ni.full.model)
Anova(Ni.full.model, test="LR")
with(summary(Ni.full.model), 1 - deviance/null.deviance)

Ni.step.model <- Ni.full.model %>% stepAIC(trace = T)
coef(Ni.step.model)
Anova(Ni.step.model)
summary(Ni.step.model)
with(summary(Ni.step.model), 1 - deviance/null.deviance)

Ni.final.model <- glm(Ni ~ (OM)^2, data = Ni_DF, family = gaussian)
coef(Ni.final.model)
Anova(Ni.final.model, test="LR")
summary(Ni.final.model)
with(summary(Ni.final.model), 1 - deviance/null.deviance)

plot(allEffects(Ni.final.model))

plot(Ni_DF$Ni_TRD, Ni_DF$Ni)

plot(Ni_DF$Ni_COL, Ni_DF$Ni)

plot(Ni_DF$Ni_SPM, Ni_DF$Ni)

plot(Ni_DF$TURNOVER, Ni_DF$Ni)

plot(Ni_DF$OM, Ni_DF$Ni)

calc.relimp(Ni.final.model)

