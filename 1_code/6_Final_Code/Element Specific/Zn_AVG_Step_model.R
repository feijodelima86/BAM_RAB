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

#### Zn ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Zn_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Zn","colloidal_Zn","trulydissolved_Zn","OM.AREA.g.m2","TURNOVER","Zn")])

names(Zn_DF)<-c("COMP","Zn_SPM","Zn_COL","Zn_TRD","OM","TURNOVER", "Zn")

sum(Zn_DF$Zn_SPM, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_COL, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_TRD, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))

Zn_DF<-Zn_DF[Zn_DF$OM>0,]

#Zn_DF$OM <- log10(Zn_DF$OM+1)

Zn_DF$TURNOVER <- sqrt(Zn_DF$TURNOVER)

Zn_DF <- Zn_DF[complete.cases(Zn_DF),]

max.VIF<-5

Zn.RED <- data.frame(scale(Zn_DF[, c(2:(ncol(Zn_DF)-1))],center=T, scale=T))

Zn.RED <- replace_outliers_with_na(Zn.RED)

VIF.Zn<-vif(Zn.RED)

repeat {   
  if (max(VIF.Zn[,2])>max.VIF) {
    VIF.Zn<-vif(Zn.RED)
    Zn.RED <- within(Zn.RED, rm(list=VIF.Zn[which.max(VIF.Zn[,2]),1]))
  } else {
    break
  }
}

Zn_DF<-data.frame(Zn_DF[,1],Zn.RED,Zn_DF[,ncol(Zn_DF)])

Zn_DF <- Zn_DF[complete.cases(Zn_DF),]

names(Zn_DF)[1]<-"COMP"

names(Zn_DF)[ncol(Zn_DF)]<-"Zn"

Zn.full.model <- glm(Zn ~. ^2, data = Zn_DF, family = gaussian)
coef(Zn.full.model)
Anova(Zn.full.model, test="LR")
with(summary(Zn.full.model), 1 - deviance/null.deviance)

Zn.step.model <- Zn.full.model %>% stepAIC(trace = T)
coef(Zn.step.model)
Anova(Zn.step.model)
summary(Zn.step.model)
with(summary(Zn.step.model), 1 - deviance/null.deviance)

Zn.final.model <- glm(Zn ~ (Zn_TRD+ OM), data = Zn_DF, family = gaussian)
coef(Zn.final.model)
Anova(Zn.final.model, test="LR")
with(summary(Zn.final.model), 1 - deviance/null.deviance)

plot(allEffects(Zn.final.model))

plot(Zn_DF$Zn_TRD, Zn_DF$Zn)

plot(Zn_DF$Zn_COL, Zn_DF$Zn)

plot(Zn_DF$Zn_SPM, Zn_DF$Zn)

plot(Zn_DF$TURNOVER, Zn_DF$Zn)

plot(Zn_DF$OM, Zn_DF$Zn)

calc.relimp(Zn.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Zn_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Zn","colloidal_Zn","trulydissolved_Zn","OM.AREA.g.m2","TURNOVER","Zn")])

names(Zn_DF)<-c("COMP","Zn_SPM","Zn_COL","Zn_TRD","OM","TURNOVER", "Zn")

sum(Zn_DF$Zn_SPM, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_COL, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_TRD, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))

Zn_DF<-Zn_DF[Zn_DF$OM>0,]

#Zn_DF$OM <- log10(Zn_DF$OM+1)

Zn_DF$TURNOVER <- sqrt(Zn_DF$TURNOVER)

Zn_DF <- Zn_DF[complete.cases(Zn_DF),]

max.VIF<-5

Zn.RED <- data.frame(scale(Zn_DF[, c(2:(ncol(Zn_DF)-1))],center=T, scale=T))

Zn.RED <- replace_outliers_with_na(Zn.RED)

VIF.Zn<-vif(Zn.RED)

repeat {   
  if (max(VIF.Zn[,2])>max.VIF) {
    VIF.Zn<-vif(Zn.RED)
    Zn.RED <- within(Zn.RED, rm(list=VIF.Zn[which.max(VIF.Zn[,2]),1]))
  } else {
    break
  }
}

Zn_DF<-data.frame(Zn_DF[,1],Zn.RED,Zn_DF[,ncol(Zn_DF)])

Zn_DF <- Zn_DF[complete.cases(Zn_DF),]

names(Zn_DF)[1]<-"COMP"

names(Zn_DF)[ncol(Zn_DF)]<-"Zn"

Zn_DF <- Zn_DF[,-1]

Zn.full.model <- glm(Zn ~. ^2, data = Zn_DF, family = gaussian)
coef(Zn.full.model)
Anova(Zn.full.model, test="LR")
with(summary(Zn.full.model), 1 - deviance/null.deviance)

Zn.step.model <- Zn.full.model %>% stepAIC(trace = T)
coef(Zn.step.model)
Anova(Zn.step.model)
summary(Zn.step.model)
with(summary(Zn.step.model), 1 - deviance/null.deviance)

Zn.final.model <- glm(Zn ~ (OM)^2, data = Zn_DF, family = gaussian)
coef(Zn.final.model)
Anova(Zn.final.model, test="LR")
summary(Zn.final.model)
with(summary(Zn.final.model), 1 - deviance/null.deviance)

plot(allEffects(Zn.final.model))

plot(Zn_DF$Zn_TRD, Zn_DF$Zn)

plot(Zn_DF$Zn_COL, Zn_DF$Zn)

plot(Zn_DF$Zn_SPM, Zn_DF$Zn)

plot(Zn_DF$TURNOVER, Zn_DF$Zn)

plot(Zn_DF$OM, Zn_DF$Zn)

calc.relimp(Zn.final.model)

