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

#### Pb ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Pb_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Pb_NH3","colloidal_Pb_NH3","trulydissolved_Pb_NH3","OM.AREA.g.m2","TURNOVER","Pb")])

names(Pb_DF)<-c("COMP","Pb_SPM","Pb_COL","Pb_TRD","OM","TURNOVER", "Pb")

sum(Pb_DF$Pb_SPM, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_COL, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_TRD, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))

Pb_DF<-Pb_DF[Pb_DF$OM>0,]

#Pb_DF$OM <- log10(Pb_DF$OM+1)

Pb_DF$TURNOVER <- sqrt(Pb_DF$TURNOVER)

Pb_DF <- Pb_DF[complete.cases(Pb_DF),]

max.VIF<-5

Pb.RED <- data.frame(scale(Pb_DF[, c(2:(ncol(Pb_DF)-1))],center=T, scale=T))

Pb.RED <- replace_outliers_with_na(Pb.RED)

VIF.Pb<-vif(Pb.RED)

repeat {   
  if (max(VIF.Pb[,2])>max.VIF) {
    VIF.Pb<-vif(Pb.RED)
    Pb.RED <- within(Pb.RED, rm(list=VIF.Pb[which.max(VIF.Pb[,2]),1]))
  } else {
    break
  }
}

Pb_DF<-data.frame(Pb_DF[,1],Pb.RED,Pb_DF[,ncol(Pb_DF)])

Pb_DF <- Pb_DF[complete.cases(Pb_DF),]

names(Pb_DF)[1]<-"COMP"

names(Pb_DF)[ncol(Pb_DF)]<-"Pb"

Pb.full.model <- glm(Pb ~. ^2, data = Pb_DF, family = gaussian)
coef(Pb.full.model)
Anova(Pb.full.model, test="LR")
with(summary(Pb.full.model), 1 - deviance/null.deviance)

Pb.step.model <- Pb.full.model %>% stepAIC(trace = T)
coef(Pb.step.model)
Anova(Pb.step.model)
summary(Pb.step.model)
with(summary(Pb.step.model), 1 - deviance/null.deviance)

Pb.final.model <- glm(Pb ~ (Pb_TRD+OM+TURNOVER)^2, data = Pb_DF, family = gaussian)
coef(Pb.final.model)
Anova(Pb.final.model, test="LR")
summary(Pb.final.model)
with(summary(Pb.final.model), 1 - deviance/null.deviance)

plot(allEffects(Pb.final.model))

plot(Pb_DF$Pb_TRD, Pb_DF$Pb)

plot(Pb_DF$Pb_COL, Pb_DF$Pb)

plot(Pb_DF$Pb_SPM, Pb_DF$Pb)

plot(Pb_DF$TURNOVER, Pb_DF$Pb)

plot(Pb_DF$OM, Pb_DF$Pb)

calc.relimp(Pb.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Pb_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Pb_NH3","colloidal_Pb_NH3","trulydissolved_Pb_NH3","OM.AREA.g.m2","TURNOVER","Pb")])

names(Pb_DF)<-c("COMP","Pb_SPM","Pb_COL","Pb_TRD","OM","TURNOVER", "Pb")

sum(Pb_DF$Pb_SPM, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_COL, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_TRD, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))

Pb_DF<-Pb_DF[Pb_DF$OM>0,]

#Pb_DF$OM <- log10(Pb_DF$OM+1)

Pb_DF$TURNOVER <- sqrt(Pb_DF$TURNOVER)

Pb_DF <- Pb_DF[complete.cases(Pb_DF),]

max.VIF<-5

Pb.RED <- data.frame(scale(Pb_DF[, c(2:(ncol(Pb_DF)-1))],center=T, scale=T))

Pb.RED <- replace_outliers_with_na(Pb.RED)

VIF.Pb<-vif(Pb.RED)

repeat {   
  if (max(VIF.Pb[,2])>max.VIF) {
    VIF.Pb<-vif(Pb.RED)
    Pb.RED <- within(Pb.RED, rm(list=VIF.Pb[which.max(VIF.Pb[,2]),1]))
  } else {
    break
  }
}

Pb_DF<-data.frame(Pb_DF[,1],Pb.RED,Pb_DF[,ncol(Pb_DF)])

Pb_DF <- Pb_DF[complete.cases(Pb_DF),]

names(Pb_DF)[1]<-"COMP"

names(Pb_DF)[ncol(Pb_DF)]<-"Pb"

Pb_DF <- Pb_DF[,-1]

Pb_DF <- Pb_DF[,-2]

Pb.full.model <- glm(Pb ~. ^2, data = Pb_DF, family = gaussian)
coef(Pb.full.model)
Anova(Pb.full.model, test="LR")
with(summary(Pb.full.model), 1 - deviance/null.deviance)

Pb.step.model <- Pb.full.model %>% stepAIC(trace = T)
coef(Pb.step.model)
Anova(Pb.step.model)
summary(Pb.step.model)
with(summary(Pb.step.model), 1 - deviance/null.deviance)

Pb.final.model <- glm(Pb ~ (Pb_TRD + OM + TURNOVER + OM:TURNOVER), data = Pb_DF, family = gaussian)
coef(Pb.final.model)
Anova(Pb.final.model, test="LR")
summary(Pb.final.model)
with(summary(Pb.final.model), 1 - deviance/null.deviance)

plot(allEffects(Pb.final.model))

plot(Pb_DF$Pb_TRD, Pb_DF$Pb)

plot(Pb_DF$Pb_COL, Pb_DF$Pb)

plot(Pb_DF$Pb_SPM, Pb_DF$Pb)

plot(Pb_DF$TURNOVER, Pb_DF$Pb)

plot(Pb_DF$OM, Pb_DF$Pb)

calc.relimp(Pb.final.model)

