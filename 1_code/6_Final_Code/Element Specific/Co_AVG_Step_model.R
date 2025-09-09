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

#### Co ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Co_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Co","F_Co","UF_Co","OM.AREA.g.m2","TURNOVER","Co")])

names(Co_DF)<-c("COMP","Co_SPM","Co_COL","Co_TRD","OM","TURNOVER", "Co")

sum(Co_DF$Co_SPM, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_COL, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_TRD, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))

Co_DF<-Co_DF[Co_DF$OM>0,]

#Co_DF$OM <- log10(Co_DF$OM+1)

Co_DF$TURNOVER <- sqrt(Co_DF$TURNOVER)

Co_DF <- Co_DF[complete.cases(Co_DF),]

max.VIF<-5

Co.RED <- data.frame(scale(Co_DF[, c(2:(ncol(Co_DF)-1))],center=T, scale=T))

Co.RED <- replace_outliers_with_na(Co.RED)

VIF.Co<-vif(Co.RED)

repeat {   
  if (max(VIF.Co[,2])>max.VIF) {
    VIF.Co<-vif(Co.RED)
    Co.RED <- within(Co.RED, rm(list=VIF.Co[which.max(VIF.Co[,2]),1]))
  } else {
    break
  }
}

Co_DF<-data.frame(Co_DF[,1],Co.RED,Co_DF[,ncol(Co_DF)])

Co_DF <- Co_DF[complete.cases(Co_DF),]

names(Co_DF)[1]<-"COMP"

names(Co_DF)[ncol(Co_DF)]<-"Co"

Co.full.model <- glm(Co ~. ^2, data = Co_DF, family = gaussian)
coef(Co.full.model)
Anova(Co.full.model, test="LR")
with(summary(Co.full.model), 1 - deviance/null.deviance)

Co.step.model <- Co.full.model %>% stepAIC(trace = T)
coef(Co.step.model)
Anova(Co.step.model)
summary(Co.step.model)
with(summary(Co.step.model), 1 - deviance/null.deviance)

Co.final.model <- glm(Co ~ (Co_COL + OM)^2, data = Co_DF, family = gaussian)
coef(Co.final.model)
Anova(Co.final.model, test="LR")
with(summary(Co.final.model), 1 - deviance/null.deviance)

plot(allEffects(Co.final.model))

plot(Co_DF$Co_TRD, Co_DF$Co)

plot(Co_DF$Co_COL, Co_DF$Co)

plot(Co_DF$Co_SPM, Co_DF$Co)

plot(Co_DF$TURNOVER, Co_DF$Co)

plot(Co_DF$OM, Co_DF$Co)

calc.relimp(Co.final.model)

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Co_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Co","colloidal_Co","trulydissolved_Co","OM.AREA.g.m2","TURNOVER","Co")])

names(Co_DF)<-c("COMP","Co_SPM","Co_COL","Co_TRD","OM","TURNOVER", "Co")

sum(Co_DF$Co_SPM, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_COL, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_TRD, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))

Co_DF<-Co_DF[Co_DF$OM>0,]

#Co_DF$OM <- log10(Co_DF$OM+1)

Co_DF$TURNOVER <- sqrt(Co_DF$TURNOVER)

Co_DF <- Co_DF[complete.cases(Co_DF),]

max.VIF<-5

Co.RED <- data.frame(scale(Co_DF[, c(2:(ncol(Co_DF)-1))],center=T, scale=T))

Co.RED <- replace_outliers_with_na(Co.RED)

VIF.Co<-vif(Co.RED)

repeat {   
  if (max(VIF.Co[,2])>max.VIF) {
    VIF.Co<-vif(Co.RED)
    Co.RED <- within(Co.RED, rm(list=VIF.Co[which.max(VIF.Co[,2]),1]))
  } else {
    break
  }
}

Co_DF<-data.frame(Co_DF[,1],Co.RED,Co_DF[,ncol(Co_DF)])

Co_DF <- Co_DF[complete.cases(Co_DF),]

names(Co_DF)[1]<-"COMP"

names(Co_DF)[ncol(Co_DF)]<-"Co"

Co_DF <- Co_DF[,-1]

#Co_DF <- Co_DF[,-2]

Co.full.model <- glm(Co ~. ^2, data = Co_DF, family = gaussian)
coef(Co.full.model)
Anova(Co.full.model, test="LR")
with(summary(Co.full.model), 1 - deviance/null.deviance)

Co.step.model <- Co.full.model %>% stepAIC(trace = T)
coef(Co.step.model)
Anova(Co.step.model)
summary(Co.step.model)
with(summary(Co.step.model), 1 - deviance/null.deviance)

Co.final.model <- glm(Co ~ (OM)^2, data = Co_DF, family = gaussian)
coef(Co.final.model)
Anova(Co.final.model, test="LR")
summary(Co.final.model)
with(summary(Co.final.model), 1 - deviance/null.deviance)

plot(allEffects(Co.final.model))

plot(Co_DF$Co_TRD, Co_DF$Co)

plot(Co_DF$Co_COL, Co_DF$Co)

plot(Co_DF$Co_SPM, Co_DF$Co)

plot(Co_DF$TURNOVER, Co_DF$Co)

plot(Co_DF$OM, Co_DF$Co)

calc.relimp(Co.final.model)

