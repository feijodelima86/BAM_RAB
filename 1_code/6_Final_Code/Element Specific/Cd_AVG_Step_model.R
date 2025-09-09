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

#### Cd ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

names(COMPARTMENTS_AVG)

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Cd_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","OM.AREA.g.m2","TURNOVER","Cd")])

names(Cd_DF)<-c("COMP","Cd_SPM","Cd_COL","Cd_TRD","OM","TURNOVER", "Cd")

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

#Cd_DF$OM <- log10(Cd_DF$OM+1)

#Cd_DF$TURNOVER <- sqrt(Cd_DF$TURNOVER)

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

Cd_DF$OM<-log10(Cd_DF$OM+1)

Cd_DF$TURNOVER <- sqrt(Cd_DF$TURNOVER)

max.VIF<-5

Cd.RED <- data.frame(scale(Cd_DF[, c(2:(ncol(Cd_DF)-1))],center=T, scale=T))

Cd.RED <- replace_outliers_with_na(Cd.RED)

VIF.Cd<-vif(Cd.RED)

repeat {   
  if (max(VIF.Cd[,2])>max.VIF) {
    VIF.Cd<-vif(Cd.RED)
    Cd.RED <- within(Cd.RED, rm(list=VIF.Cd[which.max(VIF.Cd[,2]),1]))
  } else {
    break
  }
}

Cd_DF<-data.frame(Cd_DF[,1],Cd.RED,Cd_DF[,ncol(Cd_DF)])

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

names(Cd_DF)[1]<-"COMP"

names(Cd_DF)[ncol(Cd_DF)]<-"Cd"

cd.full.model <- glm(Cd ~. ^2, data = Cd_DF, family = gaussian)
coef(cd.full.model)
Anova(cd.full.model, test="LR")
with(summary(cd.full.model), 1 - deviance/null.deviance)

cd.step.model <- cd.full.model %>% stepAIC(trace = T)
coef(cd.step.model)
Anova(cd.step.model)
summary(cd.step.model)
with(summary(cd.step.model), 1 - deviance/null.deviance)

cd.final.model <- glm(Cd ~ (Cd_TRD), data = Cd_DF, family = gaussian)
coef(cd.final.model)
Anova(cd.final.model, test="LR")
with(summary(cd.final.model), 1 - deviance/null.deviance)

plot(allEffects(cd.final.model))

plot(Cd_DF$Cd_TRD, Cd_DF$Cd)

plot(Cd_DF$Cd_COL, Cd_DF$Cd)

plot(Cd_DF$Cd_SPM, Cd_DF$Cd)

plot(Cd_DF$TURNOVER, Cd_DF$Cd)

plot(Cd_DF$OM, Cd_DF$Cd)

calc.relimp(cd.final.model)

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Rank variable for colour assignment
Cd_DF$order = findInterval(Cd_DF$TURNOVER, sort(Cd_DF$TURNOVER))
# Make plot
plot(Cd ~ Cd_TRD, Cd_DF, pch=19, col=pal(nrow(Cd_DF))[Cd_DF$order])
# Add a simple legend
legend("topright", col=pal(2), pch=19,
       legend=c(round(range(df$Temp), 1)))



#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

Cd_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","OM.AREA.g.m2","TURNOVER","Cd")])

names(Cd_DF)<-c("COMP","Cd_SPM","Cd_COL","Cd_TRD","OM","TURNOVER", "Cd")

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

Cd_DF<-Cd_DF[Cd_DF$OM>0,]

#Cd_DF$OM <- log10(Cd_DF$OM+1)

Cd_DF$TURNOVER <- sqrt(Cd_DF$TURNOVER)

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

max.VIF<-5

Cd.RED <- data.frame(scale(Cd_DF[, c(2:(ncol(Cd_DF)-1))],center=T, scale=T))

Cd.RED <- replace_outliers_with_na(Cd.RED)

VIF.Cd<-vif(Cd.RED)

repeat {   
  if (max(VIF.Cd[,2])>max.VIF) {
    VIF.Cd<-vif(Cd.RED)
    Cd.RED <- within(Cd.RED, rm(list=VIF.Cd[which.max(VIF.Cd[,2]),1]))
  } else {
    break
  }
}

Cd_DF<-data.frame(Cd_DF[,1],Cd.RED,Cd_DF[,ncol(Cd_DF)])

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

names(Cd_DF)[1]<-"COMP"

names(Cd_DF)[ncol(Cd_DF)]<-"Cd"

Cd_DF <- Cd_DF[,-1]

Cd_DF <- Cd_DF[,-2]

cd.full.model <- glm(Cd ~. ^2, data = Cd_DF, family = gaussian)
coef(cd.full.model)
Anova(cd.full.model, test="LR")
with(summary(cd.full.model), 1 - deviance/null.deviance)

cd.step.model <- cd.full.model %>% stepAIC(trace = T)
coef(cd.step.model)
Anova(cd.step.model)
summary(cd.step.model)
with(summary(cd.step.model), 1 - deviance/null.deviance)

cd.final.model <- glm(Cd ~ (OM), data = Cd_DF, family = gaussian)
coef(cd.final.model)
Anova(cd.final.model, test="LR")
summary(cd.final.model)
with(summary(cd.final.model), 1 - deviance/null.deviance)

plot(allEffects(cd.final.model))

plot(Cd_DF$Cd_TRD, Cd_DF$Cd)

plot(Cd_DF$Cd_COL, Cd_DF$Cd)

plot(Cd_DF$Cd_SPM, Cd_DF$Cd)

plot(Cd_DF$TURNOVER, Cd_DF$Cd)

plot(Cd_DF$OM, Cd_DF$Cd)

calc.relimp(cd.final.model)

