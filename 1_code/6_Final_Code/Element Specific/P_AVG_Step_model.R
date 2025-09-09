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

#### P ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

P_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_P","colloidal_P","trulydissolved_P","OM.AREA.g.m2","TURNOVER","P")])

names(P_DF)<-c("COMP","P_SPM","P_COL","P_TRD","OM","TURNOVER", "P")

sum(P_DF$P_SPM, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_COL, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_TRD, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))

P_DF<-P_DF[P_DF$OM>0,]

#P_DF$OM <- log10(P_DF$OM+1)

P_DF$TURNOVER <- sqrt(P_DF$TURNOVER)

P_DF <- P_DF[complete.cases(P_DF),]

max.VIF<-5

P.RED <- data.frame(scale(P_DF[, c(2:(ncol(P_DF)-1))],center=T, scale=T))

P.RED <- replace_outliers_with_na(P.RED)

VIF.P<-vif(P.RED)

repeat {   
  if (max(VIF.P[,2])>max.VIF) {
    VIF.P<-vif(P.RED)
    P.RED <- within(P.RED, rm(list=VIF.P[which.max(VIF.P[,2]),1]))
  } else {
    break
  }
}

P_DF<-data.frame(P_DF[,1],P.RED,P_DF[,ncol(P_DF)])

P_DF <- P_DF[complete.cases(P_DF),]

names(P_DF)[1]<-"COMP"

names(P_DF)[ncol(P_DF)]<-"P"

P.full.model <- glm(P ~. ^2, data = P_DF, family = gaussian)
coef(P.full.model)
Anova(P.full.model, test="LR")
with(summary(P.full.model), 1 - deviance/null.deviance)

P.step.model <- P.full.model %>% stepAIC(trace = T)
coef(P.step.model)
Anova(P.step.model)
summary(P.step.model)
with(summary(P.step.model), 1 - deviance/null.deviance)

P.final.model <- glm(P ~ (P_TRD + COMP + OM), data = P_DF, family = gaussian)
coef(P.final.model)
summary(P.final.model)
Anova(P.final.model, test="LR")
with(summary(P.final.model), 1 - deviance/null.deviance)

plot(allEffects(P.final.model))

plot(P_DF$P_TRD, P_DF$P)

plot(P_DF$P_COL, P_DF$P)

plot(P_DF$P_SPM, P_DF$P)

plot(P_DF$TURNOVER, P_DF$P)

plot(P_DF$OM, P_DF$P)

calc.relimp(P.final.model)

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Rank variable for colour assignment
P_DF$order = findInterval(P_DF$COMP, sort(P_DF$COMP))
# Make plot
plot(P ~ P_TRD, P_DF, pch=19, col=pal(nrow(P_DF))[P_DF$order])
# Add a simple legend
legend("topright", col=pal(2), pch=19,
       legend=c(round(range(df$Temp), 1)))



#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

P_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_P","colloidal_P","trulydissolved_P","OM.AREA.g.m2","TURNOVER","P")])

names(P_DF)<-c("COMP","P_SPM","P_COL","P_TRD","OM","TURNOVER", "P")

sum(P_DF$P_SPM, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_COL, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_TRD, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))

P_DF<-P_DF[P_DF$OM>0,]

#P_DF$OM <- log10(P_DF$OM+1)

P_DF$TURNOVER <- sqrt(P_DF$TURNOVER)

P_DF <- P_DF[complete.cases(P_DF),]

max.VIF<-5

P.RED <- data.frame(scale(P_DF[, c(2:(ncol(P_DF)-1))],center=T, scale=T))

P.RED <- replace_outliers_with_na(P.RED)

VIF.P<-vif(P.RED)

repeat {   
  if (max(VIF.P[,2])>max.VIF) {
    VIF.P<-vif(P.RED)
    P.RED <- within(P.RED, rm(list=VIF.P[which.max(VIF.P[,2]),1]))
  } else {
    break
  }
}

P_DF<-data.frame(P_DF[,1],P.RED,P_DF[,ncol(P_DF)])

P_DF <- P_DF[complete.cases(P_DF),]

names(P_DF)[1]<-"COMP"

names(P_DF)[ncol(P_DF)]<-"P"

P_DF <- P_DF[,-1]

P.full.model <- glm(P ~. ^2, data = P_DF, family = gaussian)
coef(P.full.model)
Anova(P.full.model, test="LR")
with(summary(P.full.model), 1 - deviance/null.deviance)

P.step.model <- P.full.model %>% stepAIC(trace = T)
coef(P.step.model)
Anova(P.step.model)
summary(P.step.model)
with(summary(P.step.model), 1 - deviance/null.deviance)

P.final.model <- glm(P ~ (P_SPM +TURNOVER)^2, data = P_DF, family = gaussian)
coef(P.final.model)
Anova(P.final.model, test="LR")
summary(P.final.model)
with(summary(P.final.model), 1 - deviance/null.deviance)

plot(allEffects(P.final.model))

plot(P_DF$P_TRD, P_DF$P)

plot(P_DF$P_COL, P_DF$P)

plot(P_DF$P_SPM, P_DF$P)

plot(P_DF$TURNOVER, P_DF$P)

plot(P_DF$OM, P_DF$P)

calc.relimp(P.final.model)
