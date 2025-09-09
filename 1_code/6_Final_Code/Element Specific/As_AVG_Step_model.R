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
library(lubridate)
library(anytime)

#### Housekeeping

replace_outliers_with_na <- function(x, threshold = 3) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}

## Data Onborading ###

#### As ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

As_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","OM.AREA.g.m2","TURNOVER","trulydissolved_P","As")])

names(As_DF)<-c("COMP","As_SPM","As_COL","As_TRD","OM","TURNOVER","P","As")

sum(As_DF$As_SPM, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_COL, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_TRD, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))

As_DF$OM<-log10(As_DF$OM+1)

As_DF$TURNOVER <- sqrt(As_DF$TURNOVER)

As_DF <- As_DF[complete.cases(As_DF),]

max.VIF<-5

As.RED <- data.frame(scale(As_DF[, c(2:(ncol(As_DF)-1))],center=T, scale=T))

As.RED <- replace_outliers_with_na(As.RED)

VIF.As<-vif(As.RED)

repeat {   
  if (max(VIF.As[,2])>max.VIF) {
    VIF.As<-vif(As.RED)
    As.RED <- within(As.RED, rm(list=VIF.As[which.max(VIF.As[,2]),1]))
  } else {
    break
  }
}

As_DF<-data.frame(As_DF[,1],As.RED,As_DF[,ncol(As_DF)])

As_DF <- As_DF[complete.cases(As_DF),]

names(As_DF)[1]<-"COMP"

names(As_DF)[ncol(As_DF)]<-"As"

As.full.model <- glm(As ~. ^2, data = As_DF, family = gaussian)
coef(As.full.model)
summary(As.full.model)
Anova(As.full.model, test="LR")
with(summary(As.full.model), 1 - deviance/null.deviance)

As.step.model <- As.full.model %>% stepAIC(trace = T)
coef(As.step.model)
Anova(As.step.model)
summary(As.step.model)
with(summary(As.step.model), 1 - deviance/null.deviance)

As.final.model <- glm(As ~ (COMP + As_TRD + P)^2, data = As_DF, family = gaussian)
coef(As.final.model)
Anova(As.final.model)
with(summary(As.final.model), 1 - deviance/null.deviance)

plot(allEffects(As.final.model))

plot(As_DF$As_TRD, As_DF$As)

plot(As_DF$P, As_DF$As)

plot(As_DF$As_SPM, As_DF$As)

plot(As_DF$TURNOVER, As_DF$As)

plot(As_DF$OM, As_DF$As)

calc.relimp(As.final.model)

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Rank variable for colour assignment
As_DF$order = findInterval(As_DF$P, sort(As_DF$P))
# Make plot
plot(As ~ As_TRD, As_DF, pch=19, col=pal(nrow(As_DF))[As_DF$order])
# Add a simple legend
legend("topright", col=pal(2), pch=19,
       legend=c(round(range(df$Temp), 1)))

#### Fila

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

As_DF <- (COMPARTMENTS_AVG[,c("spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","OM.AREA.g.m2","TURNOVER","As")])

names(As_DF)<-c("As_SPM","As_COL","As_TRD","OM","TURNOVER", "As")

sum(As_DF$As_SPM, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_COL, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_TRD, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))

As_DF$OM<-log10(As_DF$OM+1)

As_DF$TURNOVER <- sqrt(As_DF$TURNOVER)

As_DF <- As_DF[complete.cases(As_DF),]

As_DF$DATE<-yday(anydate(As_DF$DATE))

max.VIF<-5

As.RED <- data.frame(scale(As_DF[, c(2:(ncol(As_DF)-1))],center=T, scale=T))

As.RED <- replace_outliers_with_na(As.RED)

VIF.As<-vif(As.RED)

repeat {   
  if (max(VIF.As[,2])>max.VIF) {
    VIF.As<-vif(As.RED)
    As.RED <- within(As.RED, rm(list=VIF.As[which.max(VIF.As[,2]),1]))
  } else {
    break
  }
}

As_DF<-data.frame(As_DF[,1],As.RED,As_DF[,ncol(As_DF)])

As_DF <- As_DF[complete.cases(As_DF),]

names(As_DF)[1]<-"DATE"

names(As_DF)[ncol(As_DF)]<-"As"

As.full.model <- glm(As ~. , data = As_DF, family = gaussian)
coef(As.full.model)
summary(As.full.model)
Anova(As.full.model, test="LR")
with(summary(As.full.model), 1 - deviance/null.deviance)

As.step.model <- As.full.model %>% stepAIC(trace = T)
coef(As.step.model)
Anova(As.step.model)
summary(As.step.model)
with(summary(As.step.model), 1 - deviance/null.deviance)

As.final.model <- glm(As ~ As_TRD * DATE, data = As_DF, family = gaussian)
coef(As.final.model)
Anova(As.final.model)
summary(As.final.model)
with(summary(As.final.model), 1 - deviance/null.deviance)

plot(allEffects(As.final.model))

plot(As_DF$As_TRD, As_DF$As)

plot(As_DF$As_COL, As_DF$As)

plot(As_DF$As_SPM, As_DF$As)

plot(As_DF$TURNOVER, As_DF$As)

plot(As_DF$OM, As_DF$As)

calc.relimp(As.final.model)

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Rank variable for colour assignment
As_DF$order = findInterval(As_DF$COMP, sort(As_DF$COMP))
# Make plot
plot(As ~ As_TRD, As_DF, pch=19, col=pal(nrow(As_DF))[As_DF$order])




