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
library(MCdS)
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

#### Cd ####

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

### Epis ##

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

Cd_DF <- (COMPARTMENTS_AVG[,c("SAMPLING_DATE","SAMPLE_DESCRIPTOR","spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","OM.AREA.g.m2","TURNOVER","Cd")])

names(Cd_DF)<-c("DATE","COMP","Cd_SPM","Cd_COL","Cd_TRD","OM","TURNOVER", "Cd")

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

Cd_DF$OM<-log10(Cd_DF$OM+1)

Cd_DF$TURNOVER <- sqrt(Cd_DF$TURNOVER)

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

Cd_DF$DATE<-yday(anydate(Cd_DF$DATE))

max.VIF<-5

Cd.RED <- data.frame(scale(Cd_DF[, c(3:(ncol(Cd_DF)-1))],center=T, scale=T))

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

Cd_DF<-data.frame(Cd_DF[,1],Cd_DF[,2],Cd.RED,Cd_DF[,ncol(Cd_DF)])

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

names(Cd_DF)[1]<-"DATE"
names(Cd_DF)[2]<-"COMP"

names(Cd_DF)[ncol(Cd_DF)]<-"Cd"

Cd.full.model <- glm(Cd ~. ^2, data = Cd_DF, family = gaussian)
coef(Cd.full.model)
summary(Cd.full.model)
Anova(Cd.full.model, test="LR")
with(summary(Cd.full.model), 1 - deviance/null.deviance)

Cd.step.model <- Cd.full.model %>% stepAIC(trace = T)
coef(Cd.step.model)
Anova(Cd.step.model)
summary(Cd.step.model)
with(summary(Cd.step.model), 1 - deviance/null.deviance)

Cd.final.model <- glm(Cd ~ Cd_TRD + DATE, data = Cd_DF, family = gaussian)
coef(Cd.final.model)
Anova(Cd.final.model)
with(summary(Cd.final.model), 1 - deviance/null.deviance)

plot(allEffects(Cd.final.model))

plot(Cd_DF$Cd_TRD, Cd_DF$Cd)

plot(Cd_DF$Cd_COL, Cd_DF$Cd)

plot(Cd_DF$Cd_SPM, Cd_DF$Cd)

plot(Cd_DF$TURNOVER, Cd_DF$Cd)

plot(Cd_DF$OM, Cd_DF$Cd)

calc.relimp(Cd.final.model)

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Rank variable for colour Cdsignment
Cd_DF$order = findInterval(Cd_DF$DATE, sort(Cd_DF$DATE))
# Make plot
plot(Cd ~ Cd_TRD, Cd_DF, pch=19, col=pal(nrow(Cd_DF))[Cd_DF$order])
# Add a simple legend
legend("topright", col=pal(2), pch=19,
       legend=c(round(range(df$Temp), 1)))




