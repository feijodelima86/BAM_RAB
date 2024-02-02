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
library(MuMIn)


## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

replace_outliers_with_na <- function(x, threshold = 3) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}


#### Cd ####

Cd_DF <- cbind(COMPARTMENTS_AVG[,c("SITE","FIELD.REP","spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","TURNOVER","Cd"),],(paste(COMPARTMENTS_AVG$SAMPLING_DATE ,COMPARTMENTS_AVG$SITE)))

Cd_DF

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

### Cd_TRD ###

Cd_TRD <- (Cd_DF[,c(1,2,5,6,7,8,9)])

names(Cd_TRD)<- c("SITE", "REP", "Cd_TRD", "TURNOVER", "Cd","DATESITE", "cat")

Cd_TRD<-Cd_TRD[Cd_TRD$Cd_TRD>0,]

#Cd_TRD$Cd_TRD <- log10(Cd_TRD$Cd_TRD+1)

for (col in colnames(Cd_TRD[,c(1)])) {
  Cd_TRD[[col]] <- replace_outliers_with_na(Cd_TRD[[col]])
}

Cd_TRD <- Cd_TRD[complete.cases(Cd_TRD), ]

LMERFIT <- glmer(formula = Cd ~ 1 + Cd_TRD+ (Cd_TRD|SITE) ,
                 data = Cd_TRD,
                 family=gaussian(link = "log"),
                 nAGQ = 0
                 #                        nAGQ=25
)

par(mfrow=c(1,2))

plot(glmm.int_Cd_TRD)

plot(fitted(glmm.int_Cd_TRD),resid(glmm.int_Cd_TRD))

abline(h=0,lty=2,col="red")

qqnorm(resid(glmm.int_Cd_TRD))
qqline(resid(glmm.int_Cd_TRD))

summary(glmm.int_Cd_TRD)

r.squaredGLMM(glmm.int_Cd_TRD)

Anova(glmm.int_Cd_TRD, test="Chisq")

RIaS <-unlist( ranef(LMERFIT)) #Random Intercepts and Slopes
FixedEff <- fixef(LMERFIT)    # Fixed Intercept and Slope

par(mfrow=c(1,2))

plot(Cd_TRD$Cd_TRD,Cd_TRD$Cd,xlab="x",ylab="readings")

#, xlim=c(0,3), ylim=c(min(Y)-1,max(Y)+1)

plot  (Cd_TRD[Cd_TRD$SITE == 'WS', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'WS', ]$Cd, pch=16,xlab="x",ylab="readings",
       ylim = c(0, max(Cd_TRD$Cd,na.rm=TRUE)),
       xlim = c(min(Cd_TRD$Cd_TRD,na.rm=TRUE), max(Cd_TRD$Cd_TRD,na.rm=TRUE)),
       )

points(Cd_TRD[Cd_TRD$SITE == 'DL', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'DL', ]$Cd, col='red', pch=16)
points(Cd_TRD[Cd_TRD$SITE == 'GR', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'GR', ]$Cd, col='green', pch=16)
points(Cd_TRD[Cd_TRD$SITE == 'GC', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'GC', ]$Cd, col='blue', pch=16)
points(Cd_TRD[Cd_TRD$SITE == 'BG', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'BG', ]$Cd, col='orange', pch=16)
points(Cd_TRD[Cd_TRD$SITE == 'BN', ]$Cd_TRD,Cd_TRD[Cd_TRD$SITE == 'BN', ]$Cd, col='brown', pch=16)

lines(sort(Cd_TRD[Cd_TRD$SITE == 'WS', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[7]  + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'WS', ]$Cd_TRD) + RIaS[1]), col='black' , lwd=4)
lines(sort(Cd_TRD[Cd_TRD$SITE == 'DL', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[8]  + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'DL', ]$Cd_TRD) + RIaS[2]), col='red'   , lwd=4)
lines(sort(Cd_TRD[Cd_TRD$SITE == 'GR', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[9]  + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'GR', ]$Cd_TRD) + RIaS[3]), col='green' , lwd=4)
lines(sort(Cd_TRD[Cd_TRD$SITE == 'GC', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[10] + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'GC', ]$Cd_TRD) + RIaS[4]), col='blue'  , lwd=4)
lines(sort(Cd_TRD[Cd_TRD$SITE == 'BG', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[11] + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'BG', ]$Cd_TRD) + RIaS[5]), col='orange', lwd=4) 
lines(sort(Cd_TRD[Cd_TRD$SITE == 'BN', ]$Cd_TRD), exp(FixedEff[1]+ (RIaS[12] + FixedEff[2]) * sort(Cd_TRD[Cd_TRD$SITE == 'BN', ]$Cd_TRD) + RIaS[6]), col='brown' , lwd=4) 

#abline(v=(seq(-1,4 ,1)), col="lightgray", lty="dotted");        
#abline(h=(seq( -1,25 ,1)), col="lightgray", lty="dotted")   

