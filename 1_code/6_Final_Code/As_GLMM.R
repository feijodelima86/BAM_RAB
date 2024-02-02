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

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

replace_outliers_with_na <- function(x, threshold = 1.5) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}


#### As ####

As_DF <- cbind(COMPARTMENTS_AVG[,c("SITE","FIELD.REP","spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","TURNOVER","As"),],(paste(COMPARTMENTS_AVG$SAMPLING_DATE ,COMPARTMENTS_AVG$SITE)))

As_DF

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

### As_TRD ###

As_TRD <- (As_DF[,c(1,2,5,6,7,8,9)])

names(As_TRD)<- c("SITE", "REP", "As_TRD", "TURNOVER", "As","DATESITE", "cat")

As_TRD<-As_TRD[As_TRD$As_TRD>0,]

As_TRD$As_TRD <- log(As_TRD$As_TRD+1)

#As_TRD$As <- sqrt(As_TRD$As)

for (col in colnames(As_TRD[,c(1)])) {
  As_TRD[[col]] <- replace_outliers_with_na(As_TRD[[col]])
}

As_TRD <- As_TRD[complete.cases(As_TRD), ]

glmm.int_As_TRD <- glmer(formula = As ~ 1 + As_TRD * TURNOVER  + (As_TRD|DATESITE) ,
                 data = As_TRD,
                 family=gaussian(link = "log"),
                 nAGQ = 0
                 #                        nAGQ=25
)

par(mfrow=c(1,2))

plot(glmm.int_As_TRD)

plot(fitted(glmm.int_As_TRD),resid(glmm.int_As_TRD))

abline(h=0,lty=2,col="red")

qqnorm(resid(glmm.int_As_TRD))
qqline(resid(glmm.int_As_TRD))

summary(glmm.int_As_TRD)

r.squaredGLMM(glmm.int_As_TRD)

Anova(glmm.int_As_TRD, test="Chisq")

RIaS <-unlist( ranef(glmm.int_As_TRD)) #Random Intercepts and Slopes
FixedEff <- fixef(glmm.int_As_TRD)    # Fixed Intercept and Slope

par(mfrow=c(1,2))

plot(As_TRD$As_TRD,As_TRD$As,xlab="x",ylab="readings")

#, xlim=c(0,3), ylim=c(min(Y)-1,max(Y)+1)

plot  (As_TRD[As_TRD$cat == '0.21 - 0.46', ]$As_TRD,As_TRD[As_TRD$cat == '0.21 - 0.46', ]$As, pch=16,xlab="x",ylab="readings",
       ylim = c(0, max(As_TRD$As,na.rm=TRUE)),
       xlim = c(min(As_TRD$As_TRD,na.rm=TRUE), max(As_TRD$As_TRD,na.rm=TRUE)),
)

points(As_TRD[As_TRD$cat == '0.46 - 0.50', ]$As_TRD,As_TRD[As_TRD$cat == '0.46 - 0.50', ]$As, col='blue', pch=16)
points(As_TRD[As_TRD$cat == '0.50 - 0.54', ]$As_TRD,As_TRD[As_TRD$cat == '0.50 - 0.54', ]$As, col='orange', pch=16)
points(As_TRD[As_TRD$cat == '0.54 - 0.59', ]$As_TRD,As_TRD[As_TRD$cat == '0.54 - 0.59', ]$As, col='brown', pch=16)
points(As_TRD[As_TRD$cat == '0.59 - 0.76', ]$As_TRD,As_TRD[As_TRD$cat == '0.59 - 0.76', ]$As, col='red', pch=16)

lines(sort(As_TRD[As_TRD$cat == '0.21 - 0.46', ]$As_TRD), exp(FixedEff[1]+ (RIaS[7]  + FixedEff[2]) * sort(As_TRD[As_TRD$cat == '0.21 - 0.46', ]$As_TRD) + RIaS[1]), col='black' , lwd=4)
lines(sort(As_TRD[As_TRD$SITE == 'DL', ]$As_TRD), exp(FixedEff[1]+ (RIaS[8]  + FixedEff[2]) * sort(As_TRD[As_TRD$SITE == 'DL', ]$As_TRD) + RIaS[2]), col='red'   , lwd=4)
lines(sort(As_TRD[As_TRD$SITE == 'GR', ]$As_TRD), exp(FixedEff[1]+ (RIaS[9]  + FixedEff[2]) * sort(As_TRD[As_TRD$SITE == 'GR', ]$As_TRD) + RIaS[3]), col='green' , lwd=4)
lines(sort(As_TRD[As_TRD$SITE == 'GC', ]$As_TRD), exp(FixedEff[1]+ (RIaS[10] + FixedEff[2]) * sort(As_TRD[As_TRD$SITE == 'GC', ]$As_TRD) + RIaS[4]), col='blue'  , lwd=4)
lines(sort(As_TRD[As_TRD$SITE == 'BG', ]$As_TRD), exp(FixedEff[1]+ (RIaS[11] + FixedEff[2]) * sort(As_TRD[As_TRD$SITE == 'BG', ]$As_TRD) + RIaS[5]), col='orange', lwd=4) 
lines(sort(As_TRD[As_TRD$SITE == 'BN', ]$As_TRD), exp(FixedEff[1]+ (RIaS[12] + FixedEff[2]) * sort(As_TRD[As_TRD$SITE == 'BN', ]$As_TRD) + RIaS[6]), col='brown' , lwd=4) 

#abline(v=(seq(-1,4 ,1)), col="lightgray", lty="dotted");        
#abline(h=(seq( -1,25 ,1)), col="lightgray", lty="dotted")   

