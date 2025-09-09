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
library(glmnet)


## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_SEP_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

COMPARTMENTS_AVG[,"TURNOVER"]<-COMPARTMENTS_AVG[,"NEW"]/(2*COMPARTMENTS_AVG[,"OM.AREA.g.m2"])

replace_outliers_with_na <- function(x, threshold = 3) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}


#### As ####

As_DF <- (COMPARTMENTS_AVG[,c("spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","OM.AREA.g.m2","NEW","As")])

As_DF <- As_DF[complete.cases(As_DF),]

As_DF<-As_DF[As_DF$OM.AREA.g.m2>0,]

As_DF$OM.AREA.g.m2 <- log10(As_DF$OM.AREA.g.m2+1)

names(As_DF)<-c("As_SPM","As_COL","As_TRD","OM","TURNOVER", "As")

sum(As_DF$As_SPM, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_COL, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_TRD, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))

As_DF <- As_DF[,-c(1,2)]

#As_DF$TURNOVER <- sqrt(As_DF$TURNOVER)

As.full.model <- glm(As ~. ^2, data = As_DF, family = gaussian)
coef(As.full.model)

As.step.model <- As.full.model %>% stepAIC(trace = T)
coef(As.step.model)

Anova(As.step.model, test="LR")
summary(As.step.model)
with(summary(As.step.model), 1 - deviance/null.deviance)

plot(allEffects(As.step.model))

calc.relimp(As.step.model, rela=T)

y <- As_DF$As

x <- data.matrix(As_DF[, c(1:ncol(As_DF)-1)])

x <- scale(x,center=T, scale=T)

y <- scale(y,center=T, scale=T)

cv_model <- cv.glmnet(x, y, alpha = 1)

best_lambda <- cv_model$lambda.min

best_lambda

plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)

coef(best_model)

y_predicted <- predict(best_model, s = best_lambda, newx = x)

plot(As_DF$As_TRD, y)

plot(As_DF$OM , y)

plot(As_DF$TURNOVER, y)

#find SST and SSE

sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared

rsq <- 1 - sse/sst
rsq

#### Cd ####

Cd_DF <- (COMPARTMENTS_AVG[,c("spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","OM.AREA.g.m2","NEW","Cd")])

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

Cd_DF<-Cd_DF[Cd_DF$OM.AREA.g.m2>0,]
Cd_DF$OM.AREA.g.m2 <- log10(Cd_DF$OM.AREA.g.m2+1)

names(Cd_DF)<-c("Cd_SPM","Cd_COL","Cd_TRD","OM","TURNOVER", "Cd")

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

#Cd_DF$TURNOVER <- sqrt(Cd_DF$TURNOVER)

cd.full.model <- glm(Cd ~. ^2, data = Cd_DF, family = gaussian)
coef(cd.full.model)

cd.step.model <- cd.full.model %>% stepAIC(trace = T)
coef(cd.step.model)

Anova(cd.step.model, test="LR")
summary(cd.step.model)
with(summary(cd.step.model), 1 - deviance/null.deviance)

plot(allEffects(cd.step.model))

calc.relimp(cd.step.model)

y <- Cd_DF$Cd

x <- data.matrix(Cd_DF[, c(1:ncol(Cd_DF)-1)])

x <- scale(x,center=T, scale=T)

y <- scale(y,center=T, scale=T)

cv_model <- cv.glmnet(x, y, alpha = 1)

best_lambda <- cv_model$lambda.min

best_lambda

plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)

coef(best_model)

y_predicted <- predict(best_model, s = best_lambda, newx = x)

#find SST and SSE

sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared

rsq <- 1 - sse/sst
rsq






