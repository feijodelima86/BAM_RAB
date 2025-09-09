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

As_DF <- (COMPARTMENTS_AVG[,c("W_As_Oshift","F_As_Oshift","trulydissolved_As_Oshift","OM.AREA.g.m2","TURNOVER","As")])

As_DF<-As_DF[As_DF$OM.AREA.g.m2>0,]

As_DF$OM.AREA.g.m2 <- log10(As_DF$OM.AREA.g.m2+1)

As_DF <- As_DF[complete.cases(As_DF),]

#As_DF$TURNOVER <- sqrt(As_DF$NEW)

quantile <- quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### As_SPM ###

As_SPM <- (As_DF[,c(1,4,5,6,7)])

names(As_SPM)<- c("As_SPM", "TURNOVER", "OM" ,"As", "cat")

As_SPM<-As_SPM[As_SPM$As_SPM>0,]

As_SPM$As_SPM <- log10(As_SPM$As_SPM+1)

for (col in colnames(As_SPM[,c(1)])) {
  As_SPM[[col]] <- replace_outliers_with_na(As_SPM[[col]])
}

As_SPM <- As_SPM[complete.cases(As_SPM), ]

### As_COL ###

As_COL <- (As_DF[,c(2,4,5,6,7)])

names(As_COL)<- c("As_COL", "TURNOVER", "OM","As", "cat")

As_COL<-As_COL[As_COL$As_COL>0,]

As_COL$As_COL <- log10(As_COL$As_COL+1)

for (col in colnames(As_COL[,c(1)])) {
  As_COL[[col]] <- replace_outliers_with_na(As_COL[[col]])
}

As_COL <- As_COL[complete.cases(As_COL), ]

### As_TRD ###

As_TRD <- (As_DF[,c(3,4,5,6,7)])

names(As_TRD)<- c("As_TRD", "TURNOVER", "OM" ,"As", "cat")

As_TRD<-As_TRD[As_TRD$As_TRD>0,]

As_TRD$As_TRD <- log10(As_TRD$As_TRD+1)

for (col in colnames(As_TRD[,c(1)])) {
  As_TRD[[col]] <- replace_outliers_with_na(As_TRD[[col]])
}

As_TRD <- As_TRD[complete.cases(As_TRD), ]

#dev.new()
plot(density(As_SPM$As_SPM, bw=0.09), ylim=c(0,3), xlim=c(0,3), col="Purple", lwd=2)
lines(density(As_COL$As_COL, bw=0.09), col="Blue", lwd=2)
lines(density(As_TRD$As_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_As_SPM <- glm(As ~ (As_SPM + TURNOVER + OM), data = As_SPM, family=gaussian)

mix.int_As_COL <- glm(As ~ (As_COL + TURNOVER + OM), data = As_COL, family=gaussian)

mix.int_As_TRD <- glm(As ~ (As_TRD + TURNOVER + OM), data = As_TRD, family=gaussian)

plot(allEffects(mix.int_As_SPM))

plot(allEffects(mix.int_As_COL))

plot(allEffects(mix.int_As_TRD))

Anova(mix.int_As_SPM, test="LR")
summary(mix.int_As_SPM)
with(summary(mix.int_As_SPM), 1 - deviance/null.deviance)

Anova(mix.int_As_COL, test="LR")
summary(mix.int_As_COL)
with(summary(mix.int_As_COL), 1 - deviance/null.deviance)

Anova(mix.int_As_TRD, test="LR")
summary(mix.int_As_TRD)
with(summary(mix.int_As_TRD), 1 - deviance/null.deviance)

### As Asl Plots ###

As_p_SPM=ggplot(data=As_SPM,aes(As_SPM,As))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=As_SPM, y=As, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "As  (SPM, mg/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

As_p_COL=ggplot(data=As_COL,aes(As_COL,As))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=As_COL, y=As, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "As  (COL, mg/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

As_p_TRD=ggplot(data=As_TRD,aes(As_TRD,As))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=As_TRD, y=As, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "As  (TD, mg/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(As_p_SPM,As_p_COL,As_p_TRD,nrow=1)

#### % Total ####

names(As_DF)<-c("As_SPM","As_COL","As_TRD","OM","TURNOVER", "Al", "cat")
As_DF[,c(1:3)][As_DF[,c(1:3)] <= 0] <- 0

sum(As_DF$As_SPM, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_COL, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_TRD, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))

