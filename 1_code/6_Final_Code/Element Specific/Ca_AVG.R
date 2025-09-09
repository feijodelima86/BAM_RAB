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


#### Ca ####

Ca_DF <- (COMPARTMENTS_AVG[,c("spm_Ca","colloidal_Ca","trulydissolved_Ca","OM.AREA.g.m2","TURNOVER","Car")])

Ca_DF<-Ca_DF[Ca_DF$OM.AREA.g.m2>0,]

Ca_DF$OM.AREA.g.m2 <- log10(Ca_DF$OM.AREA.g.m2+1)

Ca_DF <- Ca_DF[complete.cases(Ca_DF),]

#Ca_DF$TURNOVER <- sqrt(Ca_DF$NEW)

quantile <- quantile(Ca_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ca_DF$cat <- cut(Ca_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### Ca_SPM ###

Ca_SPM <- (Ca_DF[,c(1,4,5,6,7)])

names(Ca_SPM)<- c("Ca_SPM", "TURNOVER", "OM" ,"Ca", "cat")

Ca_SPM<-Ca_SPM[Ca_SPM$Ca_SPM>0,]

Ca_SPM$Ca_SPM <- log10(Ca_SPM$Ca_SPM+1)

for (col in colnames(Ca_SPM[,c(1)])) {
  Ca_SPM[[col]] <- replace_outliers_with_na(Ca_SPM[[col]])
}

Ca_SPM <- Ca_SPM[complete.cases(Ca_SPM), ]

### Ca_COL ###

Ca_COL <- (Ca_DF[,c(2,4,5,6,7)])

names(Ca_COL)<- c("Ca_COL", "TURNOVER", "OM","Ca", "cat")

Ca_COL<-Ca_COL[Ca_COL$Ca_COL>0,]

Ca_COL$Ca_COL <- log10(Ca_COL$Ca_COL+1)

for (col in colnames(Ca_COL[,c(1)])) {
  Ca_COL[[col]] <- replace_outliers_with_na(Ca_COL[[col]])
}

Ca_COL <- Ca_COL[complete.cases(Ca_COL), ]

### Ca_TRD ###

Ca_TRD <- (Ca_DF[,c(3,4,5,6,7)])

names(Ca_TRD)<- c("Ca_TRD", "TURNOVER", "OM" ,"Ca", "cat")

Ca_TRD<-Ca_TRD[Ca_TRD$Ca_TRD>0,]

Ca_TRD$Ca_TRD <- log10(Ca_TRD$Ca_TRD+1)

for (col in colnames(Ca_TRD[,c(1)])) {
  Ca_TRD[[col]] <- replace_outliers_with_na(Ca_TRD[[col]])
}

Ca_TRD <- Ca_TRD[complete.cases(Ca_TRD), ]

#dev.new()
plot(density(Ca_SPM$Ca_SPM, bw=0.09), ylim=c(0,5), xlim=c(0,1), col="Purple", lwd=2)
lines(density(Ca_COL$Ca_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Ca_TRD$Ca_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Ca_SPM <- glm(Ca ~ (Ca_SPM + TURNOVER + OM), data = Ca_SPM, family=gaussian)

mix.int_Ca_COL <- glm(Ca ~ (Ca_COL + TURNOVER + OM), data = Ca_COL, family=gaussian)

mix.int_Ca_TRD <- glm(Ca ~ (Ca_TRD + TURNOVER + OM), data = Ca_TRD, family=gaussian)

plot(allEffects(mix.int_Ca_SPM))

plot(allEffects(mix.int_Ca_COL))

plot(allEffects(mix.int_Ca_TRD))

Anova(mix.int_Ca_SPM, test="LR")
summary(mix.int_Ca_SPM)
with(summary(mix.int_Ca_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Ca_COL, test="LR")
summary(mix.int_Ca_COL)
with(summary(mix.int_Ca_COL), 1 - deviance/null.deviance)

Anova(mix.int_Ca_TRD, test="LR")
summary(mix.int_Ca_TRD)
with(summary(mix.int_Ca_TRD), 1 - deviance/null.deviance)

### Ca Cal Plots ###

Ca_p_SPM=ggplot(data=Ca_SPM,aes(Ca_SPM,Ca))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ca_SPM, y=Ca, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ca  (SPM, mg/l)",    
       y= "Ca Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ca_p_COL=ggplot(data=Ca_COL,aes(Ca_COL,Ca))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ca_COL, y=Ca, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ca  (COL, mg/l)",    
       y= "Ca Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ca_p_TRD=ggplot(data=Ca_TRD,aes(Ca_TRD,Ca))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ca_TRD, y=Ca, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ca  (TD, mg/l)",    
       y= "Ca Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Ca_p_SPM,Ca_p_COL,Ca_p_TRD,nrow=1)

#### % Total ####

names(Ca_DF)<-c("Ca_SPM","Ca_COL","Ca_TRD","OM","TURNOVER", "Al", "cat")
Ca_DF[,c(1:3)][Ca_DF[,c(1:3)] <= 0] <- 0

sum(Ca_DF$Ca_SPM, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))
sum(Ca_DF$Ca_COL, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))
sum(Ca_DF$Ca_TRD, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))

