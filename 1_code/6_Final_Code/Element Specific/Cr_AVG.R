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


#### Cr ####

Cr_DF <- (COMPARTMENTS_AVG[,c("W_Cr_NH3","F_Cr_NH3","trulydissolved_Cr_NH3","OM.AREA.g.m2","TURNOVER","Cr")])

Cr_DF<-Cr_DF[Cr_DF$OM.AREA.g.m2>0,]

Cr_DF$OM.AREA.g.m2 <- log10(Cr_DF$OM.AREA.g.m2+1)

Cr_DF <- Cr_DF[complete.cases(Cr_DF),]

#Cr_DF$TURNOVER <- sqrt(Cr_DF$NEW)

quantile <- quantile(Cr_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cr_DF$cat <- cut(Cr_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### Cr_SPM ###

Cr_SPM <- (Cr_DF[,c(1,4,5,6,7)])

names(Cr_SPM)<- c("Cr_SPM", "TURNOVER", "OM" ,"Cr", "cat")

Cr_SPM<-Cr_SPM[Cr_SPM$Cr_SPM>0,]

Cr_SPM$Cr_SPM <- log10(Cr_SPM$Cr_SPM+1)

for (col in colnames(Cr_SPM[,c(1)])) {
  Cr_SPM[[col]] <- replace_outliers_with_na(Cr_SPM[[col]])
}

Cr_SPM <- Cr_SPM[complete.cases(Cr_SPM), ]

### Cr_COL ###

Cr_COL <- (Cr_DF[,c(2,4,5,6,7)])

names(Cr_COL)<- c("Cr_COL", "TURNOVER", "OM","Cr", "cat")

Cr_COL<-Cr_COL[Cr_COL$Cr_COL>0,]

Cr_COL$Cr_COL <- log10(Cr_COL$Cr_COL+1)

for (col in colnames(Cr_COL[,c(1)])) {
  Cr_COL[[col]] <- replace_outliers_with_na(Cr_COL[[col]])
}

Cr_COL <- Cr_COL[complete.cases(Cr_COL), ]

### Cr_TRD ###

Cr_TRD <- (Cr_DF[,c(3,4,5,6,7)])

names(Cr_TRD)<- c("Cr_TRD", "TURNOVER", "OM" ,"Cr", "cat")

Cr_TRD<-Cr_TRD[Cr_TRD$Cr_TRD>0,]

Cr_TRD$Cr_TRD <- log10(Cr_TRD$Cr_TRD+1)

for (col in colnames(Cr_TRD[,c(1)])) {
  Cr_TRD[[col]] <- replace_outliers_with_na(Cr_TRD[[col]])
}

Cr_TRD <- Cr_TRD[complete.cases(Cr_TRD), ]

#dev.new()
plot(density(Cr_SPM$Cr_SPM, bw=0.09), ylim=c(0,5), xlim=c(0,1), col="Purple", lwd=2)
lines(density(Cr_COL$Cr_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Cr_TRD$Cr_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cr_SPM <- glm(Cr ~ (Cr_SPM + TURNOVER + OM), data = Cr_SPM, family=gaussian)

mix.int_Cr_COL <- glm(Cr ~ (Cr_COL + TURNOVER + OM), data = Cr_COL, family=gaussian)

mix.int_Cr_TRD <- glm(Cr ~ (TURNOVER + OM), data = Cr_TRD, family=gaussian)

plot(allEffects(mix.int_Cr_SPM))

plot(allEffects(mix.int_Cr_COL))

plot(allEffects(mix.int_Cr_TRD))

Anova(mix.int_Cr_SPM, test="LR")
summary(mix.int_Cr_SPM)
with(summary(mix.int_Cr_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Cr_COL, test="LR")
summary(mix.int_Cr_COL)
with(summary(mix.int_Cr_COL), 1 - deviance/null.deviance)

Anova(mix.int_Cr_TRD, test="LR")
summary(mix.int_Cr_TRD)
with(summary(mix.int_Cr_TRD), 1 - deviance/null.deviance)

### Cr Crl Plots ###

Cr_p_SPM=ggplot(data=Cr_SPM,aes(Cr_SPM,Cr))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cr_SPM, y=Cr, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cr  (SPM, mg/l)",    
       y= "Cr Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cr_p_COL=ggplot(data=Cr_COL,aes(Cr_COL,Cr))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cr_COL, y=Cr, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cr  (COL, mg/l)",    
       y= "Cr Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cr_p_TRD=ggplot(data=Cr_TRD,aes(Cr_TRD,Cr))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cr_TRD, y=Cr, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cr  (TD, mg/l)",    
       y= "Cr Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Cr_p_SPM,Cr_p_COL,Cr_p_TRD,nrow=1)

#### % Total ####

names(Cr_DF)<-c("Cr_SPM","Cr_COL","Cr_TRD","OM","TURNOVER", "Al", "cat")
Cr_DF[,c(1:3)][Cr_DF[,c(1:3)] <= 0] <- 0

sum(Cr_DF$Cr_SPM, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))
sum(Cr_DF$Cr_COL, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))
sum(Cr_DF$Cr_TRD, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))

