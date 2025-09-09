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


#### Cd ####

Cd_DF <- (COMPARTMENTS_AVG[,c("W_Cd114","F_Cd114","trulydissolved_Cd114","OM.AREA.g.m2","TURNOVER","Cd")])

Cd_DF<-Cd_DF[Cd_DF$OM.AREA.g.m2>0,]

Cd_DF$OM.AREA.g.m2 <- log10(Cd_DF$OM.AREA.g.m2+1)

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

#Cd_DF$TURNOVER <- sqrt(Cd_DF$NEW)

quantile <- quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### Cd_SPM ###

Cd_SPM <- (Cd_DF[,c(1,4,5,6,7)])

names(Cd_SPM)<- c("Cd_SPM", "TURNOVER", "OM" ,"Cd", "cat")

Cd_SPM<-Cd_SPM[Cd_SPM$Cd_SPM>0,]

Cd_SPM$Cd_SPM <- log10(Cd_SPM$Cd_SPM+1)

for (col in colnames(Cd_SPM[,c(1)])) {
  Cd_SPM[[col]] <- replace_outliers_with_na(Cd_SPM[[col]])
}

Cd_SPM <- Cd_SPM[complete.cases(Cd_SPM), ]

### Cd_COL ###

Cd_COL <- (Cd_DF[,c(2,4,5,6,7)])

names(Cd_COL)<- c("Cd_COL", "TURNOVER", "OM","Cd", "cat")

Cd_COL<-Cd_COL[Cd_COL$Cd_COL>0,]

Cd_COL$Cd_COL <- log10(Cd_COL$Cd_COL+1)

for (col in colnames(Cd_COL[,c(1)])) {
  Cd_COL[[col]] <- replace_outliers_with_na(Cd_COL[[col]])
}

Cd_COL <- Cd_COL[complete.cases(Cd_COL), ]

### Cd_TRD ###

Cd_TRD <- (Cd_DF[,c(3,4,5,6,7)])

names(Cd_TRD)<- c("Cd_TRD", "TURNOVER", "OM" ,"Cd", "cat")

Cd_TRD<-Cd_TRD[Cd_TRD$Cd_TRD>0,]

Cd_TRD$Cd_TRD <- log10(Cd_TRD$Cd_TRD+1)

for (col in colnames(Cd_TRD[,c(1)])) {
  Cd_TRD[[col]] <- replace_outliers_with_na(Cd_TRD[[col]])
}

Cd_TRD <- Cd_TRD[complete.cases(Cd_TRD), ]

#dev.new()
plot(density(Cd_SPM$Cd_SPM, bw=0.09), ylim=c(0,1), xlim=c(0,1), col="Purple", lwd=2)
lines(density(Cd_COL$Cd_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Cd_TRD$Cd_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cd_SPM <- glm(Cd ~ (Cd_SPM + TURNOVER + OM), data = Cd_SPM, family=gaussian)

mix.int_Cd_COL <- glm(Cd ~ (Cd_COL + TURNOVER + OM), data = Cd_COL, family=gaussian)

mix.int_Cd_TRD <- glm(Cd ~ (Cd_TRD + TURNOVER + OM), data = Cd_TRD, family=gaussian)




plot(allEffects(mix.int_Cd_SPM))

plot(allEffects(mix.int_Cd_COL))

plot(allEffects(mix.int_Cd_TRD))

Anova(mix.int_Cd_SPM, test="LR")
summary(mix.int_Cd_SPM)
with(summary(mix.int_Cd_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Cd_COL, test="LR")
summary(mix.int_Cd_COL)
with(summary(mix.int_Cd_COL), 1 - deviance/null.deviance)

Anova(mix.int_Cd_TRD, test="LR")
summary(mix.int_Cd_TRD)
with(summary(mix.int_Cd_TRD), 1 - deviance/null.deviance)

### Cd Cdl Plots ###

Cd_p_SPM=ggplot(data=Cd_SPM,aes(Cd_SPM,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_SPM, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (SPM, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cd_p_COL=ggplot(data=Cd_COL,aes(Cd_COL,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_COL, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (COL, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cd_p_TRD=ggplot(data=Cd_TRD,aes(Cd_TRD,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_TRD, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (TD, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Cd_p_SPM,Cd_p_COL,Cd_p_TRD,nrow=1)

#### % Total ####

names(Cd_DF)<-c("Cd_SPM","Cd_COL","Cd_TRD","OM","TURNOVER", "Al", "cat")
Cd_DF[,c(1:3)][Cd_DF[,c(1:3)] <= 0] <- 0

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

