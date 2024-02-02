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

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

replace_outliers_with_na <- function(x, threshold = 3) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}



#### Al ####

Al_DF <- (COMPARTMENTS_AVG[,c("spm_Al","colloidal_Al","trulydissolved_Al","TURNOVER","Al")])

Al_DF <- Al_DF[complete.cases(Al_DF),]

Al_DF$TURNOVER <- asin(sqrt(Al_DF$TURNOVER ))

quantile(Al_DF$TURNOVER, probs = seq(0, 1, 1/5))

Al_DF$cat <- cut(Al_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Al_SPM ###

Al_SPM <- (Al_DF[,c(1,4,5,6)])

names(Al_SPM)<- c("Al_SPM", "TURNOVER", "Al", "cat")

Al_SPM<-Al_SPM[Al_SPM$Al_SPM>0,]

Al_SPM$Al_SPM <- log10(Al_SPM$Al_SPM+1)

for (col in colnames(Al_SPM[,c(1)])) {
  Al_SPM[[col]] <- replace_outliers_with_na(Al_SPM[[col]])
}

Al_SPM <- Al_SPM[complete.cases(Al_SPM), ]

### Al_COL ###

Al_COL <- (Al_DF[,c(2,4,5,6)])

names(Al_COL)<- c("Al_COL", "TURNOVER", "Al", "cat")

Al_COL<-Al_COL[Al_COL$Al_COL>0,]

Al_COL$Al_COL <- log10(Al_COL$Al_COL+1)

for (col in colnames(Al_COL[,c(1)])) {
  Al_COL[[col]] <- replace_outliers_with_na(Al_COL[[col]])
}

Al_COL <- Al_COL[complete.cases(Al_COL), ]

### Al_TRD ###

Al_TRD <- (Al_DF[,c(3,4,5,6)])

names(Al_TRD)<- c("Al_TRD", "TURNOVER", "Al", "cat")

Al_TRD<-Al_TRD[Al_TRD$Al_TRD>0,]

Al_TRD$Al_TRD <- log10(Al_TRD$Al_TRD+1)

for (col in colnames(Al_TRD[,c(1)])) {
  Al_TRD[[col]] <- replace_outliers_with_na(Al_TRD[[col]])
}

Al_TRD <- Al_TRD[complete.cases(Al_TRD), ]

#dev.new()
plot(density(Al_SPM$Al_SPM))
plot(density(Al_SPM$Al_SPM, bw=0.09), ylim=c(0,3), xlim=c(0,3), col="Purple", lwd=2)
lines(density(Al_COL$Al_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Al_TRD$Al_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Al_SPM <- glm(Al ~ (Al_SPM * TURNOVER), data = Al_SPM, family=gaussian)

mix.int_Al_COL <- glm(Al ~ (Al_COL * TURNOVER), data = Al_COL, family=gaussian)

mix.int_Al_TRD <- glm(Al ~ (Al_TRD * TURNOVER), data = Al_TRD, family=gaussian)

plot(allEffects(mix.int_Al_SPM))

plot(allEffects(mix.int_Al_COL))

plot(allEffects(mix.int_Al_TRD))

Anova(mix.int_Al_SPM, test="LR")
summary(mix.int_Al_SPM)
with(summary(mix.int_Al_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Al_COL, test="LR")
summary(mix.int_Al_COL)
with(summary(mix.int_Al_COL), 1 - deviance/null.deviance)

Anova(mix.int_Al_TRD, test="LR")
summary(mix.int_Al_TRD)
with(summary(mix.int_Al_TRD), 1 - deviance/null.deviance)

### Al All Plots ###

Al_p_SPM=ggplot(data=Al_SPM,aes(Al_SPM,Al))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Al_SPM, y=Al, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Al  (SPM, mg/l)",    
       y= "Al Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Al_p_COL=ggplot(data=Al_COL,aes(Al_COL,Al))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Al_COL, y=Al, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Al  (COL, mg/l)",    
       y= "Al Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Al_p_TRD=ggplot(data=Al_TRD,aes(Al_TRD,Al))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Al_TRD, y=Al, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Al  (TD, mg/l)",    
       y= "Al Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Al_p_SPM,Al_p_COL,Al_p_TRD,nrow=1)


#dev.new()
#### As ####

As_DF <- (COMPARTMENTS_AVG[,c("spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","TURNOVER","As")])

As_DF

As_DF <- As_DF[complete.cases(As_DF),]

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### As_SPM ###

As_SPM <- (As_DF[,c(1,4,5,6)])

names(As_SPM)<- c("As_SPM", "TURNOVER", "As", "cat")

As_SPM<-As_SPM[As_SPM$As_SPM>0,]

As_SPM$As_SPM <- log10(As_SPM$As_SPM+1)

for (col in colnames(As_SPM[,c(1)])) {
  As_SPM[[col]] <- replace_outliers_with_na(As_SPM[[col]])
}

As_SPM <- As_SPM[complete.cases(As_SPM), ]

### As_COL ###

As_COL <- (As_DF[,c(2,4,5,6)])

names(As_COL)<- c("As_COL", "TURNOVER", "As", "cat")

As_COL<-As_COL[As_COL$As_COL>0,]

As_COL$As_COL <- log10(As_COL$As_COL+1)

for (col in colnames(As_COL[,c(1)])) {
  As_COL[[col]] <- replace_outliers_with_na(As_COL[[col]])
}

As_COL <- As_COL[complete.cases(As_COL), ]

### As_TRD ###

As_TRD <- (As_DF[,c(3,4,5,6)])

names(As_TRD)<- c("As_TRD", "TURNOVER", "As", "cat")

As_TRD<-As_TRD[As_TRD$As_TRD>0,]

As_TRD$As_TRD <- log10(As_TRD$As_TRD+1)

for (col in colnames(As_TRD[,c(1)])) {
  As_TRD[[col]] <- replace_outliers_with_na(As_TRD[[col]])
}

As_TRD <- As_TRD[complete.cases(As_TRD), ]

#dev.new()

plot(density(As_SPM$As_SPM, bw=0.051), ylim=c(0,4), xlim=c(0,1.4), col="Purple", lwd=2)
lines(density(As_COL$As_COL, bw=0.051), col="Blue", lwd=2)
lines(density(As_TRD$As_TRD, bw=0.051), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_As_SPM <- glm(As ~ (As_SPM * TURNOVER), data = As_SPM, family=gaussian)

mix.int_As_COL <- glm(As ~ (As_COL * TURNOVER), data = As_COL, family=gaussian)

mix.int_As_TRD <- glm(As ~ (As_TRD * TURNOVER), data = As_TRD, family=gaussian)

plot(allEffects(mix.int_As_SPM))

plot(allEffects(mix.int_As_COL))

plot(allEffects(mix.int_As_TRD))

Anova(mix.int_As_SPM)
summary(mix.int_As_SPM)
with(summary(mix.int_As_SPM), 1 - deviance/null.deviance)

Anova(mix.int_As_COL)
summary(mix.int_As_COL)
with(summary(mix.int_As_COL), 1 - deviance/null.deviance)

Anova(mix.int_As_TRD)
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


#dev.new()


#### Ca ####

Ca_DF <- (COMPARTMENTS_AVG[,c("spm_Ca","colloidal_Ca","trulydissolved_Ca","TURNOVER","Car")])

Ca_DF <- Ca_DF[complete.cases(Ca_DF),]

Ca_DF$TURNOVER <- asin(sqrt(Ca_DF$TURNOVER ))

quantile(Ca_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ca_DF$cat <- cut(Ca_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Ca_SPM ###

Ca_SPM <- (Ca_DF[,c(1,4,5,6)])

names(Ca_SPM)<- c("Ca_SPM", "TURNOVER", "Ca", "cat")

Ca_SPM<-Ca_SPM[Ca_SPM$Ca_SPM>0,]

Ca_SPM$Ca_SPM <- log10(Ca_SPM$Ca_SPM+1)

for (col in colnames(Ca_SPM[,c(1)])) {
  Ca_SPM[[col]] <- replace_outliers_with_na(Ca_SPM[[col]])
}

Ca_SPM <- Ca_SPM[complete.cases(Ca_SPM), ]

### Ca_COL ###

Ca_COL <- (Ca_DF[,c(2,4,5,6)])

names(Ca_COL)<- c("Ca_COL", "TURNOVER", "Ca", "cat")

Ca_COL<-Ca_COL[Ca_COL$Ca_COL>0,]

Ca_COL$Ca_COL <- log10(Ca_COL$Ca_COL+1)

for (col in colnames(Ca_COL[,c(1)])) {
  Ca_COL[[col]] <- replace_outliers_with_na(Ca_COL[[col]])
}

Ca_COL <- Ca_COL[complete.cases(Ca_COL), ]

### Ca_TRD ###

Ca_TRD <- (Ca_DF[,c(3,4,5,6)])

names(Ca_TRD)<- c("Ca_TRD", "TURNOVER", "Ca", "cat")

Ca_TRD<-Ca_TRD[Ca_TRD$Ca_TRD>0,]

Ca_TRD$Ca_TRD <- log10(Ca_TRD$Ca_TRD+1)

for (col in colnames(Ca_TRD[,c(1)])) {
  Ca_TRD[[col]] <- replace_outliers_with_na(Ca_TRD[[col]])
}

Ca_TRD <- Ca_TRD[complete.cases(Ca_TRD), ]

#dev.new()

plot(density(Ca_SPM$Ca_SPM, bw=0.1726), ylim=c(0,2.5), xlim=c(1.7,5.5), col="Purple", lwd=2)
lines(density(Ca_COL$Ca_COL, bw=0.1726), col="Blue", lwd=2)
lines(density(Ca_TRD$Ca_TRD, bw=0.1726), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Ca_SPM <- glm(Ca ~ (Ca_SPM * TURNOVER), data = Ca_SPM, family=gaussian)

mix.int_Ca_COL <- glm(Ca ~ (Ca_COL * TURNOVER), data = Ca_COL, family=gaussian)

mix.int_Ca_TRD <- glm(Ca ~ (Ca_TRD * TURNOVER), data = Ca_TRD, family=gaussian)

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


#dev.new()

#### Cd ####

Cd_DF <- (COMPARTMENTS_AVG[,c("spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","TURNOVER","Cd")])

Cd_DF <- Cd_DF[complete.cases(Cd_DF),]

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Cd_SPM ###

Cd_SPM <- (Cd_DF[,c(1,4,5,6)])

names(Cd_SPM)<- c("Cd_SPM", "TURNOVER", "Cd", "cat")

Cd_SPM<-Cd_SPM[Cd_SPM$Cd_SPM>0,]

Cd_SPM$Cd_SPM <- log10(Cd_SPM$Cd_SPM+1)

for (col in colnames(Cd_SPM[,c(1)])) {
  Cd_SPM[[col]] <- replace_outliers_with_na(Cd_SPM[[col]])
}

Cd_SPM <- Cd_SPM[complete.cases(Cd_SPM), ]

### Cd_COL ###

Cd_COL <- (Cd_DF[,c(2,4,5,6)])

names(Cd_COL)<- c("Cd_COL", "TURNOVER", "Cd", "cat")

Cd_COL<-Cd_COL[Cd_COL$Cd_COL>0,]

Cd_COL$Cd_COL <- log10(Cd_COL$Cd_COL+1)

for (col in colnames(Cd_COL[,c(1)])) {
  Cd_COL[[col]] <- replace_outliers_with_na(Cd_COL[[col]])
}

Cd_COL <- Cd_COL[complete.cases(Cd_COL), ]

### Cd_TRD ###

Cd_TRD <- (Cd_DF[,c(3,4,5,6)])

names(Cd_TRD)<- c("Cd_TRD", "TURNOVER", "Cd", "cat")

Cd_TRD<-Cd_TRD[Cd_TRD$Cd_TRD>0,]

Cd_TRD$Cd_TRD <- log10(Cd_TRD$Cd_TRD+1)

for (col in colnames(Cd_TRD[,c(1)])) {
  Cd_TRD[[col]] <- replace_outliers_with_na(Cd_TRD[[col]])
}

Cd_TRD <- Cd_TRD[complete.cases(Cd_TRD), ]

#dev.new()

plot(density(Cd_SPM$Cd_SPM, bw=0.0029), ylim=c(0,150), xlim=c(0,0.04), col="Purple", lwd=2)
lines(density(Cd_COL$Cd_COL, bw=0.0029), col="Blue", lwd=2)
lines(density(Cd_TRD$Cd_TRD, bw=0.0029), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cd_SPM <- glm(Cd ~ (Cd_SPM * TURNOVER), data = Cd_SPM, family=gaussian)

mix.int_Cd_COL <- glm(Cd ~ (Cd_COL * TURNOVER), data = Cd_COL, family=gaussian)

mix.int_Cd_TRD <- glm(Cd ~ (Cd_TRD * TURNOVER), data = Cd_TRD, family=gaussian)

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


#dev.new()


#### Co ####

Co_DF <- (COMPARTMENTS_AVG[,c("spm_Co","colloidal_Co","trulydissolved_Co","TURNOVER","Co")])

Co_DF <- Co_DF[complete.cases(Co_DF),]

Co_DF$TURNOVER <- asin(sqrt(Co_DF$TURNOVER ))

quantile(Co_DF$TURNOVER, probs = seq(0, 1, 1/5))

Co_DF$cat <- cut(Co_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Co_SPM ###

Co_SPM <- (Co_DF[,c(1,4,5,6)])

names(Co_SPM)<- c("Co_SPM", "TURNOVER", "Co", "cat")

Co_SPM<-Co_SPM[Co_SPM$Co_SPM>0,]

Co_SPM$Co_SPM <- log10(Co_SPM$Co_SPM+1)

for (col in colnames(Co_SPM[,c(1)])) {
  Co_SPM[[col]] <- replace_outliers_with_na(Co_SPM[[col]])
}

Co_SPM <- Co_SPM[complete.cases(Co_SPM), ]

### Co_COL ###

Co_COL <- (Co_DF[,c(2,4,5,6)])

names(Co_COL)<- c("Co_COL", "TURNOVER", "Co", "cat")

Co_COL<-Co_COL[Co_COL$Co_COL>0,]

Co_COL$Co_COL <- log10(Co_COL$Co_COL+1)

for (col in colnames(Co_COL[,c(1)])) {
  Co_COL[[col]] <- replace_outliers_with_na(Co_COL[[col]])
}

Co_COL <- Co_COL[complete.cases(Co_COL), ]

### Co_TRD ###

Co_TRD <- (Co_DF[,c(3,4,5,6)])

names(Co_TRD)<- c("Co_TRD", "TURNOVER", "Co", "cat")

Co_TRD<-Co_TRD[Co_TRD$Co_TRD>0,]

Co_TRD$Co_TRD <- log10(Co_TRD$Co_TRD+1)

for (col in colnames(Co_TRD[,c(1)])) {
  Co_TRD[[col]] <- replace_outliers_with_na(Co_TRD[[col]])
}

Co_TRD <- Co_TRD[complete.cases(Co_TRD), ]

#dev.new()

plot(density(Co_SPM$Co_SPM, bw=0.0029), ylim=c(0,100), xlim=c(0,0.04), col="Purple", lwd=2)
lines(density(Co_COL$Co_COL, bw=0.0029), col="Blue", lwd=2)
lines(density(Co_TRD$Co_TRD, bw=0.0029), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Co_SPM <- glm(Co ~ (Co_SPM * TURNOVER), data = Co_SPM, family=gaussian)

mix.int_Co_COL <- glm(Co ~ (Co_COL * TURNOVER), data = Co_COL, family=gaussian)

mix.int_Co_TRD <- glm(Co ~ (Co_TRD * TURNOVER), data = Co_TRD, family=gaussian)

plot(allEffects(mix.int_Co_SPM))

plot(allEffects(mix.int_Co_COL))

plot(allEffects(mix.int_Co_TRD))

Anova(mix.int_Co_SPM, test="LR")
summary(mix.int_Co_SPM)
with(summary(mix.int_Co_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Co_COL, test="LR")
summary(mix.int_Co_COL)
with(summary(mix.int_Co_COL), 1 - deviance/null.deviance)

Anova(mix.int_Co_TRD, test="LR")
summary(mix.int_Co_TRD)
with(summary(mix.int_Co_TRD), 1 - deviance/null.deviance)

### Co Col Plots ###

Co_p_SPM=ggplot(data=Co_SPM,aes(Co_SPM,Co))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Co_SPM, y=Co, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Co  (SPM, mg/l)",    
       y= "Co Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Co_p_COL=ggplot(data=Co_COL,aes(Co_COL,Co))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Co_COL, y=Co, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Co  (COL, mg/l)",    
       y= "Co Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Co_p_TRD=ggplot(data=Co_TRD,aes(Co_TRD,Co))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Co_TRD, y=Co, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Co  (TD, mg/l)",    
       y= "Co Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Co_p_SPM,Co_p_COL,Co_p_TRD,nrow=1)


#dev.new()
#### Cr ####

Cr_DF <- (COMPARTMENTS_AVG[,c("spm_Cr_NH3","colloidal_Cr_NH3","trulydissolved_Cr_NH3","TURNOVER","Cr")])

Cr_DF <- Cr_DF[complete.cases(Cr_DF),]

Cr_DF$TURNOVER <- asin(sqrt(Cr_DF$TURNOVER ))

quantile(Cr_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cr_DF$cat <- cut(Cr_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Cr_SPM ###

Cr_SPM <- (Cr_DF[,c(1,4,5,6)])

names(Cr_SPM)<- c("Cr_SPM", "TURNOVER", "Cr", "cat")

Cr_SPM<-Cr_SPM[Cr_SPM$Cr_SPM>0,]

Cr_SPM$Cr_SPM <- log10(Cr_SPM$Cr_SPM+1)

for (col in colnames(Cr_SPM[,c(1)])) {
  Cr_SPM[[col]] <- replace_outliers_with_na(Cr_SPM[[col]])
}

Cr_SPM <- Cr_SPM[complete.cases(Cr_SPM), ]

### Cr_COL ###

Cr_COL <- (Cr_DF[,c(2,4,5,6)])

names(Cr_COL)<- c("Cr_COL", "TURNOVER", "Cr", "cat")

Cr_COL<-Cr_COL[Cr_COL$Cr_COL>0,]

Cr_COL$Cr_COL <- log10(Cr_COL$Cr_COL+1)

for (col in colnames(Cr_COL[,c(1)])) {
  Cr_COL[[col]] <- replace_outliers_with_na(Cr_COL[[col]])
}

Cr_COL <- Cr_COL[complete.cases(Cr_COL), ]

### Cr_TRD ###

Cr_TRD <- (Cr_DF[,c(3,4,5,6)])

names(Cr_TRD)<- c("Cr_TRD", "TURNOVER", "Cr", "cat")

Cr_TRD<-Cr_TRD[Cr_TRD$Cr_TRD>0,]

Cr_TRD$Cr_TRD <- log10(Cr_TRD$Cr_TRD+1)

for (col in colnames(Cr_TRD[,c(1)])) {
  Cr_TRD[[col]] <- replace_outliers_with_na(Cr_TRD[[col]])
}

Cr_TRD <- Cr_TRD[complete.cases(Cr_TRD), ]

#dev.new()

plot(density(Cr_SPM$Cr_SPM, bw=0.0249), ylim=c(0,13), xlim=c(0,0.8), col="Purple", lwd=2)
lines(density(Cr_COL$Cr_COL, bw=0.0249), col="Blue", lwd=2)
lines(density(Cr_TRD$Cr_TRD, bw=0.0249), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cr_SPM <- glm(Cr ~ (Cr_SPM * TURNOVER), data = Cr_SPM, family=gaussian)

mix.int_Cr_COL <- glm(Cr ~ (Cr_COL * TURNOVER), data = Cr_COL, family=gaussian)

mix.int_Cr_TRD <- glm(Cr ~ (Cr_TRD * TURNOVER), data = Cr_TRD, family=gaussian)

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


#dev.new()
#### Cu ####

Cu_DF <- (COMPARTMENTS_AVG[,c("spm_Cu_NH3","colloidal_Cu_NH3","trulydissolved_Cu_NH3","TURNOVER","Cu")])

Cu_DF <- Cu_DF[complete.cases(Cu_DF),]

Cu_DF$TURNOVER <- asin(sqrt(Cu_DF$TURNOVER ))

quantile(Cu_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cu_DF$cat <- cut(Cu_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Cu_SPM ###

Cu_SPM <- (Cu_DF[,c(1,4,5,6)])

names(Cu_SPM)<- c("Cu_SPM", "TURNOVER", "Cu", "cat")

Cu_SPM<-Cu_SPM[Cu_SPM$Cu_SPM>0,]

Cu_SPM$Cu_SPM <- log10(Cu_SPM$Cu_SPM+1)

for (col in colnames(Cu_SPM[,c(1)])) {
  Cu_SPM[[col]] <- replace_outliers_with_na(Cu_SPM[[col]])
}

Cu_SPM <- Cu_SPM[complete.cases(Cu_SPM), ]

### Cu_COL ###

Cu_COL <- (Cu_DF[,c(2,4,5,6)])

names(Cu_COL)<- c("Cu_COL", "TURNOVER", "Cu", "cat")

Cu_COL<-Cu_COL[Cu_COL$Cu_COL>0,]

Cu_COL$Cu_COL <- log10(Cu_COL$Cu_COL+1)

for (col in colnames(Cu_COL[,c(1)])) {
  Cu_COL[[col]] <- replace_outliers_with_na(Cu_COL[[col]])
}

Cu_COL <- Cu_COL[complete.cases(Cu_COL), ]

### Cu_TRD ###

Cu_TRD <- (Cu_DF[,c(3,4,5,6)])

names(Cu_TRD)<- c("Cu_TRD", "TURNOVER", "Cu", "cat")

Cu_TRD<-Cu_TRD[Cu_TRD$Cu_TRD>0,]

Cu_TRD$Cu_TRD <- log10(Cu_TRD$Cu_TRD+1)

for (col in colnames(Cu_TRD[,c(1)])) {
  Cu_TRD[[col]] <- replace_outliers_with_na(Cu_TRD[[col]])
}

Cu_TRD <- Cu_TRD[complete.cases(Cu_TRD), ]

#dev.new()

plot(density(Cu_SPM$Cu_SPM, bw=0.09), ylim=c(0,3), xlim=c(0,1.5), col="Purple", lwd=2)
lines(density(Cu_COL$Cu_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Cu_TRD$Cu_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cu_SPM <- glm(Cu ~ (Cu_SPM * TURNOVER), data = Cu_SPM, family=gaussian)

mix.int_Cu_COL <- glm(Cu ~ (Cu_COL * TURNOVER), data = Cu_COL, family=gaussian)

mix.int_Cu_TRD <- glm(Cu ~ (Cu_TRD * TURNOVER), data = Cu_TRD, family=gaussian)

plot(allEffects(mix.int_Cu_SPM))

plot(allEffects(mix.int_Cu_COL))

plot(allEffects(mix.int_Cu_TRD))

Anova(mix.int_Cu_SPM, test="LR")
summary(mix.int_Cu_SPM)
with(summary(mix.int_Cu_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Cu_COL, test="LR")
summary(mix.int_Cu_COL)
with(summary(mix.int_Cu_COL), 1 - deviance/null.deviance)

Anova(mix.int_Cu_TRD, test="LR")
summary(mix.int_Cu_TRD)
with(summary(mix.int_Cu_TRD), 1 - deviance/null.deviance)

### Cu Cul Plots ###

Cu_p_SPM=ggplot(data=Cu_SPM,aes(Cu_SPM,Cu))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cu_SPM, y=Cu, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cu  (SPM, mg/l)",    
       y= "Cu Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cu_p_COL=ggplot(data=Cu_COL,aes(Cu_COL,Cu))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cu_COL, y=Cu, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cu  (COL, mg/l)",    
       y= "Cu Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cu_p_TRD=ggplot(data=Cu_TRD,aes(Cu_TRD,Cu))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cu_TRD, y=Cu, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cu  (TD, mg/l)",    
       y= "Cu Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Cu_p_SPM,Cu_p_COL,Cu_p_TRD,nrow=1)


#dev.new()

#### Fe ####

Fe_DF <- (COMPARTMENTS_AVG[,c("spm_Fe_NH3","colloidal_Fe_NH3","trulydissolved_Fe_NH3","TURNOVER","Fer.1")])

Fe_DF <- Fe_DF[complete.cases(Fe_DF),]

Fe_DF$TURNOVER <- asin(sqrt(Fe_DF$TURNOVER ))

quantile(Fe_DF$TURNOVER, probs = seq(0, 1, 1/5))

Fe_DF$cat <- cut(Fe_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Fe_SPM ###

Fe_SPM <- (Fe_DF[,c(1,4,5,6)])

names(Fe_SPM)<- c("Fe_SPM", "TURNOVER", "Fe", "cat")

Fe_SPM<-Fe_SPM[Fe_SPM$Fe_SPM>0,]

Fe_SPM$Fe_SPM <- log10(Fe_SPM$Fe_SPM+1)

for (col in colnames(Fe_SPM[,c(1)])) {
  Fe_SPM[[col]] <- replace_outliers_with_na(Fe_SPM[[col]])
}

Fe_SPM <- Fe_SPM[complete.cases(Fe_SPM), ]

### Fe_COL ###

Fe_COL <- (Fe_DF[,c(2,4,5,6)])

names(Fe_COL)<- c("Fe_COL", "TURNOVER", "Fe", "cat")

Fe_COL<-Fe_COL[Fe_COL$Fe_COL>0,]

Fe_COL$Fe_COL <- log10(Fe_COL$Fe_COL+1)

for (col in colnames(Fe_COL[,c(1)])) {
  Fe_COL[[col]] <- replace_outliers_with_na(Fe_COL[[col]])
}

Fe_COL <- Fe_COL[complete.cases(Fe_COL), ]

### Fe_TRD ###

Fe_TRD <- (Fe_DF[,c(3,4,5,6)])

names(Fe_TRD)<- c("Fe_TRD", "TURNOVER", "Fe", "cat")

Fe_TRD<-Fe_TRD[Fe_TRD$Fe_TRD>0,]

Fe_TRD$Fe_TRD <- log10(Fe_TRD$Fe_TRD+1)

for (col in colnames(Fe_TRD[,c(1)])) {
  Fe_TRD[[col]] <- replace_outliers_with_na(Fe_TRD[[col]])
}

Fe_TRD <- Fe_TRD[complete.cases(Fe_TRD), ]

#dev.new()

plot(density(Fe_SPM$Fe_SPM, bw=0.09), ylim=c(0,3), xlim=c(0,3), col="Purple", lwd=2)
lines(density(Fe_COL$Fe_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Fe_TRD$Fe_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Fe_SPM <- glm(Fe ~ (Fe_SPM * TURNOVER), data = Fe_SPM, family=gaussian)

mix.int_Fe_COL <- glm(Fe ~ (Fe_COL * TURNOVER), data = Fe_COL, family=gaussian)

mix.int_Fe_TRD <- glm(Fe ~ (Fe_TRD * TURNOVER), data = Fe_TRD, family=gaussian)

plot(allEffects(mix.int_Fe_SPM))

plot(allEffects(mix.int_Fe_COL))

plot(allEffects(mix.int_Fe_TRD))

Anova(mix.int_Fe_SPM, test="LR")
summary(mix.int_Fe_SPM)
with(summary(mix.int_Fe_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Fe_COL, test="LR")
summary(mix.int_Fe_COL)
with(summary(mix.int_Fe_COL), 1 - deviance/null.deviance)

Anova(mix.int_Fe_TRD, test="LR")
summary(mix.int_Fe_TRD)
with(summary(mix.int_Fe_TRD), 1 - deviance/null.deviance)

### Fe Fel Plots ###

Fe_p_SPM=ggplot(data=Fe_SPM,aes(Fe_SPM,Fe))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Fe_SPM, y=Fe, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Fe  (SPM, mg/l)",    
       y= "Fe Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Fe_p_COL=ggplot(data=Fe_COL,aes(Fe_COL,Fe))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Fe_COL, y=Fe, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Fe  (COL, mg/l)",    
       y= "Fe Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Fe_p_TRD=ggplot(data=Fe_TRD,aes(Fe_TRD,Fe))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Fe_TRD, y=Fe, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Fe  (TD, mg/l)",    
       y= "Fe Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Fe_p_SPM,Fe_p_COL,Fe_p_TRD,nrow=1)


#dev.new()
#### Mo ####

Mo_DF <- (COMPARTMENTS_AVG[,c("spm_Mo_Oshift","colloidal_Mo_Oshift","trulydissolved_Mo_Oshift","TURNOVER","Mo")])

Mo_DF <- Mo_DF[complete.cases(Mo_DF),]

Mo_DF$TURNOVER <- asin(sqrt(Mo_DF$TURNOVER ))

quantile(Mo_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mo_DF$cat <- cut(Mo_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Mo_SPM ###

Mo_SPM <- (Mo_DF[,c(1,4,5,6)])

names(Mo_SPM)<- c("Mo_SPM", "TURNOVER", "Mo", "cat")

Mo_SPM<-Mo_SPM[Mo_SPM$Mo_SPM>0,]

Mo_SPM$Mo_SPM <- log10(Mo_SPM$Mo_SPM+1)

for (col in colnames(Mo_SPM[,c(1)])) {
  Mo_SPM[[col]] <- replace_outliers_with_na(Mo_SPM[[col]])
}

Mo_SPM <- Mo_SPM[complete.cases(Mo_SPM), ]

### Mo_COL ###

Mo_COL <- (Mo_DF[,c(2,4,5,6)])

names(Mo_COL)<- c("Mo_COL", "TURNOVER", "Mo", "cat")

Mo_COL<-Mo_COL[Mo_COL$Mo_COL>0,]

Mo_COL$Mo_COL <- log10(Mo_COL$Mo_COL+1)

for (col in colnames(Mo_COL[,c(1)])) {
  Mo_COL[[col]] <- replace_outliers_with_na(Mo_COL[[col]])
}

Mo_COL <- Mo_COL[complete.cases(Mo_COL), ]

### Mo_TRD ###

Mo_TRD <- (Mo_DF[,c(3,4,5,6)])

names(Mo_TRD)<- c("Mo_TRD", "TURNOVER", "Mo", "cat")

Mo_TRD<-Mo_TRD[Mo_TRD$Mo_TRD>0,]

Mo_TRD$Mo_TRD <- log10(Mo_TRD$Mo_TRD+1)

for (col in colnames(Mo_TRD[,c(1)])) {
  Mo_TRD[[col]] <- replace_outliers_with_na(Mo_TRD[[col]])
}

Mo_TRD <- Mo_TRD[complete.cases(Mo_TRD), ]

#dev.new()

plot(density(Mo_SPM$Mo_SPM, bw=0.1079), ylim=c(0,3), xlim=c(0,1.8), col="Purple", lwd=2)
lines(density(Mo_COL$Mo_COL, bw=0.1079), col="Blue", lwd=2)
lines(density(Mo_TRD$Mo_TRD, bw=0.1079), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Mo_SPM <- glm(Mo ~ (Mo_SPM * TURNOVER), data = Mo_SPM, family=gaussian)

mix.int_Mo_COL <- glm(Mo ~ (Mo_COL * TURNOVER), data = Mo_COL, family=gaussian)

mix.int_Mo_TRD <- glm(Mo ~ (Mo_TRD * TURNOVER), data = Mo_TRD, family=gaussian)

plot(allEffects(mix.int_Mo_SPM))

plot(allEffects(mix.int_Mo_COL))

plot(allEffects(mix.int_Mo_TRD))

Anova(mix.int_Mo_SPM, test="LR")
summary(mix.int_Mo_SPM)
with(summary(mix.int_Mo_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Mo_COL, test="LR")
summary(mix.int_Mo_COL)
with(summary(mix.int_Mo_COL), 1 - deviance/null.deviance)

Anova(mix.int_Mo_TRD, test="LR")
summary(mix.int_Mo_TRD)
with(summary(mix.int_Mo_TRD), 1 - deviance/null.deviance)

### Mo Mol Plots ###

Mo_p_SPM=ggplot(data=Mo_SPM,aes(Mo_SPM,Mo))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mo_SPM, y=Mo, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mo  (SPM, mg/l)",    
       y= "Mo Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mo_p_COL=ggplot(data=Mo_COL,aes(Mo_COL,Mo))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mo_COL, y=Mo, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mo  (COL, mg/l)",    
       y= "Mo Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mo_p_TRD=ggplot(data=Mo_TRD,aes(Mo_TRD,Mo))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mo_TRD, y=Mo, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mo  (TD, mg/l)",    
       y= "Mo Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Mo_p_SPM,Mo_p_COL,Mo_p_TRD,nrow=1)


#dev.new()

#### Mn ####

Mn_DF <- (COMPARTMENTS_AVG[,c("spm_Mn","colloidal_Mn","trulydissolved_Mn","TURNOVER","Mn")])

Mn_DF <- Mn_DF[complete.cases(Mn_DF),]

Mn_DF$TURNOVER <- asin(sqrt(Mn_DF$TURNOVER ))

quantile(Mn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mn_DF$cat <- cut(Mn_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Mn_SPM ###

Mn_SPM <- (Mn_DF[,c(1,4,5,6)])

names(Mn_SPM)<- c("Mn_SPM", "TURNOVER", "Mn", "cat")

Mn_SPM<-Mn_SPM[Mn_SPM$Mn_SPM>0,]

Mn_SPM$Mn_SPM <- log10(Mn_SPM$Mn_SPM+1)

for (col in colnames(Mn_SPM[,c(1)])) {
  Mn_SPM[[col]] <- replace_outliers_with_na(Mn_SPM[[col]])
}

Mn_SPM <- Mn_SPM[complete.cases(Mn_SPM), ]

### Mn_COL ###

Mn_COL <- (Mn_DF[,c(2,4,5,6)])

names(Mn_COL)<- c("Mn_COL", "TURNOVER", "Mn", "cat")

Mn_COL<-Mn_COL[Mn_COL$Mn_COL>0,]

Mn_COL$Mn_COL <- log10(Mn_COL$Mn_COL+1)

for (col in colnames(Mn_COL[,c(1)])) {
  Mn_COL[[col]] <- replace_outliers_with_na(Mn_COL[[col]])
}

Mn_COL <- Mn_COL[complete.cases(Mn_COL), ]

### Mn_TRD ###

Mn_TRD <- (Mn_DF[,c(3,4,5,6)])

names(Mn_TRD)<- c("Mn_TRD", "TURNOVER", "Mn", "cat")

Mn_TRD<-Mn_TRD[Mn_TRD$Mn_TRD>0,]

Mn_TRD$Mn_TRD <- log10(Mn_TRD$Mn_TRD+1)

for (col in colnames(Mn_TRD[,c(1)])) {
  Mn_TRD[[col]] <- replace_outliers_with_na(Mn_TRD[[col]])
}

Mn_TRD <- Mn_TRD[complete.cases(Mn_TRD), ]

#dev.new()

plot(density(Mn_SPM$Mn_SPM, bw=.1205), ylim=c(0,2), xlim=c(0,2.4), col="Purple", lwd=2)
lines(density(Mn_COL$Mn_COL, bw=.1205), col="Blue", lwd=2)
lines(density(Mn_TRD$Mn_TRD, bw=.1205), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Mn_SPM <- glm(Mn ~ (Mn_SPM * TURNOVER), data = Mn_SPM, family=gaussian)

mix.int_Mn_COL <- glm(Mn ~ (Mn_COL * TURNOVER), data = Mn_COL, family=gaussian)

mix.int_Mn_TRD <- glm(Mn ~ (Mn_TRD * TURNOVER), data = Mn_TRD, family=gaussian)

plot(allEffects(mix.int_Mn_SPM))

plot(allEffects(mix.int_Mn_COL))

plot(allEffects(mix.int_Mn_TRD))

Anova(mix.int_Mn_SPM, test="LR")
summary(mix.int_Mn_SPM)
with(summary(mix.int_Mn_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Mn_COL, test="LR")
summary(mix.int_Mn_COL)
with(summary(mix.int_Mn_COL), 1 - deviance/null.deviance)

Anova(mix.int_Mn_TRD, test="LR")
summary(mix.int_Mn_TRD)
with(summary(mix.int_Mn_TRD), 1 - deviance/null.deviance)

### Mn Mnl Plots ###

Mn_p_SPM=ggplot(data=Mn_SPM,aes(Mn_SPM,Mn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mn_SPM, y=Mn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mn  (SPM, mg/l)",    
       y= "Mn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mn_p_COL=ggplot(data=Mn_COL,aes(Mn_COL,Mn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mn_COL, y=Mn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mn  (COL, mg/l)",    
       y= "Mn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mn_p_TRD=ggplot(data=Mn_TRD,aes(Mn_TRD,Mn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Mn_TRD, y=Mn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mn  (TD, mg/l)",    
       y= "Mn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Mn_p_SPM,Mn_p_COL,Mn_p_TRD,nrow=1)


#### Ni ####

Ni_DF <- (COMPARTMENTS_AVG[,c("spm_Ni","colloidal_Ni","trulydissolved_Ni","TURNOVER","Ni")])

Ni_DF <- Ni_DF[complete.cases(Ni_DF),]

Ni_DF$TURNOVER <- asin(sqrt(Ni_DF$TURNOVER ))

quantile(Ni_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ni_DF$cat <- cut(Ni_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Ni_SPM ###

Ni_SPM <- (Ni_DF[,c(1,4,5,6)])

names(Ni_SPM)<- c("Ni_SPM", "TURNOVER", "Ni", "cat")

Ni_SPM<-Ni_SPM[Ni_SPM$Ni_SPM>0,]

Ni_SPM$Ni_SPM <- log10(Ni_SPM$Ni_SPM+1)

for (col in colnames(Ni_SPM[,c(1)])) {
  Ni_SPM[[col]] <- replace_outliers_with_na(Ni_SPM[[col]])
}

Ni_SPM <- Ni_SPM[complete.cases(Ni_SPM), ]

### Ni_COL ###

Ni_COL <- (Ni_DF[,c(2,4,5,6)])

names(Ni_COL)<- c("Ni_COL", "TURNOVER", "Ni", "cat")

Ni_COL<-Ni_COL[Ni_COL$Ni_COL>0,]

Ni_COL$Ni_COL <- log10(Ni_COL$Ni_COL+1)

for (col in colnames(Ni_COL[,c(1)])) {
  Ni_COL[[col]] <- replace_outliers_with_na(Ni_COL[[col]])
}

Ni_COL <- Ni_COL[complete.cases(Ni_COL), ]

### Ni_TRD ###

Ni_TRD <- (Ni_DF[,c(3,4,5,6)])

names(Ni_TRD)<- c("Ni_TRD", "TURNOVER", "Ni", "cat")

Ni_TRD<-Ni_TRD[Ni_TRD$Ni_TRD>0,]

Ni_TRD$Ni_TRD <- log10(Ni_TRD$Ni_TRD+1)

for (col in colnames(Ni_TRD[,c(1)])) {
  Ni_TRD[[col]] <- replace_outliers_with_na(Ni_TRD[[col]])
}

Ni_TRD <- Ni_TRD[complete.cases(Ni_TRD), ]

#dev.new()

plot(density(Ni_SPM$Ni_SPM, bw=0.0029), ylim=c(0,150), xlim=c(0,0.04), col="Purple", lwd=2)
lines(density(Ni_COL$Ni_COL, bw=0.0029), col="Blue", lwd=2)
lines(density(Ni_TRD$Ni_TRD, bw=0.0029), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Ni_SPM <- glm(Ni ~ (Ni_SPM * TURNOVER), data = Ni_SPM, family=gaussian)

mix.int_Ni_COL <- glm(Ni ~ (Ni_COL * TURNOVER), data = Ni_COL, family=gaussian)

mix.int_Ni_TRD <- glm(Ni ~ (Ni_TRD * TURNOVER), data = Ni_TRD, family=gaussian)

plot(allEffects(mix.int_Ni_SPM))

plot(allEffects(mix.int_Ni_COL))

plot(allEffects(mix.int_Ni_TRD))

Anova(mix.int_Ni_SPM, test="LR")
summary(mix.int_Ni_SPM)
with(summary(mix.int_Ni_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Ni_COL, test="LR")
summary(mix.int_Ni_COL)
with(summary(mix.int_Ni_COL), 1 - deviance/null.deviance)

Anova(mix.int_Ni_TRD, test="LR")
summary(mix.int_Ni_TRD)
with(summary(mix.int_Ni_TRD), 1 - deviance/null.deviance)

### Ni Nil Plots ###

Ni_p_SPM=ggplot(data=Ni_SPM,aes(Ni_SPM,Ni))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ni_SPM, y=Ni, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ni  (SPM, mg/l)",    
       y= "Ni Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ni_p_COL=ggplot(data=Ni_COL,aes(Ni_COL,Ni))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ni_COL, y=Ni, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ni  (COL, mg/l)",    
       y= "Ni Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ni_p_TRD=ggplot(data=Ni_TRD,aes(Ni_TRD,Ni))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Ni_TRD, y=Ni, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ni  (TD, mg/l)",    
       y= "Ni Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Ni_p_SPM,Ni_p_COL,Ni_p_TRD,nrow=1)


#dev.new()
#### P ####

P_DF <- (COMPARTMENTS_AVG[,c("spm_P","colloidal_P","trulydissolved_P","TURNOVER","P")])

P_DF <- P_DF[complete.cases(P_DF),]

P_DF$TURNOVER <- asin(sqrt(P_DF$TURNOVER ))

quantile(P_DF$TURNOVER, probs = seq(0, 1, 1/5))

P_DF$cat <- cut(P_DF$TURNOVER, 
                breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### P_SPM ###

P_SPM <- (P_DF[,c(1,4,5,6)])

names(P_SPM)<- c("P_SPM", "TURNOVER", "P", "cat")

P_SPM<-P_SPM[P_SPM$P_SPM>0,]

P_SPM$P_SPM <- log10(P_SPM$P_SPM+1)

for (col in colnames(P_SPM[,c(1)])) {
  P_SPM[[col]] <- replace_outliers_with_na(P_SPM[[col]])
}

P_SPM <- P_SPM[complete.cases(P_SPM), ]

### P_COL ###

P_COL <- (P_DF[,c(2,4,5,6)])

names(P_COL)<- c("P_COL", "TURNOVER", "P", "cat")

P_COL<-P_COL[P_COL$P_COL>0,]

P_COL$P_COL <- log10(P_COL$P_COL+1)

for (col in colnames(P_COL[,c(1)])) {
  P_COL[[col]] <- replace_outliers_with_na(P_COL[[col]])
}

P_COL <- P_COL[complete.cases(P_COL), ]

### P_TRD ###

P_TRD <- (P_DF[,c(3,4,5,6)])

names(P_TRD)<- c("P_TRD", "TURNOVER", "P", "cat")

P_TRD<-P_TRD[P_TRD$P_TRD>0,]

P_TRD$P_TRD <- log10(P_TRD$P_TRD+1)

for (col in colnames(P_TRD[,c(1)])) {
  P_TRD[[col]] <- replace_outliers_with_na(P_TRD[[col]])
}

P_TRD <- P_TRD[complete.cases(P_TRD), ]

#dev.new()

plot(density(P_SPM$P_SPM, bw=0.1637), ylim=c(0,1.2), xlim=c(0,2.1), col="Purple", lwd=2)
lines(density(P_COL$P_COL, bw=0.1637), col="Blue", lwd=2)
lines(density(P_TRD$P_TRD, bw=0.1637), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_P_SPM <- glm(P ~ (P_SPM * TURNOVER), data = P_SPM, family=gaussian)

mix.int_P_COL <- glm(P ~ (P_COL * TURNOVER), data = P_COL, family=gaussian)

mix.int_P_TRD <- glm(P ~ (P_TRD * TURNOVER), data = P_TRD, family=gaussian)

plot(allEffects(mix.int_P_SPM))

plot(allEffects(mix.int_P_COL))

plot(allEffects(mix.int_P_TRD))

Anova(mix.int_P_SPM, test="LR")
summary(mix.int_P_SPM)
with(summary(mix.int_P_SPM), 1 - deviance/null.deviance)

Anova(mix.int_P_COL, test="LR")
summary(mix.int_P_COL)
with(summary(mix.int_P_COL), 1 - deviance/null.deviance)

Anova(mix.int_P_TRD, test="LR")
summary(mix.int_P_TRD)
with(summary(mix.int_P_TRD), 1 - deviance/null.deviance)

### P Pl Plots ###

P_p_SPM=ggplot(data=P_SPM,aes(P_SPM,P))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=P_SPM, y=P, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "P  (SPM, mg/l)",    
       y= "P Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

P_p_COL=ggplot(data=P_COL,aes(P_COL,P))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=P_COL, y=P, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "P  (COL, mg/l)",    
       y= "P Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

P_p_TRD=ggplot(data=P_TRD,aes(P_TRD,P))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=P_TRD, y=P, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "P  (TD, mg/l)",    
       y= "P Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(P_p_SPM,P_p_COL,P_p_TRD,nrow=1)


#dev.new()

#### Pb ####

Pb_DF <- (COMPARTMENTS_AVG[,c("spm_Pb_NH3","colloidal_Pb_NH3","trulydissolved_Pb_NH3","TURNOVER","Pb")])

Pb_DF <- Pb_DF[complete.cases(Pb_DF),]

Pb_DF$TURNOVER <- asin(sqrt(Pb_DF$TURNOVER ))

quantile(Pb_DF$TURNOVER, probs = seq(0, 1, 1/5))

Pb_DF$cat <- cut(Pb_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Pb_SPM ###

Pb_SPM <- (Pb_DF[,c(1,4,5,6)])

names(Pb_SPM)<- c("Pb_SPM", "TURNOVER", "Pb", "cat")

Pb_SPM<-Pb_SPM[Pb_SPM$Pb_SPM>0,]

Pb_SPM$Pb_SPM <- log10(Pb_SPM$Pb_SPM+1)

for (col in colnames(Pb_SPM[,c(1)])) {
  Pb_SPM[[col]] <- replace_outliers_with_na(Pb_SPM[[col]])
}

Pb_SPM <- Pb_SPM[complete.cases(Pb_SPM), ]

### Pb_COL ###

Pb_COL <- (Pb_DF[,c(2,4,5,6)])

names(Pb_COL)<- c("Pb_COL", "TURNOVER", "Pb", "cat")

Pb_COL<-Pb_COL[Pb_COL$Pb_COL>0,]

Pb_COL$Pb_COL <- log10(Pb_COL$Pb_COL+1)

for (col in colnames(Pb_COL[,c(1)])) {
  Pb_COL[[col]] <- replace_outliers_with_na(Pb_COL[[col]])
}

Pb_COL <- Pb_COL[complete.cases(Pb_COL), ]

### Pb_TRD ###

Pb_TRD <- (Pb_DF[,c(3,4,5,6)])

names(Pb_TRD)<- c("Pb_TRD", "TURNOVER", "Pb", "cat")

Pb_TRD<-Pb_TRD[Pb_TRD$Pb_TRD>0,]

Pb_TRD$Pb_TRD <- log10(Pb_TRD$Pb_TRD+1)

for (col in colnames(Pb_TRD[,c(1)])) {
  Pb_TRD[[col]] <- replace_outliers_with_na(Pb_TRD[[col]])
}

Pb_TRD <- Pb_TRD[complete.cases(Pb_TRD), ]

#dev.new()

plot(density(Pb_SPM$Pb_SPM, bw=0.05112), ylim=c(0,7), xlim=c(0,0.7), col="Purple", lwd=2)
lines(density(Pb_COL$Pb_COL, bw=0.05112), col="Blue", lwd=2)
lines(density(Pb_TRD$Pb_TRD, bw=0.05112), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Pb_SPM <- glm(Pb ~ (Pb_SPM * TURNOVER), data = Pb_SPM, family=gaussian)

mix.int_Pb_COL <- glm(Pb ~ (Pb_COL * TURNOVER), data = Pb_COL, family=gaussian)

mix.int_Pb_TRD <- glm(Pb ~ (Pb_TRD * TURNOVER), data = Pb_TRD, family=gaussian)

plot(allEffects(mix.int_Pb_SPM))

plot(allEffects(mix.int_Pb_COL))

plot(allEffects(mix.int_Pb_TRD))

Anova(mix.int_Pb_SPM, test="LR")
summary(mix.int_Pb_SPM)
with(summary(mix.int_Pb_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Pb_COL, test="LR")
summary(mix.int_Pb_COL)
with(summary(mix.int_Pb_COL), 1 - deviance/null.deviance)

Anova(mix.int_Pb_TRD, test="LR")
summary(mix.int_Pb_TRD)
with(summary(mix.int_Pb_TRD), 1 - deviance/null.deviance)

### Pb Pbl Plots ###

Pb_p_SPM=ggplot(data=Pb_SPM,aes(Pb_SPM,Pb))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Pb_SPM, y=Pb, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Pb  (SPM, mg/l)",    
       y= "Pb Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Pb_p_COL=ggplot(data=Pb_COL,aes(Pb_COL,Pb))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Pb_COL, y=Pb, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Pb  (COL, mg/l)",    
       y= "Pb Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Pb_p_TRD=ggplot(data=Pb_TRD,aes(Pb_TRD,Pb))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Pb_TRD, y=Pb, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Pb  (TD, mg/l)",    
       y= "Pb Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Pb_p_SPM,Pb_p_COL,Pb_p_TRD,nrow=1)


#dev.new()

#### S ####

S_DF <- (COMPARTMENTS_AVG[,c("spm_S","colloidal_S","trulydissolved_S","TURNOVER","S")])

S_DF <- S_DF[complete.cases(S_DF),]

S_DF$TURNOVER <- asin(sqrt(S_DF$TURNOVER ))

quantile(S_DF$TURNOVER, probs = seq(0, 1, 1/5))

S_DF$cat <- cut(S_DF$TURNOVER, 
                breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### S_SPM ###

S_SPM <- (S_DF[,c(1,4,5,6)])

names(S_SPM)<- c("S_SPM", "TURNOVER", "S", "cat")

S_SPM<-S_SPM[S_SPM$S_SPM>0,]

S_SPM$S_SPM <- log10(S_SPM$S_SPM+1)

for (col in colnames(S_SPM[,c(1)])) {
  S_SPM[[col]] <- replace_outliers_with_na(S_SPM[[col]])
}

S_SPM <- S_SPM[complete.cases(S_SPM), ]

### S_COL ###

S_COL <- (S_DF[,c(2,4,5,6)])

names(S_COL)<- c("S_COL", "TURNOVER", "S", "cat")

S_COL<-S_COL[S_COL$S_COL>0,]

S_COL$S_COL <- log10(S_COL$S_COL+1)

for (col in colnames(S_COL[,c(1)])) {
  S_COL[[col]] <- replace_outliers_with_na(S_COL[[col]])
}

S_COL <- S_COL[complete.cases(S_COL), ]

### S_TRD ###

S_TRD <- (S_DF[,c(3,4,5,6)])

names(S_TRD)<- c("S_TRD", "TURNOVER", "S", "cat")

S_TRD<-S_TRD[S_TRD$S_TRD>0,]

S_TRD$S_TRD <- log10(S_TRD$S_TRD+1)

for (col in colnames(S_TRD[,c(1)])) {
  S_TRD[[col]] <- replace_outliers_with_na(S_TRD[[col]])
}

S_TRD <- S_TRD[complete.cases(S_TRD), ]

#dev.new()

plot(density(S_SPM$S_SPM, bw=0.1565), ylim=c(0,2.5), xlim=c(2,5), col="Purple", lwd=2)
lines(density(S_COL$S_COL, bw=0.1565), col="Blue", lwd=2)
lines(density(S_TRD$S_TRD, bw=0.1565), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_S_SPM <- glm(S ~ (S_SPM * TURNOVER), data = S_SPM, family=gaussian)

mix.int_S_COL <- glm(S ~ (S_COL * TURNOVER), data = S_COL, family=gaussian)

mix.int_S_TRD <- glm(S ~ (S_TRD * TURNOVER), data = S_TRD, family=gaussian)

plot(allEffects(mix.int_S_SPM))

plot(allEffects(mix.int_S_COL))

plot(allEffects(mix.int_S_TRD))

Anova(mix.int_S_SPM, test="LR")
summary(mix.int_S_SPM)
with(summary(mix.int_S_SPM), 1 - deviance/null.deviance)

Anova(mix.int_S_COL, test="LR")
summary(mix.int_S_COL)
with(summary(mix.int_S_COL), 1 - deviance/null.deviance)

Anova(mix.int_S_TRD, test="LR")
summary(mix.int_S_TRD)
with(summary(mix.int_S_TRD), 1 - deviance/null.deviance)

### S Sl Plots ###

S_p_SPM=ggplot(data=S_SPM,aes(S_SPM,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_SPM, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (SPM, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

S_p_COL=ggplot(data=S_COL,aes(S_COL,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_COL, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (COL, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

S_p_TRD=ggplot(data=S_TRD,aes(S_TRD,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_TRD, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (TD, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(S_p_SPM,S_p_COL,S_p_TRD,nrow=1)


#dev.new()
#### Se ####

Se_DF <- (COMPARTMENTS_AVG[,c("spm_Se82_Oshift","colloidal_Se82_Oshift","trulydissolved_Se82_Oshift","TURNOVER","Se")])

Se_DF <- Se_DF[complete.cases(Se_DF),]

Se_DF$TURNOVER <- asin(sqrt(Se_DF$TURNOVER ))

quantile(Se_DF$TURNOVER, probs = seq(0, 1, 1/5))

Se_DF$cat <- cut(Se_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Se_SPM ###

Se_SPM <- (Se_DF[,c(1,4,5,6)])

names(Se_SPM)<- c("Se_SPM", "TURNOVER", "Se", "cat")

Se_SPM<-Se_SPM[Se_SPM$Se_SPM>0,]

Se_SPM$Se_SPM <- log10(Se_SPM$Se_SPM+1)

for (col in colnames(Se_SPM[,c(1)])) {
  Se_SPM[[col]] <- replace_outliers_with_na(Se_SPM[[col]])
}

Se_SPM <- Se_SPM[complete.cases(Se_SPM), ]

### Se_COL ###

Se_COL <- (Se_DF[,c(2,4,5,6)])

names(Se_COL)<- c("Se_COL", "TURNOVER", "Se", "cat")

Se_COL<-Se_COL[Se_COL$Se_COL>0,]

Se_COL$Se_COL <- log10(Se_COL$Se_COL+1)

for (col in colnames(Se_COL[,c(1)])) {
  Se_COL[[col]] <- replace_outliers_with_na(Se_COL[[col]])
}

Se_COL <- Se_COL[complete.cases(Se_COL), ]

### Se_TRD ###

Se_TRD <- (Se_DF[,c(3,4,5,6)])

names(Se_TRD)<- c("Se_TRD", "TURNOVER", "Se", "cat")

Se_TRD<-Se_TRD[Se_TRD$Se_TRD>0,]

Se_TRD$Se_TRD <- log10(Se_TRD$Se_TRD+1)

for (col in colnames(Se_TRD[,c(1)])) {
  Se_TRD[[col]] <- replace_outliers_with_na(Se_TRD[[col]])
}

Se_TRD <- Se_TRD[complete.cases(Se_TRD), ]

#dev.new()

plot(density(Se_SPM$Se_SPM, bw=0.01358), ylim=c(0,20), xlim=c(0,0.20), col="Purple", lwd=2)
lines(density(Se_COL$Se_COL, bw=0.01358), col="Blue", lwd=2)
lines(density(Se_TRD$Se_TRD, bw=0.01358), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Se_SPM <- glm(Se ~ (Se_SPM * TURNOVER), data = Se_SPM, family=gaussian)

mix.int_Se_COL <- glm(Se ~ (Se_COL * TURNOVER), data = Se_COL, family=gaussian)

mix.int_Se_TRD <- glm(Se ~ (Se_TRD * TURNOVER), data = Se_TRD, family=gaussian)

plot(allEffects(mix.int_Se_SPM))

plot(allEffects(mix.int_Se_COL))

plot(allEffects(mix.int_Se_TRD))

Anova(mix.int_Se_SPM, test="LR")
summary(mix.int_Se_SPM)
with(summary(mix.int_Se_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Se_COL, test="LR")
summary(mix.int_Se_COL)
with(summary(mix.int_Se_COL), 1 - deviance/null.deviance)

Anova(mix.int_Se_TRD, test="LR")
summary(mix.int_Se_TRD)
with(summary(mix.int_Se_TRD), 1 - deviance/null.deviance)

### Se Sel Plots ###

Se_p_SPM=ggplot(data=Se_SPM,aes(Se_SPM,Se))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Se_SPM, y=Se, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Se  (SPM, mg/l)",    
       y= "Se Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Se_p_COL=ggplot(data=Se_COL,aes(Se_COL,Se))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Se_COL, y=Se, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Se  (COL, mg/l)",    
       y= "Se Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Se_p_TRD=ggplot(data=Se_TRD,aes(Se_TRD,Se))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Se_TRD, y=Se, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Se  (TD, mg/l)",    
       y= "Se Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Se_p_SPM,Se_p_COL,Se_p_TRD,nrow=1)


#dev.new()
#### Zn ####

Zn_DF <- (COMPARTMENTS_AVG[,c("spm_Zn","colloidal_Zn","trulydissolved_Zn","TURNOVER","Zn")])

Zn_DF <- Zn_DF[complete.cases(Zn_DF),]

Zn_DF$TURNOVER <- asin(sqrt(Zn_DF$TURNOVER ))

quantile(Zn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Zn_DF$cat <- cut(Zn_DF$TURNOVER, 
                 breaks=c(-Inf, 0.1710476    , 0.2366313     , 0.2471030   , 0.2728182 , 0.2933796   , Inf), 
                 labels=c("NA", "0.18 - 0.24","0.24 - 0.25","0.25 - 0.27","0.27 - 0.29","0.29 - 0.34"))

### Zn_SPM ###

Zn_SPM <- (Zn_DF[,c(1,4,5,6)])

names(Zn_SPM)<- c("Zn_SPM", "TURNOVER", "Zn", "cat")

Zn_SPM<-Zn_SPM[Zn_SPM$Zn_SPM>0,]

Zn_SPM$Zn_SPM <- log10(Zn_SPM$Zn_SPM+1)

for (col in colnames(Zn_SPM[,c(1)])) {
  Zn_SPM[[col]] <- replace_outliers_with_na(Zn_SPM[[col]])
}

Zn_SPM <- Zn_SPM[complete.cases(Zn_SPM), ]

### Zn_COL ###

Zn_COL <- (Zn_DF[,c(2,4,5,6)])

names(Zn_COL)<- c("Zn_COL", "TURNOVER", "Zn", "cat")

Zn_COL<-Zn_COL[Zn_COL$Zn_COL>0,]

Zn_COL$Zn_COL <- log10(Zn_COL$Zn_COL+1)

for (col in colnames(Zn_COL[,c(1)])) {
  Zn_COL[[col]] <- replace_outliers_with_na(Zn_COL[[col]])
}

Zn_COL <- Zn_COL[complete.cases(Zn_COL), ]

### Zn_TRD ###

Zn_TRD <- (Zn_DF[,c(3,4,5,6)])

names(Zn_TRD)<- c("Zn_TRD", "TURNOVER", "Zn", "cat")

Zn_TRD<-Zn_TRD[Zn_TRD$Zn_TRD>0,]

Zn_TRD$Zn_TRD <- log10(Zn_TRD$Zn_TRD+1)

for (col in colnames(Zn_TRD[,c(1)])) {
  Zn_TRD[[col]] <- replace_outliers_with_na(Zn_TRD[[col]])
}

Zn_TRD <- Zn_TRD[complete.cases(Zn_TRD), ]

#dev.new()

plot(density(Zn_SPM$Zn_SPM, bw=0.1461), ylim=c(0,2), xlim=c(0,2), col="Purple", lwd=2)
lines(density(Zn_COL$Zn_COL, bw=0.1461), col="Blue", lwd=2)
lines(density(Zn_TRD$Zn_TRD, bw=0.1461), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Zn_SPM <- glm(Zn ~ (Zn_SPM * TURNOVER), data = Zn_SPM, family=gaussian)

mix.int_Zn_COL <- glm(Zn ~ (Zn_COL * TURNOVER), data = Zn_COL, family=gaussian)

mix.int_Zn_TRD <- glm(Zn ~ (Zn_TRD * TURNOVER), data = Zn_TRD, family=gaussian)

plot(allEffects(mix.int_Zn_SPM))

plot(allEffects(mix.int_Zn_COL))

plot(allEffects(mix.int_Zn_TRD))

Anova(mix.int_Zn_SPM, test="LR")
summary(mix.int_Zn_SPM)
with(summary(mix.int_Zn_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Zn_COL, test="LR")
summary(mix.int_Zn_COL)
with(summary(mix.int_Zn_COL), 1 - deviance/null.deviance)

Anova(mix.int_Zn_TRD, test="LR")
summary(mix.int_Zn_TRD)
with(summary(mix.int_Zn_TRD), 1 - deviance/null.deviance)

### Zn Znl Plots ###

Zn_p_SPM=ggplot(data=Zn_SPM,aes(Zn_SPM,Zn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Zn_SPM, y=Zn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Zn  (SPM, mg/l)",    
       y= "Zn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Zn_p_COL=ggplot(data=Zn_COL,aes(Zn_COL,Zn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Zn_COL, y=Zn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Zn  (COL, mg/l)",    
       y= "Zn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Zn_p_TRD=ggplot(data=Zn_TRD,aes(Zn_TRD,Zn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Zn_TRD, y=Zn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Zn  (TD, mg/l)",    
       y= "Zn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Zn_p_SPM,Zn_p_COL,Zn_p_TRD,nrow=1)

#### % Total ####

names(Al_DF)<-c("Al_SPM","Al_COL","Al_TRD","TURNOVER", "Al", "cat")
Al_DF[,c(1:3)][Al_DF[,c(1:3)] <= 0] <- 0

sum(Al_DF$Al_SPM, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))
sum(Al_DF$Al_COL, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))
sum(Al_DF$Al_TRD, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))

names(As_DF)<-c("As_SPM","As_COL","As_TRD","TURNOVER", "As", "cat")
As_DF[,c(1:3)][As_DF[,c(1:3)] <= 0] <- 0

sum(As_DF$As_SPM, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_COL, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))
sum(As_DF$As_TRD, na.rm=T)/(sum(As_DF$As_SPM,As_DF$As_COL,As_DF$As_TRD, na.rm=T))

names(Ca_DF)<-c("Ca_SPM","Ca_COL","Ca_TRD","TURNOVER", "Ca", "cat")
Ca_DF[,c(1:3)][Ca_DF[,c(1:3)] <= 0] <- 0

sum(Ca_DF$Ca_SPM, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))
sum(Ca_DF$Ca_COL, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))
sum(Ca_DF$Ca_TRD, na.rm=T)/(sum(Ca_DF$Ca_SPM,Ca_DF$Ca_COL,Ca_DF$Ca_TRD, na.rm=T))

names(Cd_DF)<-c("Cd_SPM","Cd_COL","Cd_TRD","TURNOVER", "Cd", "cat")
Cd_DF[,c(1:3)][Cd_DF[,c(1:3)] <= 0] <- 0

sum(Cd_DF$Cd_SPM, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_COL, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))
sum(Cd_DF$Cd_TRD, na.rm=T)/(sum(Cd_DF$Cd_SPM,Cd_DF$Cd_COL,Cd_DF$Cd_TRD, na.rm=T))

names(Cr_DF)<-c("Cr_SPM","Cr_COL","Cr_TRD","TURNOVER", "Cr", "cat")
Cr_DF[,c(1:3)][Cr_DF[,c(1:3)] <= 0] <- 0

sum(Cr_DF$Cr_SPM, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))
sum(Cr_DF$Cr_COL, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))
sum(Cr_DF$Cr_TRD, na.rm=T)/(sum(Cr_DF$Cr_SPM,Cr_DF$Cr_COL,Cr_DF$Cr_TRD, na.rm=T))

names(Co_DF)<-c("Co_SPM","Co_COL","Co_TRD","TURNOVER", "Co", "cat")
Co_DF[,c(1:3)][Co_DF[,c(1:3)] <= 0] <- 0

sum(Co_DF$Co_SPM, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_COL, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_TRD, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))

names(Cu_DF)<-c("Cu_SPM","Cu_COL","Cu_TRD","TURNOVER", "Cu", "cat")
Cu_DF[,c(1:3)][Cu_DF[,c(1:3)] <= 0] <- 0

sum(Cu_DF$Cu_SPM, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_COL, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))
sum(Cu_DF$Cu_TRD, na.rm=T)/(sum(Cu_DF$Cu_SPM,Cu_DF$Cu_COL,Cu_DF$Cu_TRD, na.rm=T))

names(Fe_DF)<-c("Fe_SPM","Fe_COL","Fe_TRD","TURNOVER", "Fe", "cat")
Fe_DF[,c(1:3)][Fe_DF[,c(1:3)] <= 0] <- 0

sum(Fe_DF$Fe_SPM, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_COL, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))
sum(Fe_DF$Fe_TRD, na.rm=T)/(sum(Fe_DF$Fe_SPM,Fe_DF$Fe_COL,Fe_DF$Fe_TRD, na.rm=T))

names(Mo_DF)<-c("Mo_SPM","Mo_COL","Mo_TRD","TURNOVER", "Mo", "cat")
Mo_DF[,c(1:3)][Mo_DF[,c(1:3)] <= 0] <- 0

sum(Mo_DF$Mo_SPM, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_COL, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))
sum(Mo_DF$Mo_TRD, na.rm=T)/(sum(Mo_DF$Mo_SPM,Mo_DF$Mo_COL,Mo_DF$Mo_TRD, na.rm=T))

names(Mn_DF)<-c("Mn_SPM","Mn_COL","Mn_TRD","TURNOVER", "Mn", "cat")
Mn_DF[,c(1:3)][Mn_DF[,c(1:3)] <= 0] <- 0

sum(Mn_DF$Mn_SPM, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_COL, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))
sum(Mn_DF$Mn_TRD, na.rm=T)/(sum(Mn_DF$Mn_SPM,Mn_DF$Mn_COL,Mn_DF$Mn_TRD, na.rm=T))

names(Ni_DF)<-c("Ni_SPM","Ni_COL","Ni_TRD","TURNOVER", "Ni", "cat")
Ni_DF[,c(1:3)][Ni_DF[,c(1:3)] <= 0] <- 0

sum(Ni_DF$Ni_SPM, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_COL, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))
sum(Ni_DF$Ni_TRD, na.rm=T)/(sum(Ni_DF$Ni_SPM,Ni_DF$Ni_COL,Ni_DF$Ni_TRD, na.rm=T))

names(P_DF)<-c("P_SPM","P_COL","P_TRD","TURNOVER", "P", "cat")
P_DF[,c(1:3)][P_DF[,c(1:3)] <= 0] <- 0

sum(P_DF$P_SPM, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_COL, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))
sum(P_DF$P_TRD, na.rm=T)/(sum(P_DF$P_SPM,P_DF$P_COL,P_DF$P_TRD, na.rm=T))

names(Pb_DF)<-c("Pb_SPM","Pb_COL","Pb_TRD","TURNOVER", "Pb", "cat")
Pb_DF[,c(1:3)][Pb_DF[,c(1:3)] <= 0] <- 0

sum(Pb_DF$Pb_SPM, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_COL, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))
sum(Pb_DF$Pb_TRD, na.rm=T)/(sum(Pb_DF$Pb_SPM,Pb_DF$Pb_COL,Pb_DF$Pb_TRD, na.rm=T))

names(S_DF)<-c("S_SPM","S_COL","S_TRD","TURNOVER", "S", "cat")
S_DF[,c(1:3)][S_DF[,c(1:3)] <= 0] <- 0

sum(S_DF$S_SPM, na.rm=T)/(sum(S_DF$S_SPM,S_DF$S_COL,S_DF$S_TRD, na.rm=T))
sum(S_DF$S_COL, na.rm=T)/(sum(S_DF$S_SPM,S_DF$S_COL,S_DF$S_TRD, na.rm=T))
sum(S_DF$S_TRD, na.rm=T)/(sum(S_DF$S_SPM,S_DF$S_COL,S_DF$S_TRD, na.rm=T))

names(Se_DF)<-c("Se_SPM","Se_COL","Se_TRD","TURNOVER", "Se", "cat")
Se_DF[,c(1:3)][Se_DF[,c(1:3)] <= 0] <- 0

sum(Se_DF$Se_SPM, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_COL, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))
sum(Se_DF$Se_TRD, na.rm=T)/(sum(Se_DF$Se_SPM,Se_DF$Se_COL,Se_DF$Se_TRD, na.rm=T))

names(Zn_DF)<-c("Zn_SPM","Zn_COL","Zn_TRD","TURNOVER", "Zn", "cat")
Zn_DF[,c(1:3)][Zn_DF[,c(1:3)] <= 0] <- 0

sum(Zn_DF$Zn_SPM, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_COL, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))
sum(Zn_DF$Zn_TRD, na.rm=T)/(sum(Zn_DF$Zn_SPM,Zn_DF$Zn_COL,Zn_DF$Zn_TRD, na.rm=T))

