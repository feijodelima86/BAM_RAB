library(readr)
library(effects)
library(AICcmodavg)
library(nlme)
library(car)
library(interactions)
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))
library(lme4)


COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

#### GLMs####

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

#dev.new()

### As ###

As_DF <- (COMPARTMENTS_AVG[,c("F_As_Oshift","TURNOVER","As")])

As_DF <- (COMPARTMENTS_AVG[,c("F_As_Oshift","TURNOVER","As")])

As_DF <- As_DF[complete.cases(As_DF), ]

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

As_TD <- glm(As ~ F_As_Oshift, data = As_DF, family=gaussian)

no.int_As_TD  <- glm(As ~ (F_As_Oshift + TURNOVER), data = As_DF, family=gaussian)

mix.int_As_TD <- glm(As ~ (F_As_Oshift * TURNOVER), data = As_DF, family=gaussian)

models <- list(As_TD, mix.int_As_TD, no.int_As_TD)

aictab(c(models), modnames=c("As_TD","mix.int_As_TD","no.int_As_TD"), method=ML)

plot(allEffects(mix.int_As_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_As_TD))

Anova(mix.int_As_TD)

### PLOT ###

p=ggplot(data=As_DF,aes(F_As_Oshift,As))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_As_Oshift, y=As, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "As  (TD, mg/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p

### Cd ###

Cd_DF <- (COMPARTMENTS_AVG[,c("UF_Cd114","TURNOVER","Cd")])

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

Cd_DF <- Cd_DF[complete.cases(Cd_DF), ]

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Cd_TD <- glm(Cd ~ UF_Cd114, data = Cd_DF, family=gaussian)

no.int_Cd_TD  <- glm(Cd ~ (UF_Cd114 + TURNOVER), data = Cd_DF, family=gaussian)

mix.int_Cd_TD <- glm(Cd ~ (UF_Cd114 * TURNOVER), data = Cd_DF, family=gaussian)

models <- list(Cd_TD, mix.int_Cd_TD, no.int_Cd_TD)

aictab(c(models), modnames=c("Cd_TD","mix.int_Cd_TD","no.int_Cd_TD"), method=ML)

plot(allEffects(mix.int_Cd_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Cd_TD))

Anova(mix.int_Cd_TD)

### PLOT ###

p=ggplot(data=Cd_DF,aes(UF_Cd114,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=UF_Cd114, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (TD, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p


#### Cu - Removed one outlier - v high -20


Cu_DF <- (COMPARTMENTS_AVG[,c("UF_Cu_NH3","TURNOVER","Cu")])

Cu_DF <- Cu_DF[complete.cases(Cu_DF), ]

Cu_DF$TURNOVER <- asin(sqrt(Cu_DF$TURNOVER ))

quantile(Cu_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cu_DF$cat <- cut(Cu_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Cu_TD <- glm(Cu ~ UF_Cu_NH3, data = Cu_DF, family=gaussian)

no.int_Cu_TD  <- glm(Cu ~ (UF_Cu_NH3 + TURNOVER), data = Cu_DF, family=gaussian)

mix.int_Cu_TD <- glm(Cu ~ (UF_Cu_NH3 * TURNOVER), data = Cu_DF, family=gaussian)

models <- list(Cu_TD, mix.int_Cu_TD, no.int_Cu_TD)

aictab(c(models), modnames=c("Cu_TD","mix.int_Cu_TD","no.int_Cu_TD"), method=ML)

plot(allEffects(mix.int_Cu_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Cu_TD))

Anova(mix.int_Cu_TD)

### PLOT ###

p=ggplot(data=Cu_DF,aes(UF_Cu_NH3,Cu))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=UF_Cu_NH3, y=Cu, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cu  (TD, mg/l)",    
       y= "Cu Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p


#### Fe - Removed 2 outliers - v high -28,5

Fe_DF <- (COMPARTMENTS_AVG[,c("F_Fe_NH3","TURNOVER","Fe")])

Fe_DF <- Fe_DF[complete.cases(Fe_DF), ]

Fe_DF$TURNOVER <- asin(sqrt(Fe_DF$TURNOVER ))

quantile(Fe_DF$TURNOVER, probs = seq(0, 1, 1/5))

Fe_DF$cat <- cut(Fe_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Fe_TD <- glm(Fe ~ F_Fe_NH3, data = Fe_DF, family=gaussian)

no.int_Fe_TD  <- glm(Fe ~ (F_Fe_NH3 + TURNOVER), data = Fe_DF, family=gaussian)

mix.int_Fe_TD <- glm(Fe ~ (F_Fe_NH3 * TURNOVER), data = Fe_DF, family=gaussian)

models <- list(Fe_TD, mix.int_Fe_TD, no.int_Fe_TD)

aictab(c(models), modnames=c("Fe_TD","mix.int_Fe_TD","no.int_Fe_TD"), method=ML)

plot(allEffects(mix.int_Fe_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Fe_TD))

Anova(mix.int_Fe_TD)

### PLOT ###

p=ggplot(data=Fe_DF,aes(F_Fe_NH3,Fe))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Fe_NH3, y=Fe, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Fe  (TD, mg/l)",    
       y= "Fe Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p

#### Mo - Removed 1 outlier, altough no leverage

Mo_DF <- (COMPARTMENTS_AVG[,c("F_Mo_Oshift","TURNOVER","Mo")])

Mo_DF <- Mo_DF[complete.cases(Mo_DF), ]

Mo_DF$TURNOVER <- asin(sqrt(Mo_DF$TURNOVER ))

quantile(Mo_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mo_DF$cat <- cut(Mo_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Mo_TD <- glm(Mo ~ F_Mo_Oshift, data = Mo_DF, family=gaussian)

no.int_Mo_TD  <- glm(Mo ~ (F_Mo_Oshift + TURNOVER), data = Mo_DF, family=gaussian)

mix.int_Mo_TD <- glm(Mo ~ (F_Mo_Oshift * TURNOVER), data = Mo_DF, family=gaussian)

models <- list(Mo_TD, mix.int_Mo_TD, no.int_Mo_TD)

aictab(c(models), modnames=c("Mo_TD","mix.int_Mo_TD","no.int_Mo_TD"), method=ML)

plot(allEffects(mix.int_Mo_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Mo_TD))

Anova(mix.int_Mo_TD)

### PLOT ###

p=ggplot(data=Mo_DF,aes(F_Mo_Oshift,Mo))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Mo_Oshift, y=Mo, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mo  (TD, mg/l)",    
       y= "Mo Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p

### P ###

P_DF <- (COMPARTMENTS_AVG[,c("F_P","TURNOVER","P")])

P_DF <- P_DF[complete.cases(P_DF), ]

P_DF$TURNOVER <- asin(sqrt(P_DF$TURNOVER ))

quantile(P_DF$TURNOVER, probs = seq(0, 1, 1/5))

P_DF$cat <- cut(P_DF$TURNOVER, 
                breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

P_TD <- glm(P ~ F_P, data = P_DF, family=gaussian)

no.int_P_TD  <- glm(P ~ (F_P + TURNOVER), data = P_DF, family=gaussian)

mix.int_P_TD <- glm(P ~ (F_P * TURNOVER), data = P_DF, family=gaussian)

models <- list(P_TD, mix.int_P_TD, no.int_P_TD)

aictab(c(models), modnames=c("P_TD","mix.int_P_TD","no.int_P_TD"), method=ML)

plot(allEffects(mix.int_P_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_P_TD))

Anova(mix.int_P_TD)

### PLOT ###

p=ggplot(data=P_DF,aes(F_P,P))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_P, y=P, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "P  (TD, mg/l)",    
       y= "P Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p


### Pb ###

Pb_DF <- (COMPARTMENTS_AVG[,c("UF_Pb_NH3","TURNOVER","Pb")])

Pb_DF <- Pb_DF[complete.cases(Pb_DF), ]

Pb_DF$TURNOVER <- asin(sqrt(Pb_DF$TURNOVER ))

quantile(Pb_DF$TURNOVER, probs = seq(0, 1, 1/5))

Pb_DF$cat <- cut(Pb_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Pb_TD <- glm(Pb ~ UF_Pb_NH3, data = Pb_DF, family=gaussian)

no.int_Pb_TD  <- glm(Pb ~ (UF_Pb_NH3 + TURNOVER), data = Pb_DF, family=gaussian)

mix.int_Pb_TD <- glm(Pb ~ (UF_Pb_NH3 * TURNOVER), data = Pb_DF, family=gaussian)

models <- list(Pb_TD, mix.int_Pb_TD, no.int_Pb_TD)

aictab(c(models), modnames=c("Pb_TD","mix.int_Pb_TD","no.int_Pb_TD"), method=ML)

plot(allEffects(mix.int_Pb_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Pb_TD))

Anova(mix.int_Pb_TD)

### PLOT ###

p=ggplot(data=Pb_DF,aes(UF_Pb_NH3,Pb))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=UF_Pb_NH3, y=Pb, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Pb  (TD, mg/l)",    
       y= "Pb Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p

pblm<-lm(Pb_DF$Pb ~ Pb_DF$UF_Pb_NH3)

summary(pblm)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

### Se ###

Se_DF <- (COMPARTMENTS_AVG[,c("F_Se82_Oshift","TURNOVER","Se")])

Se_DF <- Se_DF[complete.cases(Se_DF), ]

Se_DF$TURNOVER <- asin(sqrt(Se_DF$TURNOVER ))

quantile(Se_DF$TURNOVER, probs = seq(0, 1, 1/5))

Se_DF$cat <- cut(Se_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Se_TD <- glm(Se ~ F_Se82_Oshift, data = Se_DF, family=gaussian)

no.int_Se_TD  <- glm(Se ~ (F_Se82_Oshift + TURNOVER), data = Se_DF, family=gaussian)

mix.int_Se_TD <- glm(Se ~ (F_Se82_Oshift * TURNOVER), data = Se_DF, family=gaussian)

models <- list(Se_TD, mix.int_Se_TD, no.int_Se_TD)

aictab(c(models), modnames=c("Se_TD","mix.int_Se_TD","no.int_Se_TD"), method=ML)

plot(allEffects(mix.int_Se_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Se_TD))

Anova(mix.int_Se_TD)

### PLOT ###

p=ggplot(data=Se_DF,aes(F_Se82_Oshift,Se))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Se82_Oshift, y=Se, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Se  (TD, mg/l)",    
       y= "Se Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p

### Zn ###

Zn_DF <- (COMPARTMENTS_AVG[,c("F_Zn","TURNOVER","Zn")])

Zn_DF <- Zn_DF[complete.cases(Zn_DF), ]

Zn_DF$TURNOVER <- asin(sqrt(Zn_DF$TURNOVER ))

quantile(Zn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Zn_DF$cat <- cut(Zn_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2368515 , 0.2519822      , 0.2795858     , 0.3157512 , Inf), 
                 labels=c("0.17 - 0.23", "0.23 - 0.25","0.25 - 0.27","0.27 - 0.31", "0.31 - 0.38"))

Zn_TD <- glm(Zn ~ F_Zn, data = Zn_DF, family=gaussian)

no.int_Zn_TD  <- glm(Zn ~ (F_Zn + TURNOVER), data = Zn_DF, family=gaussian)

mix.int_Zn_TD <- glm(Zn ~ (F_Zn * TURNOVER), data = Zn_DF, family=gaussian)

models <- list(Zn_TD, mix.int_Zn_TD, no.int_Zn_TD)

aictab(c(models), modnames=c("Zn_TD","mix.int_Zn_TD","no.int_Zn_TD"), method=ML)

plot(allEffects(mix.int_Zn_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Zn_TD))

Anova(mix.int_Zn_TD)

### PLOT ###

p=ggplot(data=Zn_DF,aes(F_Zn,Zn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Zn, y=Zn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Zn  (TD, mg/l)",    
       y= "Zn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p



