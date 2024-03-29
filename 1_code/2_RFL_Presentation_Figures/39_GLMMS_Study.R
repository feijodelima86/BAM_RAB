library(readr)
library(dismo)
library(gbm)
library(effects)
library(AICcmodavg)
library(nlme)
library(MuMIn)

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_2.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

#### GLMs####

#COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

#dev.new()

### As ###



As_TD <- lme(As ~ UF_As_Oshift * SAMPLE_DESCRIPTOR, random=list (~1|SAMPLING_DATE, ~1|SITE), data = COMPARTMENTS_AVG, na.action=na.omit)

no.int_As_TD  <- lme(As ~ (UF_As_Oshift + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_As_TD <- lme(As ~ (UF_As_Oshift * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_As_TD2 <- lme(As ~ (UF_As_Oshift * TURNOVER) * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_As_TD)

summary(mix.int_As_TD)

summary(mix.int_As_TD2)

models <- list(As_TD, mix.int_As_TD, mix.int_As_TD2, no.int_As_TD)

aictab(c(models), modnames=c("As_TD","mix.int_As_TD","mix.int_As_TD2","no.int_As_TD"), method=ML)

plot(allEffects(no.int_As_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(no.int_As_TD)

anova(As_TD)

r.squaredGLMM(no.int_As_TD)

### Cd ###

Cd_TD <- lme(Cd ~ UF_Cd114 * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE, data = COMPARTMENTS_AVG, na.action=na.omit)

summary(Cd_TD)

anova(Cd_TD)

plot(allEffects(Cd_TD), lines = list(multiline = T), confint = list(style = "auto"))

no.int_Cd_TD  <- lme(Cd ~ (UF_Cd114 + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Cd_TD <- lme(Cd ~ (UF_Cd114 * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Cd_TD2 <- lme(Cd ~ UF_Cd114 * TURNOVER * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_Cd_TD)

summary(mix.int_Cd_TD)

summary(mix.int_Cd_TD2)

models <- list(Cd_TD, mix.int_Cd_TD, mix.int_Cd_TD2, no.int_Cd_TD)

aictab(c(models), modnames=c("Cd_TD","mix.int_Cd_TD","mix.int_Cd_TD2","no.int_Cd_TD"), method=ML)

plot(allEffects(Cd_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(Cd_TD)

anova(Cd_TD)

r.squaredGLMM(Cd_TD)

### Cu ###

Cu_TD <- lme(Cu ~ UF_Cu_NH3 * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE, data = COMPARTMENTS_AVG, na.action=na.omit)

summary(Cu_TD)

anova(Cu_TD)

plot(allEffects(Cu_TD), lines = list(multiline = T), confint = list(style = "auto"))

no.int_Cu_TD  <- lme(Cu ~ (UF_Cu_NH3 + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Cu_TD <- lme(Cu ~ (UF_Cu_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Cu_TD2 <- lme(Cu ~ UF_Cu_NH3 * TURNOVER * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_Cu_TD)

summary(mix.int_Cu_TD)

summary(mix.int_Cu_TD2)

models <- list(Cu_TD, mix.int_Cu_TD, mix.int_Cu_TD2, no.int_Cu_TD)

aictab(c(models), modnames=c("Cu_TD","mix.int_Cu_TD","mix.int_Cu_TD2","no.int_Cu_TD"), method=ML)

plot(allEffects(no.int_Cu_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(Cu_TD)

anova(Cu_TD)

r.squaredGLMM(Cu_TD)

### Fe ###

Fe_TD <- lme(Fer.1 ~ UF_Fe_NH3 * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE, data = COMPARTMENTS_AVG, na.action=na.omit)

summary(Fe_TD)

anova(Fe_TD)

plot(allEffects(Fe_TD), lines = list(multiline = T), confint = list(style = "auto"))

no.int_Fe_TD  <- lme(Fer.1 ~ (UF_Fe_NH3 + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Fe_TD <- lme(Fer.1 ~ (UF_Fe_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Fe_TD2 <- lme(Fer.1 ~ UF_Fe_NH3 * TURNOVER * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_Fe_TD)

summary(mix.int_Fe_TD)

summary(mix.int_Fe_TD2)

models <- list(Fe_TD, mix.int_Fe_TD, mix.int_Fe_TD2, no.int_Fe_TD)

aictab(c(models), modnames=c("Fe_TD","mix.int_Fe_TD","mix.int_Fe_TD2","no.int_Fe_TD"), method=ML)

plot(allEffects(mix.int_Fe_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(mix.int_Fe_TD2)

anova(mix.int_Fe_TD2)

r.squaredGLMM(mix.int_Fe_TD2)

### Mo ###

Mo_TD <- lme(Mo ~ UF_Mo_Oshift * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE, data = COMPARTMENTS_AVG, na.action=na.omit)

summary(Mo_TD)

anova(Mo_TD)

plot(allEffects(Mo_TD), lines = list(multiline = T), confint = list(style = "auto"))

no.int_Mo_TD  <- lme(Mo ~ (UF_Mo_Oshift + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Mo_TD <- lme(Mo ~ (UF_Mo_Oshift * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_Mo_TD2 <- lme(Mo ~ UF_Mo_Oshift * TURNOVER * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_Mo_TD)

summary(mix.int_Mo_TD)

summary(mix.int_Mo_TD2)

models <- list(Mo_TD, mix.int_Mo_TD, mix.int_Mo_TD2, no.int_Mo_TD)

aictab(c(models), modnames=c("Mo_TD","mix.int_Mo_TD","mix.int_Mo_TD2","no.int_Mo_TD"), method=ML)

plot(allEffects(Mo_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(Mo_TD)

anova(Mo_TD)

r.squaredGLMM(Mo_TD)

### P ###

P_TD <- lme(P ~ UF_P * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE, data = COMPARTMENTS_AVG, na.action=na.omit)

summary(P_TD)

anova(P_TD)

plot(allEffects(P_TD), lines = list(multiline = T), confint = list(style = "auto"))

no.int_P_TD  <- lme(P ~ (UF_P + TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_P_TD <- lme(P ~ (UF_P * TURNOVER) : SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

mix.int_P_TD2 <- lme(P ~ UF_P * TURNOVER * SAMPLE_DESCRIPTOR, random=~1|SAMPLING_DATE/SITE , data = COMPARTMENTS_AVG, na.action=na.omit)

summary(no.int_P_TD)

summary(mix.int_P_TD)

summary(mix.int_P_TD2)

models <- list(P_TD, mix.int_P_TD, mix.int_P_TD2, no.int_P_TD)

aictab(c(models), modnames=c("P_TD","mix.int_P_TD","mix.int_P_TD2","no.int_P_TD"), method=ML)

plot(allEffects(no.int_P_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(no.int_P_TD)

anova(no.int_P_TD)

r.squaredGLMM(no.int_P_TD)

### Pb ###

Pb_TD <- glm(Pb ~ UF_Pb_NH3 : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(Pb_TD)

with(summary(Pb_TD), 1 - deviance/null.deviance)

no.int_Pb_TD    <- glm(Pb ~ (UF_Pb_NH3 + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Pb_TD   <- glm(Pb ~ (UF_Pb_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Pb_TD2  <- glm(Pb ~ (UF_Pb_NH3 * TURNOVER) * SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(no.int_Pb_TD)

summary(mix.int_Pb_TD)

summary(mix.int_Pb_TD2)

models <- list(Pb_TD, mix.int_Pb_TD, mix.int_Pb_TD2, no.int_Pb_TD)

aictab(c(models), modnames=c("Pb_TD","mix.int_Pb_TD","mix.int_Pb_TD2","no.int_Pb_TD"))

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(mix.int_Pb_TD)

r.squaredGLMM(mix.int_Pb_TD)

### Se ###

Se_TD <- glm(Se ~ UF_Se82_Oshift : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(Se_TD)

with(summary(Pb_TD), 1 - deviance/null.deviance)

no.int_Se_TD    <- glm(Se ~ (UF_Se82_Oshift + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Se_TD   <- glm(Se ~ (UF_Se82_Oshift * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Se_TD2  <- glm(Se ~ (UF_Se82_Oshift * TURNOVER) * SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(no.int_Se_TD)

summary(mix.int_Se_TD)

summary(mix.int_Se_TD2)

models <- list(Se_TD, mix.int_Se_TD, mix.int_Se_TD2, no.int_Se_TD)

aictab(c(models), modnames=c("Se_TD","mix.int_Se_TD","mix.int_Se_TD2","no.int_Se_TD"))

with(summary(no.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(no.int_Se_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(no.int_Se_TD)

r.squaredGLMM(no.int_Se_TD)

### Zn ###

Zn_TD <- glm(Zn ~ UF_Zn : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(Zn_TD)

with(summary(Zn_TD), 1 - deviance/null.deviance)

no.int_Zn_TD    <- glm(Zn ~ (UF_Zn + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Zn_TD   <- glm(Zn ~ (UF_Zn * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

mix.int_Zn_TD2  <- glm(Zn ~ (UF_Zn * TURNOVER) * SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(no.int_Zn_TD)

summary(mix.int_Zn_TD)

summary(mix.int_Zn_TD2)

models <- list(Zn_TD, mix.int_Zn_TD, mix.int_Zn_TD2, no.int_Zn_TD)

aictab(c(models), modnames=c("Zn_TD","mix.int_Zn_TD","mix.int_Zn_TD2","no.int_Zn_TD"))

with(summary(no.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(no.int_Zn_TD), lines = list(multiline = T), confint = list(style = "auto"))

summary(no.int_Zn_TD)

r.squaredGLMM(no.int_Zn_TD)

### As: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,7,93,116,13)]),[,c(3,4,7,93,116,13)]))
### Cd: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,7,98,116,19)]),c(3,4,7,93,116,19)]))
### Cu: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,112,22)]),c(3,4,6,7,10,112,22)]))
### Fe: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,111,25)]),c(3,4,6,7,10,111,25)]))
### Mo: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,116,31)]),c(3,4,6,7,10,116,31)]))
### Pb: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,118,35)]),c(3,4,6,7,10,118,35)]))
### P:  br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,109,34)]),c(3,4,6,7,10,109,34)]))
### Zn: br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(4,7,10,7,10,113,49)]),c(3,4,6,7,10,113,49)]))

br.sub<-data.frame((COMPARTMENTS_AVG[complete.cases(COMPARTMENTS_AVG[,c(3,4,7,96,116,19)]),c(3,4,7,96,116,19)]))
br.sub$SAMPLE_DESCRIPTOR<-as.factor(br.sub$SAMPLE_DESCRIPTOR)
br.sub$SITE<-as.factor(br.sub$SITE)

#br.sub[,6]<-log10(br.sub[,6]+1)
#br.sub[,7]<-log10(br.sub[,7]+1)

#br.sub$TOTAL.MG.M2.RITCHIE<-log10(br.sub$TOTAL.MG.M2.RITCHIE+1)
#br.sub$OM.AREA.g.m2<-log10(br.sub$OM.AREA.g.m2+1)


br.sub[,ncol(br.sub)]<-br.sub[,ncol(br.sub)]*1000

UCFR.CD.tc6.lr002 <- gbm.step(data=br.sub, 
                              gbm.x = c(2,3,4,5),
                              gbm.y = ncol(br.sub),
                              family = "gaussian",
                              tree.complexity = 2,
                              learning.rate = 0.007,
                              bag.fraction = 0.5
)


predfit<-lm(UCFR.CD.tc6.lr002$fitted ~ br.sub[,6])

summary(predfit)

plot(br.sub[,6], UCFR.CD.tc6.lr002$fitted, 
     #xlim=c(0.5,3), 
     #ylim=c(0.5,3), 
     col=br.sub$SAMPLE_DESCRIPTOR, 
     xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

#dev.new()

gbm.plot(UCFR.CD.tc6.lr002, write.title = F, nplots = 5, plot.layout= c(2,2),
         las=1,
         lwd=2,
         cex.lab=1.5,
         cex.axis=1,
         smooth =F,
         rug=F,
         y.label = NA
)

gbm.plot.fits(UCFR.CD.tc6.lr002)

find.int.cd <- gbm.interactions(UCFR.CD.tc6.lr002)

find.int.cd$interactions

find.int.cd$rank.list

gbm.perspec(UCFR.CD.tc6.lr002,4,3, 
            z.range=c(0,5), 
            theta = 315,
            phi=45, 
            cex.lab = 1.2, font.lab = 2, cex.axis = 1, font.axis= 1,
            perspective = T,
            smooth=F)

##### GLM #####

require(scatterplot3d)
library(effects)
library(pscl)

names(br.sub)

tdsp<-scatterplot3d(br.sub[,4],br.sub[,5],br.sub[,6], pch = 19, color = factor(br.sub$SAMPLE_DESCRIPTOR, labels = c("darkgreen","gold","chartreuse")), main="3D Scatterplot",
                    type="h", angle = 120)

fitm <- lm(br.sub[,6] ~ br.sub[,4]+br.sub[,5]) 

summary(fitm)

tdsp$plane3d(fitm)

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

mix.int_As_TD <- glm(As ~ (TD.As + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_As_TD))

mix.int_Ca_TD <- glm(Car ~ (TD.Ca + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Ca_TD)

with(summary(mix.int_Ca_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Ca_TD))

mix.int_Cd_TD <- glm(Cd ~ (TD.Cd + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cd_TD))

mix.int_Cu_TD <- glm(Cu ~ (TD.Cu + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cu_TD))

mix.int_Fe_TD <- glm(Fe ~ (TD.Fe + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Fe_TD))

mix.int_Mo_TD <- glm(Mo ~ (TD.Mo + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Mo_TD)

plot(allEffects(mix.int_Mo_TD))

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

mix.int_P_TD <- glm(P ~ (TD.P + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_P_TD))

mix.int_Pb_TD <- glm(Pb ~ (TD.Pb + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD))

mix.int_Se_TD <- glm(Se ~ (TD.Se + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Se_TD))

mix.int_Zn_TD <- glm(Zn ~ (TD.Zn + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Zn_TD))

names(COMPARTMENTS_AVG)


names(COMPARTMENTS_AVG)

#dev.new()

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"),]

mix.int_As_TD <- glm(As ~ (UF_As + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_As_TD))

mix.int_Cd_TD <- glm(Cd ~ (UF_Cd + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cd_TD))

mix.int_Cu_TD <- glm(Cu ~ (UF_Cu + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cu_TD))

mix.int_Fe_TD <- glm(Fe ~ (UF_Fe + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Fe_TD))

mix.int_Mo_TD <- glm(Mo ~ (UF_Mo + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Mo_TD)

plot(allEffects(mix.int_Mo_TD))

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

mix.int_P_TD <- glm(P ~ (UF_P + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_P_TD))

mix.int_Pb_TD <- glm(Pb ~ (UF_Pb + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD))

mix.int_Se_TD <- glm(Se ~ (UF_Se  + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Se_TD))

mix.int_Zn_TD <- glm(Zn ~ (UF_Zn   + TURNOVER + OM.AREA.g.m2), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Zn_TD))

####

names(COMPARTMENTS_AVG)

#dev.new()

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"),]

mix.int_As_TD <- glm(As ~ (UF_As * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_As_TD))

mix.int_Cd_TD <- glm(Cd ~ (UF_Cd * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cd_TD))

mix.int_Cu_TD <- glm(Cu ~ (UF_Cu * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cu_TD))

mix.int_Fe_TD <- glm(Fe ~ (UF_Fe * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Fe_TD))

mix.int_Mo_TD <- glm(Mo ~ (UF_Mo * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Mo_TD)

plot(allEffects(mix.int_Mo_TD))

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

mix.int_P_TD <- glm(P ~ (UF_P * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_P_TD))

mix.int_Pb_TD <- glm(Pb ~ (UF_Pb * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD))

mix.int_Se_TD <- glm(Se ~ (UF_Se  * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Se_TD))

mix.int_Zn_TD <- glm(Zn ~ (UF_Zn  * TURNOVER), data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Zn_TD))


names(COMPARTMENTS_AVG)

