library(readr)
library(dismo)
library(gbm)
library(effects)

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#### GLMs####

#COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)


#dev.new()

mix.int_As_TD <- glm(As ~ (F_As_Oshift * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_As_TD))

mix.int_Cd_TD <- glm(Cd ~ (F_Cd114 * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cd_TD))

mix.int_Cu_TD <- glm(Cu ~ (F_Cu_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Cu_TD))

mix.int_Fe_TD <- glm(Fe ~ (F_Fe_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Fe_TD))

mix.int_Mo_TD <- glm(Mo ~ (F_Mo_Oshift + TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Mo_TD)

plot(allEffects(mix.int_Mo_TD))

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

mix.int_P_TD <- glm(P ~ (F_P * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_P_TD))

mix.int_Pb_TD <- glm(Pb ~ (F_Pb_NH3 * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Pb_TD))

mix.int_Se_TD <- glm(Se ~ (F_Se82_Oshift * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Se_TD))

mix.int_Zn_TD <- glm(Zn ~ (F_Zn * TURNOVER) : SAMPLE_DESCRIPTOR, data = COMPARTMENTS_AVG, family=gaussian)

summary(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Zn_TD))




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

