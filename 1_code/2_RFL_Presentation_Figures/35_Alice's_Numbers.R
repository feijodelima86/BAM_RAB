### 3: Add Alice's true numbers to assess metals turnover for different compartments. 

library(readr)
library(gbm)
library(dismo)
library(MASS)
library(usdm)
library(dplyr)
library(cowplot)

X2021_turnovers_for_Rafa_AC <- read_csv("0_data/external/2021_turnovers_for_Rafa_AC.csv")
COMPARTMENTS <- read_csv("2_incremental/COMPARTMENTS_AVG_7.csv")
colnames(COMPARTMENTS)[c(2:4)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

names(COMPARTMENTS)

EPIL <- COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]
EPIP <- COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]
FILA <- COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]

slm.EPIL<-lm(Pb ~ TD.Pb + OM.AREA.g.m2 + TURNOVER, data = EPIL)
slm.EPIP<-lm(Pb ~ TD.Pb + OM.AREA.g.m2 + TURNOVER, data = EPIP)
slm.FILA<-lm(Pb ~ TD.Pb + OM.AREA.g.m2 + TURNOVER, data = FILA)
slm.ALLS<-lm(Cd ~ TD.Cd + OM.AREA.g.m2 + TURNOVER, data = COMPARTMENTS)


summary(slm.EPIL)
summary(slm.EPIP)
summary(slm.FILA)
summary(slm.ALLS)



mix.int_Pb_TD <- glm(Pb ~ TD.Pb + OM.AREA.g.m2 + TURNOVER + SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian(link="log"))

summary(mix.int_Cd_TD)

#### BRT ####

#56, 16 Cd 17, 57

cd.sub<-data.frame((COMPARTMENTS[complete.cases(COMPARTMENTS[,c(3,4,5,10,13,36,16)]),c(3,4,5,10,13,36,16)]))
cd.sub$SAMPLE_DESCRIPTOR<-as.factor(cd.sub$SAMPLE_DESCRIPTOR)
cd.sub$SITE<-as.factor(cd.sub$SITE)
cd.sub[,6]<-log10(cd.sub[,6]+1)
cd.sub[,7]<-log10(cd.sub[,7]+1)
cd.sub$TOTAL.MG.M2.RITCHIE<-log10(cd.sub$TOTAL.MG.M2.RITCHIE+1)
cd.sub$OM.AREA.g.m2<-log10(cd.sub$OM.AREA.g.m2+1)

names(COMPARTMENTS)

UCFR.CD.tc6.lr002 <- gbm.step(data=cd.sub, 
                              gbm.x = c(2,3,4,6),
                              gbm.y = ncol(cd.sub),
                              family = "gaussian",
                              tree.complexity = 4,
                              learning.rate = 0.007,
                              bag.fraction = 0.5
)


predfit<-lm(UCFR.CD.tc6.lr002$fitted ~ cd.sub[,7])

summary(predfit)

par(mfrow=c(1,2))

plot(cd.sub[,7], UCFR.CD.tc6.lr002$fitted, 
     #xlim=c(0.5,3), 
     #ylim=c(0.5,3), 
     col=cd.sub$SAMPLE_DESCRIPTOR, 
     xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

dev.new()

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

gbm.perspec(UCFR.CD.tc6.lr002,3,2, 
            z.range=c(0.0005,0.002), 
            theta = 225,
            phi=45, 
            cex.lab = 1.2, font.lab = 2, cex.axis = 1, font.axis= 1,
            perspective = T,
            smooth=F)

##### GLM #####



