library(gbm)
library(dismo)
library(MASS)
library(usdm)
library(dplyr)
library(cowplot)
library("car")

## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG_2.csv"))

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

names(COMPARTMENTS_AVG)

Ni_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","W_Ni","F_Ni","UF_Ni","TURNOVER","Ni","TOTAL.MG.M2.RITCHIE")])

Ni_TRD <- (Ni_DF[,c(1,3,5,7,6)])

names(Ni_TRD)<- c("SAMPLE_DESCRIPTOR","Ni_TRD", "TURNOVER", "OM", "Ni")

Ni_TRD<-Ni_TRD[Ni_TRD$Ni_TRD>0,]

Ni_TRD<-Ni_TRD[Ni_TRD$OM>0,]

#Ni_TRD$Ni <- log10(Ni_TRD$Ni+1)

Ni_TRD$TURNOVER <- log10(Ni_TRD$TURNOVER+1)

Ni_TRD$Ni_TRD <- log10(Ni_TRD$Ni_TRD+1)

Ni_TRD$OM <- log10(Ni_TRD$OM+1)


for (col in colnames(Ni_TRD[,c(1)])) {
  Ni_TRD[[col]] <- replace_outliers_with_na(Ni_TRD[[col]])
}

Ni_TRD <- Ni_TRD[complete.cases(Ni_TRD), ]


UCFR.SS.tc5.lr002 <- gbm.step(data=Ni_TRD, 
                              gbm.x = c(2,3,4),
                              gbm.y = 5,
                              family = "gaussian",
                              tree.complexity = 2,
                              learning.rate = 0.002,
                              bag.fraction = 0.25
)


### Nicer plots

y.bar <- min(UCFR.SS.tc5.lr002$cv.values) 
y.min <- min(UCFR.SS.tc5.lr002$cv.values - UCFR.SS.tc5.lr002$cv.loss.ses)
y.max <- max(UCFR.SS.tc5.lr002$cv.values + UCFR.SS.tc5.lr002$cv.loss.ses) 

par(mar=c(5,5,3,1)+.1, mgp=c(3,1,0))
plot(UCFR.SS.tc5.lr002$trees.fitted, UCFR.SS.tc5.lr002$cv.values, 
     type = 'l', 
     ylab = "Holdout deviance",
     xlab = "Number of trees", 
     ylim = c(y.min,y.max), 
     lwd=1.5,        
     las=1,
     font.lab= 2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

abline(h = y.bar, col = 2)

lines(UCFR.SS.tc5.lr002$trees.fitted, UCFR.SS.tc5.lr002$cv.values + UCFR.SS.tc5.lr002$cv.loss.ses, lty=2)  
lines(UCFR.SS.tc5.lr002$trees.fitted, UCFR.SS.tc5.lr002$cv.values - UCFR.SS.tc5.lr002$cv.loss.ses, lty=2)  

target.trees <- UCFR.SS.tc5.lr002$trees.fitted[match(TRUE,UCFR.SS.tc5.lr002$cv.values == y.bar)]
abline(v = target.trees, col=3, lwd=2)
box(lwd=2)


# Predicred vs observed linear regression for standing stocks. 

CHL_DF<-data.frame(UCFR.SS.tc5.lr002$fitted, Ni_TRD$Ni)

predfit<-lm(UCFR.SS.tc5.lr002$fitted ~ Ni_TRD$Ni)

summary(predfit)

par(mfrow=c(1,1))

plot(Ni_TRD$Ni, UCFR.SS.tc5.lr002$fitted, xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

plot_4a<-gbm.plot(UCFR.SS.tc5.lr002, write.title = F, nplots = 3, plot.layout= c(2,2),
                  las=1,
                  lwd=2,
                  cex.lab=2,
                  cex.axis=1.5,
                  smooth =F,
                  rug=F,
                  y.label = NA
)

gbm.plot.fits(UCFR.SS.tc5.lr002)

names(Ni_TRD)

mix.int_Ni <- glm(Ni ~ (Ni_TRD+TURNOVER+OM), data = Ni_TRD, family=gaussian)

Anova(mix.int_Ni)

summary(mix.int_Ni)

with(summary(mix.int_Ni), 1 - deviance/null.deviance)

plot(allEffects(mix.int_Ni))

plot(Ni_TRD$OM, Ni_TRD$Ni)

