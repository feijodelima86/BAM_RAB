library(gbm)
library(dismo)
library(MASS)
library(usdm)
library(dplyr)
library(cowplot)

## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG.csv"))

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

names(COMPARTMENTS_AVG)

S_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_S","colloidal_S","trulydissolved_S","TURNOVER","S","OM.AREA.g.m2")])

S_TRD <- (S_DF[,c(1,4,5,7,6)])

names(S_TRD)<- c("SAMPLE_DESCRIPTOR","S_TRD", "TURNOVER", "OM", "S")

S_TRD<-S_TRD[S_TRD$S_TRD>0,]

S_TRD<-S_TRD[S_TRD$OM>0,]

#S_TRD$S_TRD <- log10(S_TRD$S_TRD+1)

S_TRD$S <- log10(S_TRD$S+1)

#S_TRD$OM <- log10(S_TRD$OM+1)


for (col in colnames(S_TRD[,c(1)])) {
  S_TRD[[col]] <- replace_outliers_with_na(S_TRD[[col]])
}

S_TRD <- S_TRD[complete.cases(S_TRD), ]


UCFR.SS.tc5.lr002 <- gbm.step(data=S_TRD, 
                              gbm.x = c(2,3,4),
                              gbm.y = 5,
                              family = "gaussian",
                              tree.complexity = 3,
                              learning.rate = 0.002,
                              bag.fraction = 0.75
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

CHL_DF<-data.frame(UCFR.SS.tc5.lr002$fitted, S_TRD$S)

predfit<-lm(UCFR.SS.tc5.lr002$fitted ~ S_TRD$S)

summary(predfit)

par(mfrow=c(1,1))

plot(S_TRD$S, UCFR.SS.tc5.lr002$fitted, xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

plot_4a<-gbm.plot(UCFR.SS.tc5.lr002, write.title = F, nplots = 3, plot.layout= c(1,3),
                  las=1,
                  lwd=2,
                  cex.lab=2,
                  cex.axis=1.5,
                  smooth =F,
                  rug=F,
                  y.label = NA
)

gbm.plot.fits(UCFR.SS.tc5.lr002)

