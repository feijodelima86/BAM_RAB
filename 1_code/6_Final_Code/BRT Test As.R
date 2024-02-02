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

As_DF <- (COMPARTMENTS_AVG[,c("SAMPLE_DESCRIPTOR","spm_As_Oshift","colloidal_As_Oshift","trulydissolved_As_Oshift","TURNOVER","As","OM.AREA.g.m2")])

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

As_TRD <- (As_DF[,c(1,4,5,7,6)])

names(As_TRD)<- c("SAMPLE_DESCRIPTOR","As_TRD", "TURNOVER", "OM", "As")

As_TRD<-As_TRD[As_TRD$As_TRD>0,]

As_TRD<-As_TRD[As_TRD$OM>0,]

#As_TRD$As_TRD <- log10(As_TRD$As_TRD+1)

#As_TRD$As <- log10(As_TRD$As+1)

As_TRD$OM <- log10(As_TRD$OM+1)


for (col in colnames(As_TRD[,c(1)])) {
  As_TRD[[col]] <- replace_outliers_with_na(As_TRD[[col]])
}

As_TRD <- As_TRD[complete.cases(As_TRD), ]

UCFR.SS.tc5.lr002 <- gbm.step(data=As_TRD, 
                              gbm.x = c(2,3,4),
                              gbm.y = 5,
                              family = "gaussian",
                              tree.complexity = 2,
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

CHL_DF<-data.frame(UCFR.SS.tc5.lr002$fitted, As_TRD$As)

predfit<-lm(UCFR.SS.tc5.lr002$fitted ~ As_TRD$As)

summary(predfit)

par(mfrow=c(1,1))

plot(As_TRD$As, UCFR.SS.tc5.lr002$fitted, xlab="Obs", ylab="Fitted Values", 
     abline(a=0, b=1),
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)

box(lwd=2)

plot_4a<-gbm.plot(UCFR.SS.tc5.lr002, write.title = F, nplots = 3, plot.layout= c(2,3),
                  las=1,
                  lwd=2,
                  cex.lab=2,
                  cex.axis=1.5,
                  smooth =F,
                  rug=F,
                  y.label = NA
)

gbm.plot.fits(UCFR.SS.tc5.lr002)


find.int <- gbm.interactions(UCFR.SS.tc5.lr002)

find.int$interactions

find.int$rank.list

gbm.perspec(UCFR.SS.tc5.lr002,1,2, 
            z.range=c(0,0.07),
            theta = 225,
            phi=45, 
            cex.lab = 1.2, font.lab = 2, cex.axis = 1, font.axis= 1,
            perspective = T,
            smooth=F)

