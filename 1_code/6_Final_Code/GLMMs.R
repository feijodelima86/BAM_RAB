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
library(MuMIn)


## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

replace_outliers_with_na <- function(x, threshold = 2) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}


#### Cd ####

Cd_DF <- cbind(COMPARTMENTS_AVG[,c("SITE","FIELD.REP","spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","TURNOVER","Cd"),],(paste(COMPARTMENTS_AVG$SAMPLING_DATE ,COMPARTMENTS_AVG$SITE)))

Cd_DF

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

### Cd_TRD ###

Cd_TRD <- (Cd_DF[,c(1,2,5,6,7,8,9)])

names(Cd_TRD)<- c("SITE", "REP", "Cd_TRD", "TURNOVER", "Cd","DATESITE", "cat")

Cd_TRD<-Cd_TRD[Cd_TRD$Cd_TRD>0,]

Cd_TRD$Cd_TRD <- log10(Cd_TRD$Cd_TRD+1)

for (col in colnames(Cd_TRD[,c(1)])) {
  Cd_TRD[[col]] <- replace_outliers_with_na(Cd_TRD[[col]])
}

Cd_TRD <- Cd_TRD[complete.cases(Cd_TRD), ]

glmm.int_Cd_TRD <- glmer(formula = Cd ~ (Cd_TRD * TURNOVER) + (1|DATESITE),
                        data = Cd_TRD,
                        family=gaussian(link = "log"),
#                        nAGQ=25
                        )


plot(allEffects(glmm.int_Cd_TRD))

#dev.new()

par(mfrow=c(1,2))

plot(glmm.int_Cd_TRD)

plot(fitted(glmm.int_Cd_TRD),resid(glmm.int_Cd_TRD))

abline(h=0,lty=2,col="red")

qqnorm(resid(glmm.int_Cd_TRD))
qqline(resid(glmm.int_Cd_TRD))

summary(glmm.int_Cd_TRD)

r.squaredGLMM(glmm.int_Cd_TRD)

Anova(glmm.int_Cd_TRD, test="Chisq")

#### Cd Plot ####

par(mfrow=c(1,2))

n1=4
n2=5
n3=3

#dev.new()

par(mar=c(7, 7, 3, 3))

plot(Cd_TRD[,n3], Cd_TRD[,n2], 
     ylim = c(0, max(Cd_TRD[,n2],na.rm=TRUE)),
     xlim = c(min(Cd_TRD[,n3],na.rm=TRUE), max(Cd_TRD[,n3],na.rm=TRUE)),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)

atx <- round(seq(min(Cd_TRD[,n3],na.rm=TRUE), max(Cd_TRD[,n3], na.rm=TRUE), length.out=6), digits=5)
aty <- round(seq(0, max(Cd_TRD[,n2], na.rm=TRUE), length.out=6), digits=4)

Xlabel=expression(bold(paste("Cd TRD (mg/L)")))
Ylabel=expression(bold("Cd Content (mg/g)"))

axis(side = 1, at = atx, labels=format(atx, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5, family="Calibri")
title(ylab=Ylabel, line=2.7, cex.lab=1.5, family="Calibri")

pal <- colorRampPalette(c("#1b98e0", "red")) 

Cd_TRD$order = findInterval(Cd_TRD$TURNOVER, sort(Cd_TRD$TURNOVER))

points(Cd_TRD[,n3], Cd_TRD[,n2],pch=23, cex=2,col="black", bg=pal(nrow(Cd_TRD))[Cd_TRD$order],lwd=3)

box(lwd=2)

model <- lm(Cd_TRD[,n2] ~ Cd_TRD[,n3])

summary(model)

X <- range(Cd_TRD[,n3], na.rm=TRUE)
Y <- predict(model, newdata=data.frame(x=X))
lines(x=X, y=range(Y), lwd=4)


par(mar=c(7, 7, 3, 3))

plot(Cd_TRD[,n1], Cd_TRD[,n2], 
     ylim = c(0, max(Cd_TRD[,n2],na.rm=TRUE)),
     xlim = c(min(Cd_TRD[,n1]), max(Cd_TRD[,n1],na.rm=TRUE)),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)

atx <- round(seq(min(Cd_TRD[,n1]), max(Cd_TRD[,n1], na.rm=TRUE), length.out=6), digits=2)
aty <- round(seq(0, max(Cd_TRD[,n2], na.rm=TRUE), length.out=6), digits=4)

Xlabel=expression(bold(paste("Turnover (%)")))
Ylabel=expression(bold("Cd Content (mg/g)"))

axis(side = 1, at = atx, labels=format(atx, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5, family="Calibri")
title(ylab=Ylabel, line=2.7, cex.lab=1.5, family="Calibri")

pal <- colorRampPalette(c("#1b98e0", "red")) 

Cd_TRD$order = findInterval(Cd_TRD$Cd_TRD, sort(Cd_TRD$Cd_TRD))

points(Cd_TRD[,n1], Cd_TRD[,n2],pch=23, cex=2,col="black", bg=pal(nrow(Cd_TRD))[Cd_TRD$order],lwd=3)

box(lwd=2)

model <- lm(Cd_TRD[,n2] ~ Cd_TRD[,n1])

summary(model)

X <- range(Cd_TRD[,n1])
Y <- predict(model, newdata=data.frame(x=X))
lines(x=X, y=range(Y), lwd=4)

###########

N =100; set.seed(123);


x1 = runif(N)*3; readings1 <- 2*x1 + 1.0 + rnorm(N)*.99;
x2 = runif(N)*3; readings2 <- 3*x2 + 1.5 + rnorm(N)*.99;
x3 = runif(N)*3; readings3 <- 4*x3 + 2.0 + rnorm(N)*.99;
x4 = runif(N)*3; readings4 <- 5*x4 + 2.5 + rnorm(N)*.99;
x5 = runif(N)*3; readings5 <- 6*x5 + 3.0 + rnorm(N)*.99;

X = c(x1,x2,x3,x4,x5);
Y = c(readings1,readings2,readings3,readings4,readings5)
Grouping  = c(rep(1,N),rep(2,N),rep(3,N),rep(4,N),rep(5,N))

library(lme4);
LMERFIT <- lmer(Y ~ 1+ X+ (X|Grouping))

RIaS <-unlist( ranef(LMERFIT)) #Random Intercepts and Slopes
FixedEff <- fixef(LMERFIT)    # Fixed Intercept and Slope

png('SampleLMERFIT_withRandomSlopes_and_Intercepts.png', width=800,height=450,units="px" )
par(mfrow=c(1,2))
plot(X,Y,xlab="x",ylab="readings")
plot(x1,readings1, xlim=c(0,3), ylim=c(min(Y)-1,max(Y)+1), pch=16,xlab="x",ylab="readings" )
points(x2,readings2, col='red', pch=16)
points(x3,readings3, col='green', pch=16)
points(x4,readings4, col='blue', pch=16)
points(x5,readings5, col='orange', pch=16)
abline(v=(seq(-1,4 ,1)), col="lightgray", lty="dotted");        
abline(h=(seq( -1,25 ,1)), col="lightgray", lty="dotted")   

lines(x1,FixedEff[1]+ (RIaS[6] + FixedEff[2])* x1+ RIaS[1], col='black')
lines(x2,FixedEff[1]+ (RIaS[7] + FixedEff[2])* x2+ RIaS[2], col='red')
lines(x3,FixedEff[1]+ (RIaS[8] + FixedEff[2])* x3+ RIaS[3], col='green')
lines(x4,FixedEff[1]+ (RIaS[9] + FixedEff[2])* x4+ RIaS[4], col='blue')
lines(x5,FixedEff[1]+ (RIaS[10]+ FixedEff[2])* x5+ RIaS[5], col='orange') 
legend(0, 24, c("Group1","Group2","Group3","Group4","Group5" ), lty=c(1,1), col=c('black','red', 'green','blue','orange'))
dev.off()

