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


#### Co ####

Co_DF <- (COMPARTMENTS_AVG[,c("W_Co","F_Co","trulydissolved_Co","OM.AREA.g.m2","TURNOVER","Co")])

Co_DF<-Co_DF[Co_DF$OM.AREA.g.m2>0,]

Co_DF$OM.AREA.g.m2 <- log10(Co_DF$OM.AREA.g.m2+1)

Co_DF <- Co_DF[complete.cases(Co_DF),]

#Co_DF$TURNOVER <- sqrt(Co_DF$NEW)

quantile <- quantile(Co_DF$TURNOVER, probs = seq(0, 1, 1/5))

Co_DF$cat <- cut(Co_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### Co_SPM ###

Co_SPM <- (Co_DF[,c(1,4,5,6,7)])

names(Co_SPM)<- c("Co_SPM", "TURNOVER", "OM" ,"Co", "cat")

Co_SPM<-Co_SPM[Co_SPM$Co_SPM>0,]

Co_SPM$Co_SPM <- log10(Co_SPM$Co_SPM+1)

for (col in colnames(Co_SPM[,c(1)])) {
  Co_SPM[[col]] <- replace_outliers_with_na(Co_SPM[[col]])
}

Co_SPM <- Co_SPM[complete.cases(Co_SPM), ]

### Co_COL ###

Co_COL <- (Co_DF[,c(2,4,5,6,7)])

names(Co_COL)<- c("Co_COL", "TURNOVER", "OM","Co", "cat")

Co_COL<-Co_COL[Co_COL$Co_COL>0,]

Co_COL$Co_COL <- log10(Co_COL$Co_COL+1)

for (col in colnames(Co_COL[,c(1)])) {
  Co_COL[[col]] <- replace_outliers_with_na(Co_COL[[col]])
}

Co_COL <- Co_COL[complete.cases(Co_COL), ]

### Co_TRD ###

Co_TRD <- (Co_DF[,c(3,4,5,6,7)])

names(Co_TRD)<- c("Co_TRD", "TURNOVER", "OM" ,"Co", "cat")

Co_TRD<-Co_TRD[Co_TRD$Co_TRD>0,]

Co_TRD$Co_TRD <- log10(Co_TRD$Co_TRD+1)

for (col in colnames(Co_TRD[,c(1)])) {
  Co_TRD[[col]] <- replace_outliers_with_na(Co_TRD[[col]])
}

Co_TRD <- Co_TRD[complete.cases(Co_TRD), ]

#dev.new()
plot(density(Co_SPM$Co_SPM, bw=0.09), ylim=c(0,5), xlim=c(0,1), col="Purple", lwd=2)
lines(density(Co_COL$Co_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Co_TRD$Co_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Co_SPM <- glm(Co ~ (Co_SPM + TURNOVER + OM), data = Co_SPM, family=gaussian)

mix.int_Co_COL <- glm(Co ~ (Co_COL + TURNOVER + OM), data = Co_COL, family=gaussian)

mix.int_Co_TRD <- glm(Co ~ (TURNOVER), data = Co_TRD, family=gaussian)

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

#### % Total ####

names(Co_DF)<-c("Co_SPM","Co_COL","Co_TRD","OM","TURNOVER", "Al", "cat")
Co_DF[,c(1:3)][Co_DF[,c(1:3)] <= 0] <- 0

sum(Co_DF$Co_SPM, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_COL, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))
sum(Co_DF$Co_TRD, na.rm=T)/(sum(Co_DF$Co_SPM,Co_DF$Co_COL,Co_DF$Co_TRD, na.rm=T))

