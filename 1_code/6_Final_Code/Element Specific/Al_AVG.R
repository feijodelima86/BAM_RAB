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


#### Al ####

Al_DF <- (COMPARTMENTS_AVG[,c("W_Al","F_Al","trulydissolved_Al","OM.AREA.g.m2","TURNOVER","Al")])

Al_DF<-Al_DF[Al_DF$OM.AREA.g.m2>0,]

Al_DF$OM.AREA.g.m2 <- log10(Al_DF$OM.AREA.g.m2+1)

Al_DF <- Al_DF[complete.cases(Al_DF),]

#Al_DF$TURNOVER <- sqrt(Al_DF$NEW)

quantile <- quantile(Al_DF$TURNOVER, probs = seq(0, 1, 1/5))

Al_DF$cat <- cut(Al_DF$TURNOVER, 
                 breaks=c(-Inf, quantile[2] , quantile[3] , quantile[4], quantile[5] , Inf)
)

### Al_SPM ###

Al_SPM <- (Al_DF[,c(1,4,5,6,7)])

names(Al_SPM)<- c("Al_SPM", "TURNOVER", "OM" ,"Al", "cat")

Al_SPM<-Al_SPM[Al_SPM$Al_SPM>0,]

Al_SPM$Al_SPM <- log10(Al_SPM$Al_SPM+1)

for (col in colnames(Al_SPM[,c(1)])) {
  Al_SPM[[col]] <- replace_outliers_with_na(Al_SPM[[col]])
}

Al_SPM <- Al_SPM[complete.cases(Al_SPM), ]

### Al_COL ###

Al_COL <- (Al_DF[,c(2,4,5,6,7)])

names(Al_COL)<- c("Al_COL", "TURNOVER", "OM","Al", "cat")

Al_COL<-Al_COL[Al_COL$Al_COL>0,]

Al_COL$Al_COL <- log10(Al_COL$Al_COL+1)

for (col in colnames(Al_COL[,c(1)])) {
  Al_COL[[col]] <- replace_outliers_with_na(Al_COL[[col]])
}

Al_COL <- Al_COL[complete.cases(Al_COL), ]

### Al_TRD ###

Al_TRD <- (Al_DF[,c(3,4,5,6,7)])

names(Al_TRD)<- c("Al_TRD", "TURNOVER", "OM" ,"Al", "cat")

Al_TRD<-Al_TRD[Al_TRD$Al_TRD>0,]

Al_TRD$Al_TRD <- log10(Al_TRD$Al_TRD+1)

for (col in colnames(Al_TRD[,c(1)])) {
  Al_TRD[[col]] <- replace_outliers_with_na(Al_TRD[[col]])
}

Al_TRD <- Al_TRD[complete.cases(Al_TRD), ]

#dev.new()
plot(density(Al_SPM$Al_SPM, bw=0.09), ylim=c(0,5), xlim=c(0,1), col="Purple", lwd=2)
lines(density(Al_COL$Al_COL, bw=0.09), col="Blue", lwd=2)
lines(density(Al_TRD$Al_TRD, bw=0.09), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Al_SPM <- glm(Al ~ (Al_SPM + TURNOVER + OM), data = Al_SPM, family=gaussian)

mix.int_Al_COL <- glm(Al ~ (Al_COL + TURNOVER + OM), data = Al_COL, family=gaussian)

mix.int_Al_TRD <- glm(Al ~ (Al_TRD + TURNOVER + OM), data = Al_TRD, family=gaussian)

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

#### % Total ####

names(Al_DF)<-c("Al_SPM","Al_COL","Al_TRD","OM","TURNOVER", "Al", "cat")
Al_DF[,c(1:3)][Al_DF[,c(1:3)] <= 0] <- 0

sum(Al_DF$Al_SPM, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))
sum(Al_DF$Al_COL, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))
sum(Al_DF$Al_TRD, na.rm=T)/(sum(Al_DF$Al_SPM,Al_DF$Al_COL,Al_DF$Al_TRD, na.rm=T))

