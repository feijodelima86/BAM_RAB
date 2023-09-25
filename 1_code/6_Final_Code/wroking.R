#### Cd ####

Cd_DF <- (COMPARTMENTS_AVG[,c("spm_Cd114","colloidal_Cd114","trulydissolved_Cd114","TURNOVER","Cd")])

Cd_DF

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

### Cd_SPM ###

Cd_SPM <- (Cd_DF[,c(1,4,5,6)])

names(Cd_SPM)<- c("Cd_SPM", "TURNOVER", "Cd", "cat")

Cd_SPM<-Cd_SPM[Cd_SPM$Cd_SPM>0,]

Cd_SPM$Cd_SPM <- log10(Cd_SPM$Cd_SPM+1)

for (col in colnames(Cd_SPM[,c(1)])) {
  Cd_SPM[[col]] <- replace_outliers_with_na(Cd_SPM[[col]])
}

Cd_SPM <- Cd_SPM[complete.cases(Cd_SPM), ]

### Cd_COL ###

Cd_COL <- (Cd_DF[,c(2,4,5,6)])

names(Cd_COL)<- c("Cd_COL", "TURNOVER", "Cd", "cat")

Cd_COL<-Cd_COL[Cd_COL$Cd_COL>0,]

Cd_COL$Cd_COL <- log10(Cd_COL$Cd_COL+1)

for (col in colnames(Cd_COL[,c(1)])) {
  Cd_COL[[col]] <- replace_outliers_with_na(Cd_COL[[col]])
}

Cd_COL <- Cd_COL[complete.cases(Cd_COL), ]

### Cd_TRD ###

Cd_TRD <- (Cd_DF[,c(3,4,5,6)])

names(Cd_TRD)<- c("Cd_TRD", "TURNOVER", "Cd", "cat")

Cd_TRD<-Cd_TRD[Cd_TRD$Cd_TRD>0,]

Cd_TRD$Cd_TRD <- log10(Cd_TRD$Cd_TRD+1)

for (col in colnames(Cd_TRD[,c(1)])) {
  Cd_TRD[[col]] <- replace_outliers_with_na(Cd_TRD[[col]])
}

Cd_TRD <- Cd_TRD[complete.cases(Cd_TRD), ]

#dev.new()

plot(density(Cd_SPM$Cd_SPM, bw=0.0029), ylim=c(0,150), xlim=c(0,0.04), col="Purple", lwd=2)
lines(density(Cd_COL$Cd_COL, bw=0.0029), col="Blue", lwd=2)
lines(density(Cd_TRD$Cd_TRD, bw=0.0029), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_Cd_SPM <- glm(Cd ~ (Cd_SPM * TURNOVER), data = Cd_SPM, family=gaussian)

mix.int_Cd_COL <- glm(Cd ~ (Cd_COL * TURNOVER), data = Cd_COL, family=gaussian)

mix.int_Cd_TRD <- glm(Cd ~ (Cd_TRD * TURNOVER), data = Cd_TRD, family=gaussian)

plot(allEffects(mix.int_Cd_SPM))

plot(allEffects(mix.int_Cd_COL))

plot(allEffects(mix.int_Cd_TRD))

Anova(mix.int_Cd_SPM, test="LR")
summary(mix.int_Cd_SPM)
with(summary(mix.int_Cd_SPM), 1 - deviance/null.deviance)

Anova(mix.int_Cd_COL, test="LR")
summary(mix.int_Cd_COL)
with(summary(mix.int_Cd_COL), 1 - deviance/null.deviance)

Anova(mix.int_Cd_TRD, test="LR")
summary(mix.int_Cd_TRD)
with(summary(mix.int_Cd_TRD), 1 - deviance/null.deviance)

### Cd Cdl Plots ###

Cd_p_SPM=ggplot(data=Cd_SPM,aes(Cd_SPM,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_SPM, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (SPM, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cd_p_COL=ggplot(data=Cd_COL,aes(Cd_COL,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_COL, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (COL, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cd_p_TRD=ggplot(data=Cd_TRD,aes(Cd_TRD,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=Cd_TRD, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (TD, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(Cd_p_SPM,Cd_p_COL,Cd_p_TRD,nrow=1)


#dev.new()