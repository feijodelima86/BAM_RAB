#### S ####

S_DF <- (COMPARTMENTS_AVG[,c("spm_S","colloidal_S","trulydissolved_S","TURNOVER","S")])

S_DF

S_DF$TURNOVER <- asin(sqrt(S_DF$TURNOVER ))

quantile(S_DF$TURNOVER, probs = seq(0, 1, 1/5))

S_DF$cat <- cut(S_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

### S_SPM ###

S_SPM <- (S_DF[,c(1,4,5,6)])

names(S_SPM)<- c("S_SPM", "TURNOVER", "S", "cat")

S_SPM<-S_SPM[S_SPM$S_SPM>0,]

S_SPM$S_SPM <- log10(S_SPM$S_SPM+1)

for (col in colnames(S_SPM[,c(1)])) {
  S_SPM[[col]] <- replace_outliers_with_na(S_SPM[[col]])
}

S_SPM <- S_SPM[complete.cases(S_SPM), ]

### S_COL ###

S_COL <- (S_DF[,c(2,4,5,6)])

names(S_COL)<- c("S_COL", "TURNOVER", "S", "cat")

S_COL<-S_COL[S_COL$S_COL>0,]

S_COL$S_COL <- log10(S_COL$S_COL+1)

for (col in colnames(S_COL[,c(1)])) {
  S_COL[[col]] <- replace_outliers_with_na(S_COL[[col]])
}

S_COL <- S_COL[complete.cases(S_COL), ]

### S_TRD ###

S_TRD <- (S_DF[,c(3,4,5,6)])

names(S_TRD)<- c("S_TRD", "TURNOVER", "S", "cat")

S_TRD<-S_TRD[S_TRD$S_TRD>0,]

S_TRD$S_TRD <- log10(S_TRD$S_TRD+1)

for (col in colnames(S_TRD[,c(1)])) {
  S_TRD[[col]] <- replace_outliers_with_na(S_TRD[[col]])
}

S_TRD <- S_TRD[complete.cases(S_TRD), ]

#dev.new()

plot(density(S_SPM$S_SPM, bw=0.0029), ylim=c(0,150), xlim=c(0,0.04), col="Purple", lwd=2)
lines(density(S_COL$S_COL, bw=0.0029), col="Blue", lwd=2)
lines(density(S_TRD$S_TRD, bw=0.0029), col="LightBlue", lwd=2)
box(lwd=2)

mix.int_S_SPM <- glm(S ~ (TURNOVER), data = S_SPM, family=gaussian)

mix.int_S_COL <- glm(S ~ (S_COL * TURNOVER), data = S_COL, family=gaussian)

mix.int_S_TRD <- glm(S ~ (S_TRD + TURNOVER), data = S_TRD, family=gaussian)

plot(allEffects(mix.int_S_SPM))

plot(allEffects(mix.int_S_COL))

plot(allEffects(mix.int_S_TRD))

Anova(mix.int_S_SPM, test="LR")
summary(mix.int_S_SPM)
with(summary(mix.int_S_SPM), 1 - deviance/null.deviance)

Anova(mix.int_S_COL, test="LR")
summary(mix.int_S_COL)
with(summary(mix.int_S_COL), 1 - deviance/null.deviance)

Anova(mix.int_S_TRD, test="LR")
summary(mix.int_S_TRD)
with(summary(mix.int_S_TRD), 1 - deviance/null.deviance)

### S Sl Plots ###

S_p_SPM=ggplot(data=S_SPM,aes(S_SPM,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_SPM, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (SPM, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

S_p_COL=ggplot(data=S_COL,aes(S_COL,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_COL, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (COL, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

S_p_TRD=ggplot(data=S_TRD,aes(S_TRD,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=S_TRD, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (TD, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(S_p_SPM,S_p_COL,S_p_TRD,nrow=1)


#dev.new()