library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("performance")
library("splines")

#alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
#alldata <- read_csv("2_incremental/wateralgae13_WORKING_MANUAL.csv")

#alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

#COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

#COMPARTMENTS_AVG <-aggregate(x = COMPARTMENTS[,colnames(COMPARTMENTS) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
#                             by = list(COMPARTMENTS$SAMPLING_DATE,COMPARTMENTS$SITE,COMPARTMENTS$SAMPLE_DESCRIPTOR),
#                             FUN = mean,
#                             na.rm = F)

#write.csv(COMPARTMENTS_AVG, "2_incremental/COMPARTMENTS_AVG.csv")

COMPARTMENTS <- read.csv("2_incremental/COMPARTMENTS_AVG.csv")

colnames(COMPARTMENTS)[c(2:4)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

###GLMS####

###Truly Dissolved Across Elements###

names(COMPARTMENTS)

mix.int_As <- glm(As ~ TD.As * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Ca <- glm(Car ~ TD.Ca * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cd <- glm(Cd ~ TD.Cd * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cu <- glm(Cu ~ TD.Cu * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Fe <- glm(Fer.1 ~ TD.Fe * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Mo <- glm(Mo ~ TD.Mo * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_P <- glm(P ~ TD.P * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Pb <- glm(Pb ~ TD.Pb * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Se <- glm(Se ~ TD.Se * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Zn <- glm(Znr ~ TD.Zn * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)


summary(mix.int_As)
summary(mix.int_Ca)
summary(mix.int_Cd)
summary(mix.int_Cu)
summary(mix.int_Fe)
summary(mix.int_Mo)
summary(mix.int_P)
summary(mix.int_Pb)
summary(mix.int_Se)
summary(mix.int_Zn)

dev.new()

plotGLM_As <- ggplot(mix.int_As, aes(TD.As, As, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_As

plotGLM_Ca <- ggplot(mix.int_Ca, aes(TD.Ca, Car, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Ca

plotGLM_Cd <- ggplot(mix.int_Cd, aes(TD.Cd, Cd, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Cd

plotGLM_Cu <- ggplot(mix.int_Cu, aes(TD.Cu, Cu, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Cu

plotGLM_Fe <- ggplot(mix.int_Fe, aes(TD.Fe, Fer.1, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Fe

plotGLM_Mo <- ggplot(mix.int_Mo, aes(TD.Mo, Mo, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Mo

plotGLM_P <- ggplot(mix.int_P, aes(TD.P, P, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_P

plotGLM_Pb <- ggplot(mix.int_Pb, aes(TD.Pb, Pb, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Pb

plotGLM_Se <- ggplot(mix.int_Se, aes(TD.Se, Se, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Se

plotGLM_Zn <- ggplot(mix.int_Zn, aes(TD.Zn, Znr, col = as.factor(SAMPLE_DESCRIPTOR)))+  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed")

plotGLM_Zn

###Size Fractions Across Compartments###

###Candidate Variables: As=EPIL Cd: All Cu:All Mo: All P: Epil Pb:All Se: EPIL

#### Cd  ####

### Cd SS 

EPIL.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.Cd, COMPARTMENTS$Cd, 
          xlim=c(0,max(COMPARTMENTS$SS.Cd, na.rm = T)*1.1), 
          ylim=c(0,max(COMPARTMENTS$Cd, na.rm = T)), 
     col="white",
          xlab="x", 
          ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cd, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cd, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cd, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Cd, col="chartreuse",lwd=3)
abline(EPIP.SS.Cd, col="gold",lwd=3)
abline(EPIL.SS.Cd, col="darkgreen",lwd=3)

## Cd Col

EPIL.CL.Cd<-lm(Cd ~ CL.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.CL.Cd<-lm(Cd ~ CL.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.CL.Cd<-lm(Cd ~ CL.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$CL.Cd, COMPARTMENTS$Cd, 
     xlim=c(0,max(COMPARTMENTS$CL.Cd, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Cd, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cd, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cd, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cd, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Cd, col="chartreuse",lwd=3)
abline(EPIP.CL.Cd, col="gold",lwd=3)
abline(EPIL.SS.Cd, col="darkgreen",lwd=3)


### Cd TD

EPIL.TD.Cd<-lm(Cd ~ TD.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Cd)

EPIP.TD.Cd<-lm(Cd ~ TD.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Cd)

FILA.TD.Cd<-lm(Cd ~ TD.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Cd)


plot(COMPARTMENTS$TD.Cd, COMPARTMENTS$Cd, 
     xlim=c(0,max(COMPARTMENTS$TD.Cd, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Cd, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Cd, col="chartreuse",lwd=3)
abline(EPIP.TD.Cd, col="gold",lwd=3)
abline(EPIL.TD.Cd, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cd, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cd, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Cd, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cd, pch=21, bg="darkgreen", col="black")



#### Cu  ####

### Cu SS 

EPIL.SS.Cu<-lm(Cu ~ SS.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.Cu<-lm(Cu ~ SS.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.Cu<-lm(Cu ~ SS.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.Cu, COMPARTMENTS$Cu, 
     xlim=c(0,max(COMPARTMENTS$SS.Cu, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Cu, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cu, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cu, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cu, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Cu, col="chartreuse",lwd=3)
abline(EPIP.SS.Cu, col="gold",lwd=3)
abline(EPIL.SS.Cu, col="darkgreen",lwd=3)

## Cu Col

EPIL.CL.Cu<-lm(Cu ~ CL.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.CL.Cu<-lm(Cu ~ CL.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.CL.Cu<-lm(Cu ~ CL.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$CL.Cu, COMPARTMENTS$Cu, 
     xlim=c(0,max(COMPARTMENTS$CL.Cu, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Cu, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cu, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cu, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cu, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Cu, col="chartreuse",lwd=3)
abline(EPIP.CL.Cu, col="gold",lwd=3)
abline(EPIL.SS.Cu, col="darkgreen",lwd=3)


### Cu TD

EPIL.TD.Cu<-lm(Cu ~ TD.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Cu)

EPIP.TD.Cu<-lm(Cu ~ TD.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Cu)

FILA.TD.Cu<-lm(Cu ~ TD.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Cu)


plot(COMPARTMENTS$TD.Cu, COMPARTMENTS$Cu, 
     xlim=c(0,max(COMPARTMENTS$TD.Cu, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Cu, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Cu, col="chartreuse",lwd=3)
abline(EPIP.TD.Cu, col="gold",lwd=3)
abline(EPIL.TD.Cu, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Cu, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Cu, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Cu, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Cu, pch=21, bg="darkgreen", col="black")

#### Mo  ####

### Mo SS 

EPIL.SS.Mo<-lm(Mo ~ SS.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.Mo<-lm(Mo ~ SS.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.Mo<-lm(Mo ~ SS.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.Mo, COMPARTMENTS$Mo, 
     xlim=c(0,max(COMPARTMENTS$SS.Mo, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Mo, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Mo, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Mo, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Mo, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Mo, col="chartreuse",lwd=3)
abline(EPIP.SS.Mo, col="gold",lwd=3)
abline(EPIL.SS.Mo, col="darkgreen",lwd=3)

## Mo Col

EPIL.CL.Mo<-lm(Mo ~ CL.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.Mo)

EPIP.CL.Mo<-lm(Mo ~ CL.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.Mo)

FILA.CL.Mo<-lm(Mo ~ CL.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.Mo)

plot(COMPARTMENTS$CL.Mo, COMPARTMENTS$Mo, 
     xlim=c(0,max(COMPARTMENTS$CL.Mo, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Mo, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Mo, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Mo, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Mo, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Mo, col="chartreuse",lwd=3)
abline(EPIP.CL.Mo, col="gold",lwd=3)
abline(EPIL.SS.Mo, col="darkgreen",lwd=3)


### Mo TD

EPIL.TD.Mo<-lm(Mo ~ TD.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Mo)

EPIP.TD.Mo<-lm(Mo ~ TD.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Mo)

FILA.TD.Mo<-lm(Mo ~ TD.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Mo)


plot(COMPARTMENTS$TD.Mo, COMPARTMENTS$Mo, 
     xlim=c(0,max(COMPARTMENTS$TD.Mo, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Mo, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Mo, col="chartreuse",lwd=3)
abline(EPIP.TD.Mo, col="gold",lwd=3)
abline(EPIL.TD.Mo, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Mo, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Mo, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Mo, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Mo, pch=21, bg="darkgreen", col="black")



#### Pb  ####

### Pb SS 

EPIL.SS.Pb<-lm(Pb ~ SS.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.Pb<-lm(Pb ~ SS.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.Pb<-lm(Pb ~ SS.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.Pb, COMPARTMENTS$Pb, 
     xlim=c(0,max(COMPARTMENTS$SS.Pb, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Pb, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Pb, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Pb, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Pb, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Pb, col="chartreuse",lwd=3)
abline(EPIP.SS.Pb, col="gold",lwd=3)
abline(EPIL.SS.Pb, col="darkgreen",lwd=3)

## Pb Col

EPIL.CL.Pb<-lm(Pb ~ CL.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.Pb)

EPIP.CL.Pb<-lm(Pb ~ CL.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.Pb)

FILA.CL.Pb<-lm(Pb ~ CL.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.Pb)

plot(COMPARTMENTS$CL.Pb, COMPARTMENTS$Pb, 
     xlim=c(0,max(COMPARTMENTS$CL.Pb, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Pb, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Pb, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Pb, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Pb, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Pb, col="chartreuse",lwd=3)
abline(EPIP.CL.Pb, col="gold",lwd=3)
abline(EPIL.SS.Pb, col="darkgreen",lwd=3)


### Pb TD

EPIL.TD.Pb<-lm(Pb ~ TD.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Pb)

EPIP.TD.Pb<-lm(Pb ~ TD.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Pb)

FILA.TD.Pb<-lm(Pb ~ TD.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Pb)


plot(COMPARTMENTS$TD.Pb, COMPARTMENTS$Pb, 
     xlim=c(0,max(COMPARTMENTS$TD.Pb, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Pb, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Pb, col="chartreuse",lwd=3)
abline(EPIP.TD.Pb, col="gold",lwd=3)
abline(EPIL.TD.Pb, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Pb, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Pb, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Pb, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Pb, pch=21, bg="darkgreen", col="black")

#### Fe  ####

### Fe SS 

EPIL.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.Fe, COMPARTMENTS$Fer.1, 
     xlim=c(0,max(COMPARTMENTS$SS.Fe, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Fer.1, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Fer.1, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Fer.1, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Fer.1, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Fe, col="chartreuse",lwd=3)
abline(EPIP.SS.Fe, col="gold",lwd=3)
abline(EPIL.SS.Fe, col="darkgreen",lwd=3)

## Fe Col

EPIL.CL.Fe<-lm(Fer.1 ~ CL.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.Fe)

EPIP.CL.Fe<-lm(Fer.1 ~ CL.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.Fe)

FILA.CL.Fe<-lm(Fer.1 ~ CL.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.Fe)

plot(COMPARTMENTS$CL.Fe, COMPARTMENTS$Fer.1, 
     xlim=c(0,max(COMPARTMENTS$CL.Fe, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Fer.1, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Fer.1, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Fer.1, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Fer.1, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Fe, col="chartreuse",lwd=3)
abline(EPIP.CL.Fe, col="gold",lwd=3)
abline(EPIL.SS.Fe, col="darkgreen",lwd=3)


### Fe TD

EPIL.TD.Fe<-lm(Fer.1 ~ TD.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Fe)

EPIP.TD.Fe<-lm(Fer.1 ~ TD.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Fe)

FILA.TD.Fe<-lm(Fer.1 ~ TD.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Fe)


plot(COMPARTMENTS$TD.Fe, COMPARTMENTS$Fer.1, 
     xlim=c(0,max(COMPARTMENTS$TD.Fe, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Fer.1, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Fe, col="chartreuse",lwd=3)
abline(EPIP.TD.Fe, col="gold",lwd=3)
abline(EPIL.TD.Fe, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Fer.1, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Fer.1, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Fe, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Fer.1, pch=21, bg="darkgreen", col="black")



