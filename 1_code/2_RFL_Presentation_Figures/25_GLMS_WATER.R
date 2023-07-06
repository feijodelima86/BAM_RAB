library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))
library("scales")
library("performance")
library("splines")
library("effects")
library("mgcv")

alldata <- read_csv("2_incremental/wateralgae13_WORKING_MANUAL.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG <-aggregate(x = COMPARTMENTS[,colnames(COMPARTMENTS) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
                             by = list(COMPARTMENTS$SAMPLING_DATE,COMPARTMENTS$SITE,COMPARTMENTS$SAMPLE_DESCRIPTOR),
                             FUN = mean,
                             na.rm = T)

write.csv(COMPARTMENTS_AVG, "2_incremental/COMPARTMENTS_AVG_4.csv")

COMPARTMENTS <- read.csv("2_incremental/COMPARTMENTS_AVG_5.csv")

colnames(COMPARTMENTS)[c(2:4)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

COMPARTMENTS<-COMPARTMENTS[,c(2,3,4,8,10:63)]

#COMPARTMENTS <- read_csv("2_incremental/wateralgae13_WORKING_MANUAL.csv")

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

#dev.new()

names(COMPARTMENTS)

plotGLM_As <- ggplot(mix.int_As, aes(x=TD.As, y=As, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "As  (TD, \u00b5g/l)",    
       y= "As Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
       axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.As, y=As, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_As

plotGLM_Ca <- ggplot(mix.int_Ca, aes(TD.Ca, Car, col = as.factor(SAMPLE_DESCRIPTOR)))+ 
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Ca  (TD, \u00b5g/l)",    
       y= "Ca Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Ca, y=Car, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Ca

plotGLM_Cd <- ggplot(mix.int_Cd, aes(TD.Cd, Cd, col = as.factor(SAMPLE_DESCRIPTOR)))+ 
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title=NULL,
       x= "Cd  (TD, \u00b5g/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Cd, y=Cd, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Cd

plot(allEffects(mix.int_Cd))

plotGLM_Cu <- ggplot(mix.int_Cu, aes(TD.Cu, Cu, col = as.factor(SAMPLE_DESCRIPTOR)))+ 
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title=NULL,
       x= "Cu  (TD, \u00b5g/l)",    
       y= "Cu Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Cu, y=Cu, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Cu

plotGLM_Fe <- ggplot(mix.int_Fe, aes(TD.Fe, Fer.1, col = as.factor(SAMPLE_DESCRIPTOR)))+ 
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title=NULL,
       x= "Fe  (TD, \u00b5g/l)",    
       y= "Fe Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Fe, y=Fer.1, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Fe

plotGLM_Mo <- ggplot(mix.int_Mo, aes(TD.Mo, Mo, col = as.factor(SAMPLE_DESCRIPTOR)))+ 
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  scale_x_continuous(trans='log10')+ 
  labs(title=NULL,
       x= "Mo  (TD, Log10 \u00b5g/l)",    
       y= "Mo Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Mo, y=Mo, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Mo

plotGLM_P <- ggplot(mix.int_P, aes(TD.P, P, col = as.factor(SAMPLE_DESCRIPTOR)))+  
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  scale_x_continuous(trans='log10')+ 
  labs(title=NULL,
       x= "P  (TD,Log10 \u00b5g/l)",    
       y= "P Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.P, y=P, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_P

plotGLM_Pb <- ggplot(mix.int_Pb, aes(TD.Pb, Pb, col = as.factor(SAMPLE_DESCRIPTOR)))+  
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
#  scale_x_continuous(trans='log10')+ 
  labs(title=NULL,
       x= "Pb  (TD, \u00b5g/l)",    
       y= "Pb Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Pb, y=Pb, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Pb

plotGLM_Se <- ggplot(mix.int_Se, aes(TD.Se, Se, col = as.factor(SAMPLE_DESCRIPTOR)))+  
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  #  scale_x_continuous(trans='log10')+ 
  labs(title=NULL,
       x= "Se  (TD, \u00b5g/l)",    
       y= "Se Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Se, y=Se, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 

plotGLM_Se

plotGLM_Zn <- ggplot(mix.int_Zn, aes(TD.Zn, Znr, col = as.factor(SAMPLE_DESCRIPTOR)))+  
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  #  scale_x_continuous(trans='log10')+ 
  labs(title=NULL,
       x= "Zn  (TD, \u00b5g/l)",    
       y= "Zn Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Zn, y=Znr, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) 


plotGLM_Zn

####Size Fractions Across Compartments####

###Candidate Variables: As=EPIL Cd: All Cu:All Mo: All P: Epil Pb:All Se: EPIL

#### As  ####

EPIL.SS.As<-lm(As ~ SS.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.SS.As)

EPIP.SS.As<-lm(As ~ SS.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.As)

FILA.SS.As<-lm(As ~ SS.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.As)

plot(COMPARTMENTS$SS.As, COMPARTMENTS$As, 
     xlim=c(0,max(COMPARTMENTS$SS.As, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$As, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$As, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$As, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$As, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.As, col="chartreuse",lwd=3)
abline(EPIP.SS.As, col="gold",lwd=3)
abline(EPIL.SS.As, col="darkgreen",lwd=3)

## As Col

EPIL.CL.As<-lm(As ~ CL.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.As)

EPIP.CL.As<-lm(As ~ CL.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.As)

FILA.CL.As<-lm(As ~ CL.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.As)


plot(COMPARTMENTS$CL.As, COMPARTMENTS$As, 
     xlim=c(0,max(COMPARTMENTS$CL.As, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$As, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$As, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$As, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$As, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.As, col="chartreuse",lwd=3)
abline(EPIP.CL.As, col="gold",lwd=3)
abline(EPIL.SS.As, col="darkgreen",lwd=3)


### As TD

EPIL.TD.As<-lm(As ~ TD.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.As)

EPIP.TD.As<-lm(As ~ TD.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.As)

FILA.TD.As<-lm(As ~ TD.As, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.As)


plot(COMPARTMENTS$TD.As, COMPARTMENTS$As, 
     xlim=c(0,max(COMPARTMENTS$TD.As, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$As, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.As, col="chartreuse",lwd=3)
abline(EPIP.TD.As, col="gold",lwd=3)
abline(EPIL.TD.As, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$As, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$As, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.As, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$As, pch=21, bg="darkgreen", col="black")



#### Cd  ####

### Cd SS 

EPIL.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.SS.Cd)

EPIP.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Cd)

FILA.SS.Cd<-lm(Cd ~ SS.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Cd)

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
summary(EPIL.CL.Cd)

EPIP.CL.Cd<-lm(Cd ~ CL.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIL.CL.Cd)

FILA.CL.Cd<-lm(Cd ~ CL.Cd, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.Cd)

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
summary(EPIL.SS.Cu)

EPIP.SS.Cu<-lm(Cu ~ SS.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Cu)

FILA.SS.Cu<-lm(Cu ~ SS.Cu, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Cu)


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
abline(EPIL.CL.Cu, col="darkgreen",lwd=3)


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
summary(EPIL.SS.Mo)

EPIP.SS.Mo<-lm(Mo ~ SS.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Mo)

FILA.SS.Mo<-lm(Mo ~ SS.Mo, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Mo)

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
summary(EPIL.SS.Pb)

EPIP.SS.Pb<-lm(Pb ~ SS.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Pb)

FILA.SS.Pb<-lm(Pb ~ SS.Pb, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Pb)

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


#### P ####

### P SS 

EPIL.SS.P<-lm(P ~ SS.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])

EPIP.SS.P<-lm(P ~ SS.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])

FILA.SS.P<-lm(P ~ SS.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])

plot(COMPARTMENTS$SS.P, COMPARTMENTS$P, 
     xlim=c(0,max(COMPARTMENTS$SS.P, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$P, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$P, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$P, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$P, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.P, col="chartreuse",lwd=3)
abline(EPIP.SS.P, col="gold",lwd=3)
abline(EPIL.SS.P, col="darkgreen",lwd=3)

## P Col

EPIL.CL.P<-lm(P ~ CL.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.P)

EPIP.CL.P<-lm(P ~ CL.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.P)

FILA.CL.P<-lm(P ~ CL.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.P)

plot(COMPARTMENTS$CL.P, COMPARTMENTS$P, 
     xlim=c(0,max(COMPARTMENTS$CL.P, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$P, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$P, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$P, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$P, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.P, col="chartreuse",lwd=3)
abline(EPIP.CL.P, col="gold",lwd=3)
abline(EPIL.SS.P, col="darkgreen",lwd=3)


### P TD

EPIL.TD.P<-lm(P ~ TD.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.P)

EPIP.TD.P<-lm(P ~ TD.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.P)

FILA.TD.P<-lm(P ~ TD.P, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.P)


plot(COMPARTMENTS$TD.P, COMPARTMENTS$P, 
     xlim=c(0,max(COMPARTMENTS$TD.P, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$P, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.P, col="chartreuse",lwd=3)
abline(EPIP.TD.P, col="gold",lwd=3)
abline(EPIL.TD.P, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$P, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$P, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.P, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$P, pch=21, bg="darkgreen", col="black")




#### Fe  ####

### Fe SS 

EPIL.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.SS.Fe)

EPIP.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Fe)

FILA.SS.Fe<-lm(Fer.1 ~ SS.Fe, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Fe)

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
abline(EPIL.CL.Fe, col="darkgreen",lwd=3)


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

#### Se ####

### Se SS 

EPIL.SS.Se<-lm(Se ~ SS.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.SS.Se)

EPIP.SS.Se<-lm(Se ~ SS.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.SS.Se)

FILA.SS.Se<-lm(Se ~ SS.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.SS.Se)

plot(COMPARTMENTS$SS.Se, COMPARTMENTS$Se, 
     xlim=c(0,max(COMPARTMENTS$SS.Se, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Se, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$SS.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Se, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$SS.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Se, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$SS.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Se, pch=21, bg="darkgreen", col="black")

abline(FILA.SS.Se, col="chartreuse",lwd=3)
abline(EPIP.SS.Se, col="gold",lwd=3)
abline(EPIL.SS.Se, col="darkgreen",lwd=3)

## Se Col

EPIL.CL.Se<-lm(Se ~ CL.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.CL.Se)

EPIP.CL.Se<-lm(Se ~ CL.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.CL.Se)

FILA.CL.Se<-lm(Se ~ CL.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.CL.Se)

plot(COMPARTMENTS$CL.Se, COMPARTMENTS$Se, 
     xlim=c(0,max(COMPARTMENTS$CL.Se, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Se, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$CL.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Se, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$CL.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Se, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$CL.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Se, pch=21, bg="darkgreen", col="black")

abline(FILA.CL.Se, col="chartreuse",lwd=3)
abline(EPIP.CL.Se, col="gold",lwd=3)
abline(EPIL.SS.Se, col="darkgreen",lwd=3)


### Se TD

EPIL.TD.Se<-lm(Se ~ TD.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),])
summary(EPIL.TD.Se)

EPIP.TD.Se<-lm(Se ~ TD.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),])
summary(EPIP.TD.Se)

FILA.TD.Se<-lm(Se ~ TD.Se, data=COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),])
summary(FILA.TD.Se)


plot(COMPARTMENTS$TD.Se, COMPARTMENTS$Se, 
     xlim=c(0,max(COMPARTMENTS$TD.Se, na.rm = T)*1.1), 
     ylim=c(0,max(COMPARTMENTS$Se, na.rm = T)), 
     col="white",
     xlab="x", 
     ylab="Y", 
     lwd=1.5,        
     las=1,
     font.lab=2,
     font.axis = 2,
     cex.lab=1.5,
     cex.axis=1)


abline(FILA.TD.Se, col="chartreuse",lwd=3)
abline(EPIP.TD.Se, col="gold",lwd=3)
abline(EPIL.TD.Se, col="darkgreen",lwd=3)

points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$TD.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "FILA"),]$Se, pch=21, bg="chartreuse", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$TD.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIP"),]$Se, pch=21, bg="gold", col="black")
points(COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$TD.Se, COMPARTMENTS[which(COMPARTMENTS$SAMPLE_DESCRIPTOR == "EPIL"),]$Se, pch=21, bg="darkgreen", col="black")

