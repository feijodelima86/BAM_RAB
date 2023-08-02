library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))
library("ggpubr")
library("scales")
library("performance")
library("splines")
library("effects")
library("mgcv")

COMPARTMENTS <- read.csv("2_incremental/COMPARTMENTS_AVG_6.csv")

colnames(COMPARTMENTS)[c(2:4)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

COMPARTMENTS<-COMPARTMENTS[,c(2,3,4,8,10:63)]

names(COMPARTMENTS)

#COMPARTMENTS <- read_csv("2_incremental/wateralgae13_WORKING_MANUAL.csv")

###GLMS####

#### As ####

mix.int_As_TD <- glm(As ~ TD.As  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_As_CL <- glm(As ~ CL.As  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_As_SS <- glm(As ~ SS.As  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_As_TD)
summary(mix.int_As_CL)
summary(mix.int_As_SS)

#dev.new()

plotGLM_As_TD <- ggplot(mix.int_As_TD, aes(x=TD.As, y=As, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "As  (TD, \u00b5g/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.As, y=As, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_As_CL <- ggplot(mix.int_As_CL, aes(x=CL.As, y=As, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "As  (TD, \u00b5g/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.As, y=As, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_As_SS <- ggplot(mix.int_As_SS, aes(x=SS.As, y=As, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "As  (TD, \u00b5g/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.As, y=As, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


As_figure <- ggarrange(plotGLM_As_TD, plotGLM_As_CL, plotGLM_As_SS,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)

#### Cd ####

mix.int_Cd_TD <- glm(Cd ~ TD.Cd  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cd_CL <- glm(Cd ~ CL.Cd  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cd_SS <- glm(Cd ~ SS.Cd  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Cd_TD)
summary(mix.int_Cd_CL)
summary(mix.int_Cd_SS)

#dev.new()

plotGLM_Cd_TD <- ggplot(mix.int_Cd_TD, aes(x=TD.Cd, y=Cd, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cd  (TD, \u00b5g/l)",    
       y= "Cd Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Cd, y=Cd, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Cd_CL <- ggplot(mix.int_Cd_CL, aes(x=CL.Cd, y=Cd, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cd  (TD, \u00b5g/l)",    
       y= "Cd Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Cd, y=Cd, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Cd_SS <- ggplot(mix.int_Cd_SS, aes(x=SS.Cd, y=Cd, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cd  (TD, \u00b5g/l)",    
       y= "Cd Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Cd, y=Cd, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Cd <- ggarrange(plotGLM_Cd_TD, plotGLM_Cd_CL, plotGLM_Cd_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)


#### Cu ####

mix.int_Cu_TD <- glm(Cu ~ TD.Cu  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cu_CL <- glm(Cu ~ CL.Cu  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Cu_SS <- glm(Cu ~ SS.Cu  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Cu_TD)
summary(mix.int_Cu_CL)
summary(mix.int_Cu_SS)

#dev.new()

plotGLM_Cu_TD <- ggplot(mix.int_Cu_TD, aes(x=TD.Cu, y=Cu, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cu  (TD, \u00b5g/l)",    
       y= "Cu Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Cu, y=Cu, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Cu_CL <- ggplot(mix.int_Cu_CL, aes(x=CL.Cu, y=Cu, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cu  (TD, \u00b5g/l)",    
       y= "Cu Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Cu, y=Cu, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Cu_SS <- ggplot(mix.int_Cu_SS, aes(x=SS.Cu, y=Cu, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Cu  (TD, \u00b5g/l)",    
       y= "Cu Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Cu, y=Cu, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Cu <- ggarrange(plotGLM_Cu_TD, plotGLM_Cu_CL, plotGLM_Cu_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)


#### Fe ####

mix.int_Fe_TD <- glm(Fer.1 ~ TD.Fe  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian(link="log"))
mix.int_Fe_CL <- glm(Fer.1 ~ CL.Fe  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian(link="log"))
mix.int_Fe_SS <- glm(Fer.1 ~ SS.Fe  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian(link="log"))

summary(mix.int_Fe_TD)
summary(mix.int_Fe_CL)
summary(mix.int_Fe_SS)


#dev.new()

plotGLM_Fe_TD <- ggplot(mix.int_Fe_TD, aes(x=TD.Fe, y=Fer.1, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Fe  (TD, \u00b5g/l)",    
       y= "Fe Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Fe, y=Fer.1, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian"(link="log"))) +
  theme(legend.position="none")

plotGLM_Fe_CL <- ggplot(mix.int_Fe_CL, aes(x=CL.Fe, y=Fer.1, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Fe  (TD, \u00b5g/l)",    
       y= "Fe Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Fe, y=Fer.1, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian"(link="log"))) +
  theme(legend.position="none")

plotGLM_Fe_SS <- ggplot(mix.int_Fe_SS, aes(x=SS.Fe, y=Fer.1, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Fe  (TD, \u00b5g/l)",    
       y= "Fe Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Fe, y=Fer.1, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian"(link="log"))) +
  theme(legend.position="none")


figure.Fe <- ggarrange(plotGLM_Fe_TD, plotGLM_Fe_CL, plotGLM_Fe_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)

#### Mo ####

mix.int_Mo_TD <- glm(Mo ~ TD.Mo  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Mo_CL <- glm(Mo ~ CL.Mo  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Mo_SS <- glm(Mo ~ SS.Mo  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Mo_TD)
summary(mix.int_Mo_CL)
summary(mix.int_Mo_SS)


#dev.new()

plotGLM_Mo_TD <- ggplot(mix.int_Mo_TD, aes(x=TD.Mo, y=Mo, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Mo  (TD, \u00b5g/l)",    
       y= "Mo Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Mo, y=Mo, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Mo_CL <- ggplot(mix.int_Mo_CL, aes(x=CL.Mo, y=Mo, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Mo  (TD, \u00b5g/l)",    
       y= "Mo Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Mo, y=Mo, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Mo_SS <- ggplot(mix.int_Mo_SS, aes(x=SS.Mo, y=Mo, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Mo  (TD, \u00b5g/l)",    
       y= "Mo Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Mo, y=Mo, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Mo <- ggarrange(plotGLM_Mo_TD, plotGLM_Mo_CL, plotGLM_Mo_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)

#### Pb ####

mix.int_Pb_TD <- glm(Pb ~ TD.Pb  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Pb_CL <- glm(Pb ~ CL.Pb  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Pb_SS <- glm(Pb ~ SS.Pb  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Pb_TD)
summary(mix.int_Pb_CL)
summary(mix.int_Pb_SS)

#dev.new()

plotGLM_Pb_TD <- ggplot(mix.int_Pb_TD, aes(x=TD.Pb, y=Pb, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Pb  (TD, \u00b5g/l)",    
       y= "Pb Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Pb, y=Pb, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Pb_CL <- ggplot(mix.int_Pb_CL, aes(x=CL.Pb, y=Pb, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Pb  (TD, \u00b5g/l)",    
       y= "Pb Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Pb, y=Pb, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Pb_SS <- ggplot(mix.int_Pb_SS, aes(x=SS.Pb, y=Pb, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Pb  (TD, \u00b5g/l)",    
       y= "Pb Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Pb, y=Pb, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Pb <- ggarrange(plotGLM_Pb_TD, plotGLM_Pb_CL, plotGLM_Pb_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)

#### Se ####

mix.int_Se_TD <- glm(Se ~ TD.Se  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Se_CL <- glm(Se ~ CL.Se  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Se_SS <- glm(Se ~ SS.Se  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Se_TD)
summary(mix.int_Se_CL)
summary(mix.int_Se_SS)


#dev.new()

plotGLM_Se_TD <- ggplot(mix.int_Se_TD, aes(x=TD.Se, y=Se, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Se  (TD, \u00b5g/l)",    
       y= "Se Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Se, y=Se, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Se_CL <- ggplot(mix.int_Se_CL, aes(x=CL.Se, y=Se, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Se  (TD, \u00b5g/l)",    
       y= "Se Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Se, y=Se, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Se_SS <- ggplot(mix.int_Se_SS, aes(x=SS.Se, y=Se, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Se  (TD, \u00b5g/l)",    
       y= "Se Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Se, y=Se, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Se <- ggarrange(plotGLM_Se_TD, plotGLM_Se_CL, plotGLM_Se_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)

#### Zn ####

mix.int_Zn_TD <- glm(Znr ~ TD.Zn  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Zn_CL <- glm(Znr ~ CL.Zn  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)
mix.int_Zn_SS <- glm(Znr ~ SS.Zn  * SAMPLE_DESCRIPTOR, data = COMPARTMENTS, family=gaussian)

summary(mix.int_Zn_TD)
summary(mix.int_Zn_CL)
summary(mix.int_Zn_SS)


#dev.new()

plotGLM_Zn_TD <- ggplot(mix.int_Zn_TD, aes(x=TD.Zn, y=Znr, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Zn  (TD, \u00b5g/l)",    
       y= "Zn Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=TD.Zn, y=Znr, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Zn_CL <- ggplot(mix.int_Zn_CL, aes(x=CL.Zn, y=Znr, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Zn  (TD, \u00b5g/l)",    
       y= "Zn Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=CL.Zn, y=Znr, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")

plotGLM_Zn_SS <- ggplot(mix.int_Zn_SS, aes(x=SS.Zn, y=Znr, group = SAMPLE_DESCRIPTOR)) +
  geom_point(aes(color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), shape=21, size=4, stroke=2)+
  scale_shape_manual(values=c(21))+
  scale_color_manual(values=c("black","black","black"), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+
  scale_fill_manual(values=c('darkgreen','gold', 'chartreuse'), name  ="Compartment",labels=c("Eplilithon", "Epiphites", "Filamentous"))+ 
  labs(title="Title",
       x= "Zn  (TD, \u00b5g/l)",    
       y= "Zn Content (\u00b5g/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  geom_smooth(aes(x=SS.Zn, y=Znr, group = as.factor(SAMPLE_DESCRIPTOR), color=SAMPLE_DESCRIPTOR, fill=SAMPLE_DESCRIPTOR), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  theme(legend.position="none")


figure.Zn <- ggarrange(plotGLM_Zn_TD, plotGLM_Zn_CL, plotGLM_Zn_SS,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)


