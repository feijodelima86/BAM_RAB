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

