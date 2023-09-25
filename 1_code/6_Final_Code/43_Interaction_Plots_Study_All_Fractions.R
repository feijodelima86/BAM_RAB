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

## Data Onborading ###

COMPARTMENTS_AVG <- data.frame(read.csv("2_incremental/TURNOVER_Full_Dataset_AVG.csv"))

COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR <-as.factor(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR)

COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIL"| COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "EPIP"),]

#COMPARTMENTS_AVG <- COMPARTMENTS_AVG[which(COMPARTMENTS_AVG$SAMPLE_DESCRIPTOR == "FILA"),]

COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'][COMPARTMENTS_AVG['SAMPLE_DESCRIPTOR'] == 'EPIP'] <- 'EPIL'

names(COMPARTMENTS_AVG)

COMPARTMENTS_AVG$SITE <- as.factor(COMPARTMENTS_AVG$SITE)

#### Al ####

# No OLRM, log scale for conc. []s  sig. Int, marg sig.

Al_DF <- (COMPARTMENTS_AVG[,c("colloidal_Al","TURNOVER","Alr")])

Al_DF <- Al_DF[complete.cases(Al_DF), ]

Al_DF$TURNOVER <- asin(sqrt(Al_DF$TURNOVER ))

quantile(Al_DF$TURNOVER, probs = seq(0, 1, 1/5))

Al_DF$cat <- cut(Al_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Al_DF)<- c("F_Al", "TURNOVER", "Al", "cat")

Al_DF<-Al_DF[Al_DF$F_Al>0,]

Al_DF$F_Al <- log10(Al_DF$F_Al)

Al_TD <- glm(Al ~ F_Al, data = Al_DF, family=gaussian)

Al_TD_2 <- glm(Al ~ TURNOVER, data = Al_DF, family=gaussian)

no.int_Al_TD  <- glm(Al ~ (F_Al + TURNOVER), data = Al_DF, family=gaussian)

mix.int_Al_TD <- glm(Al ~ (F_Al * TURNOVER), data = Al_DF, family=gaussian)

models <- list(Al_TD, Al_TD_2, mix.int_Al_TD, no.int_Al_TD)

aictab(c(models), modnames=c("Al_TD","Al_TD_2","mix.int_Al_TD","no.int_Al_TD"), method=ML)

plot(allEffects(mix.int_Al_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Al_TD))

Anova(mix.int_Al_TD)

with(summary(mix.int_Al_TD), 1 - deviance/null.deviance)

### Al all Plots ###

#dev.new()

Al_p=ggplot(data=Al_DF,aes(F_Al,Al))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Al, y=Al, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Al  (TD, mg/l)",    
       y= "Al Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Al_p

names(Al_DF)

Al_loess <- loess(Al ~ F_Al * TURNOVER, data = Al_DF)

Anova(mix.int_Al_TD)

with(summary(mix.int_Al_TD), 1 - deviance/null.deviance)


# Create a sequence of incrementally increAling (by 0.3 units) values for both wt and hp

quantile(Al_DF$F_Al, probs = seq(0, 1, 1/5))
quantile(Al_DF$TURNOVER, probs = seq(0, 1, 1/5))

Al_xgrid <-  seq(min(Al_DF$F_Al), max(Al_DF$F_Al), (max((Al_DF$F_Al-min(Al_DF$F_Al))/10)))
Al_ygrid <-  seq(min(Al_DF$TURNOVER), max(Al_DF$TURNOVER), (max((Al_DF$TURNOVER-min(Al_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Al.fit <-  expand.grid(F_Al = Al_xgrid, TURNOVER = Al_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Al_3d <-  predict(Al_loess, newdata = Al.fit)

# Abbreviated display of final matrix

Al_3d[1:4, 1:4]

# Transform data to long form
Al.melt <- melt(Al_3d, id.vars = c('F_Al', 'TURNOVER'), meAlure.vars = 'Al')

names(Al.melt) <- c('F_Al', 'TURNOVER', 'Al')

#Return data to numeric form

Al.melt$F_Al     <-  as.numeric(str_sub(Al.melt$F_Al, str_locate(Al.melt$F_Al, '=')[1,1] + 1))
Al.melt$TURNOVER    <-  as.numeric(str_sub(Al.melt$TURNOVER, str_locate(Al.melt$TURNOVER, '=')[1,1] + 1))


fig_Al <- plot_ly(Al.melt, x = ~F_Al, y = ~TURNOVER, z = ~Al, type = "contour",
                  width = 600, height = 500)


#### As ####

As_DF <- (COMPARTMENTS_AVG[,c("trulydissolved_As_Oshift","TURNOVER","As")])

As_DF <- As_DF[complete.cases(As_DF), ]

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(As_DF)<- c("F_As", "TURNOVER", "As", "cat")

As_DF<-As_DF[As_DF$F_As>0.03636364,]

As_DF$F_As <- log10(As_DF$F_As+1)

As_TD <- glm(As ~ F_As, data = As_DF, family=gaussian)

no.int_As_TD  <- glm(As ~ (F_As + TURNOVER), data = As_DF, family=gaussian)

mix.int_As_TD <- glm(As ~ (F_As * TURNOVER), data = As_DF, family=gaussian)

models <- list(As_TD, mix.int_As_TD, no.int_As_TD)

aictab(c(models), modnames=c("As_TD","mix.int_As_TD","no.int_As_TD"), method=ML)

plot(allEffects(mix.int_As_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_As_TD))

Anova(mix.int_As_TD)

with(summary(mix.int_As_TD), 1 - deviance/null.deviance)

### As Asl Plots ###

As_p=ggplot(data=As_DF,aes(F_As,As))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_As, y=As, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "As  (TD, mg/l)",    
       y= "As Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

As_p

names(As_DF)

As_loess <- loess(As ~ F_As * TURNOVER, data = As_DF)

# Create a sequence of incrementAsly increAsing (by 0.3 units) vAsues for both wt and hp

quantile(As_DF$F_As, probs = seq(0, 1, 1/5))
quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_xgrid <-  seq(min(As_DF$F_As)*0.9, max(As_DF$F_As)*1.1, (max((As_DF$F_As-min(As_DF$F_As))/10)))
As_ygrid <-  seq(min(As_DF$TURNOVER)*0.9, max(As_DF$TURNOVER)*1.1, (max((As_DF$TURNOVER-min(As_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

As.fit <-  expand.grid(F_As = As_xgrid, TURNOVER = As_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

As_3d <-  predict(As_loess, newdata = As.fit)

# Abbreviated display of finAs matrix

As_3d[1:4, 1:4]

# Transform data to long form
As.melt <- melt(As_3d, id.vars = c('F_As', 'TURNOVER'), meAsure.vars = 'As')

names(As.melt) <- c('F_As', 'TURNOVER', 'As')

#Return data to numeric form

As.melt$F_As     <-  as.numeric(str_sub(As.melt$F_As, str_locate(As.melt$F_As, '=')[1,1] + 1))
As.melt$TURNOVER    <-  as.numeric(str_sub(As.melt$TURNOVER, str_locate(As.melt$TURNOVER, '=')[1,1] + 1))


fig_As <- plot_ly(As.melt, x = ~F_As, y = ~TURNOVER, z = ~As, type = "contour",
                  width = 600, height = 500) %>%
  add_trace(x=~As_DF$F_As, y=~As_DF$TURNOVER, size=~As_DF$As, type = "scatter", mode="markers")


#### Cd ####

Cd_DF <- (COMPARTMENTS_AVG[,c("F_Cd114","TURNOVER","Cd")])

Cd_DF <- Cd_DF[complete.cases(Cd_DF), ]

Cd_DF$TURNOVER <- asin(sqrt(Cd_DF$TURNOVER ))

quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_DF$cat <- cut(Cd_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Cd_DF)<- c("F_Cd", "TURNOVER", "Cd", "cat")

Cd_TD <- glm(Cd ~ F_Cd, data = Cd_DF, family=gaussian)

no.int_Cd_TD  <- glm(Cd ~ (F_Cd + TURNOVER), data = Cd_DF, family=gaussian)

mix.int_Cd_TD <- glm(Cd ~ (F_Cd * TURNOVER), data = Cd_DF, family=gaussian)

models <- list(Cd_TD, mix.int_Cd_TD, no.int_Cd_TD)

aictab(c(models), modnames=c("Cd_TD","mix.int_Cd_TD","no.int_Cd_TD"), method=ML)

plot(allEffects(mix.int_Cd_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Cd_TD))

Anova(mix.int_Cd_TD)

with(summary(mix.int_Cd_TD), 1 - deviance/null.deviance)

### Cd Cdl Plots ###

Cd_p=ggplot(data=Cd_DF,aes(F_Cd,Cd))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Cd, y=Cd, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cd  (TD, mg/l)",    
       y= "Cd Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cd_p

names(Cd_DF)

Cd_loess <- loess(Cd ~ F_Cd * TURNOVER, data = Cd_DF)

# Create a sequence of incrementCdly increCding (by 0.3 units) vCdues for both wt and hp

quantile(Cd_DF$F_Cd, probs = seq(0, 1, 1/5))
quantile(Cd_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cd_xgrid <-  seq(min(Cd_DF$F_Cd)*0.9, max(Cd_DF$F_Cd)*1.1, (max((Cd_DF$F_Cd-min(Cd_DF$F_Cd))/10)))
Cd_ygrid <-  seq(min(Cd_DF$TURNOVER)*0.9, max(Cd_DF$TURNOVER)*1.1, (max((Cd_DF$TURNOVER-min(Cd_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Cd.fit <-  expand.grid(F_Cd = Cd_xgrid, TURNOVER = Cd_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Cd_3d <-  predict(Cd_loess, newdata = Cd.fit)

# Abbreviated display of finCd matrix

Cd_3d[1:4, 1:4]

# Transform data to long form
Cd.melt <- melt(Cd_3d, id.vars = c('F_Cd', 'TURNOVER'), meCdure.vars = 'Cd')

names(Cd.melt) <- c('F_Cd', 'TURNOVER', 'Cd')

#Return data to numeric form

Cd.melt$F_Cd     <-  as.numeric(str_sub(Cd.melt$F_Cd, str_locate(Cd.melt$F_Cd, '=')[1,1] + 1))
Cd.melt$TURNOVER    <-  as.numeric(str_sub(Cd.melt$TURNOVER, str_locate(Cd.melt$TURNOVER, '=')[1,1] + 1))


fig_Cd <- plot_ly(Cd.melt, x = ~F_Cd, y = ~TURNOVER, z = ~Cd, type = "contour",
                  width = 600, height = 500) %>%
  add_trace(x=~Cd_DF$F_Cd, y=~Cd_DF$TURNOVER, size=~Cd_DF$Cd, type = "scatter", mode="markers")


#### Cu ####

Cu_DF <- (COMPARTMENTS_AVG[,c("F_Cu_NH3","TURNOVER","Cu")])

Cu_DF <- Cu_DF[complete.cases(Cu_DF), ]

Cu_DF$TURNOVER <- asin(sqrt(Cu_DF$TURNOVER ))

quantile(Cu_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cu_DF$cat <- cut(Cu_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Cu_DF)<- c("F_Cu", "TURNOVER", "Cu", "cat")

Cu_DF$F_Cu <- log10(Cu_DF$F_Cu+1)

Cu_TD <- glm(Cu ~ F_Cu, data = Cu_DF, family=gaussian)

no.int_Cu_TD  <- glm(Cu ~ (F_Cu + TURNOVER), data = Cu_DF, family=gaussian)

mix.int_Cu_TD <- glm(Cu ~ (F_Cu * TURNOVER), data = Cu_DF, family=gaussian)

models <- list(Cu_TD, mix.int_Cu_TD, no.int_Cu_TD)

aictab(c(models), modnames=c("Cu_TD","mix.int_Cu_TD","no.int_Cu_TD"), method=ML)

plot(allEffects(mix.int_Cu_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Cu_TD))

Anova(mix.int_Cu_TD)

with(summary(mix.int_Cu_TD), 1 - deviance/null.deviance)

### Cu Cul Plots ###

Cu_p=ggplot(data=Cu_DF,aes(F_Cu,Cu))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Cu, y=Cu, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cu  (TD, mg/l)",    
       y= "Cu Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cu_p

names(Cu_DF)

Cu_loess <- loess(Cu ~ F_Cu * TURNOVER, data = Cu_DF)

# Create a sequence of incrementCuly increCuing (by 0.3 units) vCuues for both wt and hp

quantile(Cu_DF$F_Cu, probs = seq(0, 1, 1/5))
quantile(Cu_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cu_xgrid <-  seq(min(Cu_DF$F_Cu)*0.9, max(Cu_DF$F_Cu)*1.1, (max((Cu_DF$F_Cu-min(Cu_DF$F_Cu))/10)))
Cu_ygrid <-  seq(min(Cu_DF$TURNOVER)*0.9, max(Cu_DF$TURNOVER)*1.1, (max((Cu_DF$TURNOVER-min(Cu_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Cu.fit <-  expand.grid(F_Cu = Cu_xgrid, TURNOVER = Cu_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Cu_3d <-  predict(Cu_loess, newdata = Cu.fit)

# Abbreviated display of finCu matrix

Cu_3d[1:4, 1:4]

# Transform data to long form
Cu.melt <- melt(Cu_3d, id.vars = c('F_Cu', 'TURNOVER'), meCuure.vars = 'Cu')

names(Cu.melt) <- c('F_Cu', 'TURNOVER', 'Cu')

#Return data to numeric form

Cu.melt$F_Cu     <-  as.numeric(str_sub(Cu.melt$F_Cu, str_locate(Cu.melt$F_Cu, '=')[1,1] + 1))
Cu.melt$TURNOVER    <-  as.numeric(str_sub(Cu.melt$TURNOVER, str_locate(Cu.melt$TURNOVER, '=')[1,1] + 1))


fig_Cu <- plot_ly(Cu.melt, x = ~F_Cu, y = ~TURNOVER, z = ~Cu, type = "contour",
                  width = 600, height = 500)  %>%
  add_trace(x=~Cu_DF$F_Cu, y=~Cu_DF$TURNOVER, size=~Cu_DF$Cu, type = "scatter", mode="markers")


#### Fe ####

Fe_DF <- (COMPARTMENTS_AVG[,c("F_Fe_NH3","TURNOVER","Fe")])

Fe_DF <- Fe_DF[complete.cases(Fe_DF), ]

Fe_DF$TURNOVER <- asin(sqrt(Fe_DF$TURNOVER ))

quantile(Fe_DF$TURNOVER, probs = seq(0, 1, 1/5))

Fe_DF$cat <- cut(Fe_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Fe_DF)<- c("F_Fe", "TURNOVER", "Fe", "cat")

Fe_TD <- glm(Fe ~ F_Fe, data = Fe_DF, family=gaussian)

no.int_Fe_TD  <- glm(Fe ~ (F_Fe + TURNOVER), data = Fe_DF, family=gaussian)

mix.int_Fe_TD <- glm(Fe ~ (F_Fe * TURNOVER), data = Fe_DF, family=gaussian)

models <- list(Fe_TD, mix.int_Fe_TD, no.int_Fe_TD)

aictab(c(models), modnames=c("Fe_TD","mix.int_Fe_TD","no.int_Fe_TD"), method=ML)

plot(allEffects(mix.int_Fe_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Fe_TD))

Anova(mix.int_Fe_TD)

with(summary(mix.int_Fe_TD), 1 - deviance/null.deviance)

### Fe Fel Plots ###

Fe_p=ggplot(data=Fe_DF,aes(F_Fe,Fe))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Fe, y=Fe, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Fe  (TD, mg/l)",    
       y= "Fe Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Fe_p

names(Fe_DF)

Fe_loess <- loess(Fe ~ F_Fe * TURNOVER, data = Fe_DF)

# Create a sequence of incrementFely increFeing (by 0.3 units) vFeues for both wt and hp

quantile(Fe_DF$F_Fe, probs = seq(0, 1, 1/5))
quantile(Fe_DF$TURNOVER, probs = seq(0, 1, 1/5))

Fe_xgrid <-  seq(min(Fe_DF$F_Fe), max(Fe_DF$F_Fe), (max((Fe_DF$F_Fe-min(Fe_DF$F_Fe))/10)))
Fe_ygrid <-  seq(min(Fe_DF$TURNOVER), max(Fe_DF$TURNOVER), (max((Fe_DF$TURNOVER-min(Fe_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Fe.fit <-  expand.grid(F_Fe = Fe_xgrid, TURNOVER = Fe_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Fe_3d <-  predict(Fe_loess, newdata = Fe.fit)

# Abbreviated display of finFe matrix

Fe_3d[1:4, 1:4]

# Transform data to long form
Fe.melt <- melt(Fe_3d, id.vars = c('F_Fe', 'TURNOVER'), meFeure.vars = 'Fe')

names(Fe.melt) <- c('F_Fe', 'TURNOVER', 'Fe')

#Return data to numeric form

Fe.melt$F_Fe     <-  as.numeric(str_sub(Fe.melt$F_Fe, str_locate(Fe.melt$F_Fe, '=')[1,1] + 1))
Fe.melt$TURNOVER    <-  as.numeric(str_sub(Fe.melt$TURNOVER, str_locate(Fe.melt$TURNOVER, '=')[1,1] + 1))


fig_Fe <- plot_ly(Fe.melt, x = ~F_Fe, y = ~TURNOVER, z = ~Fe, type = "contour",
                  width = 600, height = 500)


#### Mo ####

Mo_DF <- (COMPARTMENTS_AVG[,c("F_Mo_Oshift","TURNOVER","Mo")])

Mo_DF <- Mo_DF[complete.cases(Mo_DF), ]

Mo_DF$TURNOVER <- asin(sqrt(Mo_DF$TURNOVER ))

quantile(Mo_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mo_DF$cat <- cut(Mo_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Mo_DF)<- c("F_Mo", "TURNOVER", "Mo", "cat")

Mo_TD <- glm(Mo ~ F_Mo, data = Mo_DF, family=gaussian)

no.int_Mo_TD  <- glm(Mo ~ (F_Mo + TURNOVER), data = Mo_DF, family=gaussian)

mix.int_Mo_TD <- glm(Mo ~ (F_Mo * TURNOVER), data = Mo_DF, family=gaussian)

models <- list(Mo_TD, mix.int_Mo_TD, no.int_Mo_TD)

aictab(c(models), modnames=c("Mo_TD","mix.int_Mo_TD","no.int_Mo_TD"), method=ML)

plot(allEffects(mix.int_Mo_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Mo_TD))

Anova(mix.int_Mo_TD)

with(summary(mix.int_Mo_TD), 1 - deviance/null.deviance)

### Mo Mol Plots ###

Mo_p=ggplot(data=Mo_DF,aes(F_Mo,Mo))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Mo, y=Mo, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mo  (TD, mg/l)",    
       y= "Mo Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mo_p

names(Mo_DF)

Mo_loess <- loess(Mo ~ F_Mo * TURNOVER, data = Mo_DF)

# Create a sequence of incrementMoly increMoing (by 0.3 units) vMoues for both wt and hp

quantile(Mo_DF$F_Mo, probs = seq(0, 1, 1/5))
quantile(Mo_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mo_xgrid <-  seq(min(Mo_DF$F_Mo), max(Mo_DF$F_Mo), (max((Mo_DF$F_Mo-min(Mo_DF$F_Mo))/10)))
Mo_ygrid <-  seq(min(Mo_DF$TURNOVER), max(Mo_DF$TURNOVER), (max((Mo_DF$TURNOVER-min(Mo_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Mo.fit <-  expand.grid(F_Mo = Mo_xgrid, TURNOVER = Mo_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Mo_3d <-  predict(Mo_loess, newdata = Mo.fit)

# Abbreviated display of finMo matrix

Mo_3d[1:4, 1:4]

# Transform data to long form
Mo.melt <- melt(Mo_3d, id.vars = c('F_Mo', 'TURNOVER'), meMoure.vars = 'Mo')

names(Mo.melt) <- c('F_Mo', 'TURNOVER', 'Mo')

#Return data to numeric form

Mo.melt$F_Mo     <-  as.numeric(str_sub(Mo.melt$F_Mo, str_locate(Mo.melt$F_Mo, '=')[1,1] + 1))
Mo.melt$TURNOVER    <-  as.numeric(str_sub(Mo.melt$TURNOVER, str_locate(Mo.melt$TURNOVER, '=')[1,1] + 1))


fig_Mo <- plot_ly(Mo.melt, x = ~F_Mo, y = ~TURNOVER, z = ~Mo, type = "contour",
                  width = 600, height = 500)



#### P ####

P_DF <- (COMPARTMENTS_AVG[,c("F_P","TURNOVER","P")])

P_DF <- P_DF[complete.cases(P_DF), ]

P_DF$TURNOVER <- asin(sqrt(P_DF$TURNOVER ))

quantile(P_DF$TURNOVER, probs = seq(0, 1, 1/5))

P_DF$cat <- cut(P_DF$TURNOVER, 
                breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(P_DF)<- c("F_P", "TURNOVER", "P", "cat")

P_TD <- glm(P ~ F_P, data = P_DF, family=gaussian)

no.int_P_TD  <- glm(P ~ (F_P + TURNOVER), data = P_DF, family=gaussian)

mix.int_P_TD <- glm(P ~ (F_P * TURNOVER), data = P_DF, family=gaussian)

models <- list(P_TD, mix.int_P_TD, no.int_P_TD)

aictab(c(models), modnames=c("P_TD","mix.int_P_TD","no.int_P_TD"), method=ML)

plot(allEffects(mix.int_P_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_P_TD))

Anova(mix.int_P_TD)

with(summary(mix.int_P_TD), 1 - deviance/null.deviance)

### P Pl Plots ###

P_p=ggplot(data=P_DF,aes(F_P,P))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_P, y=P, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "P  (TD, mg/l)",    
       y= "P Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

P_p

names(P_DF)

P_loess <- loess(P ~ F_P * TURNOVER, data = P_DF)

# Create a sequence of incrementPly increPing (by 0.3 units) vPues for both wt and hp

quantile(P_DF$F_P, probs = seq(0, 1, 1/5))
quantile(P_DF$TURNOVER, probs = seq(0, 1, 1/5))

P_xgrid <-  seq(min(P_DF$F_P), max(P_DF$F_P), (max((P_DF$F_P-min(P_DF$F_P))/10)))
P_ygrid <-  seq(min(P_DF$TURNOVER), max(P_DF$TURNOVER), (max((P_DF$TURNOVER-min(P_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

P.fit <-  expand.grid(F_P = P_xgrid, TURNOVER = P_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

P_3d <-  predict(P_loess, newdata = P.fit)

# Abbreviated display of finP matrix

P_3d[1:4, 1:4]

# Transform data to long form
P.melt <- melt(P_3d, id.vars = c('F_P', 'TURNOVER'), mePure.vars = 'P')

names(P.melt) <- c('F_P', 'TURNOVER', 'P')

#Return data to numeric form

P.melt$F_P     <-  as.numeric(str_sub(P.melt$F_P, str_locate(P.melt$F_P, '=')[1,1] + 1))
P.melt$TURNOVER    <-  as.numeric(str_sub(P.melt$TURNOVER, str_locate(P.melt$TURNOVER, '=')[1,1] + 1))


fig_P <- plot_ly(P.melt, x = ~F_P, y = ~TURNOVER, z = ~P, type = "contour",
                 width = 600, height = 500)


#### Pb ####

Pb_DF <- (COMPARTMENTS_AVG[,c("UF_Pb_NH3","TURNOVER","Pb")])

Pb_DF <- Pb_DF[complete.cases(Pb_DF), ]

Pb_DF$TURNOVER <- asin(sqrt(Pb_DF$TURNOVER ))

quantile(Pb_DF$TURNOVER, probs = seq(0, 1, 1/5))

Pb_DF$cat <- cut(Pb_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Pb_DF)<- c("F_Pb", "TURNOVER", "Pb", "cat")

Pb_TD <- glm(Pb ~ F_Pb, data = Pb_DF, family=gaussian)

no.int_Pb_TD  <- glm(Pb ~ (F_Pb + TURNOVER), data = Pb_DF, family=gaussian)

mix.int_Pb_TD <- glm(Pb ~ (F_Pb * TURNOVER), data = Pb_DF, family=gaussian)

models <- list(Pb_TD, mix.int_Pb_TD, no.int_Pb_TD)

aictab(c(models), modnames=c("Pb_TD","mix.int_Pb_TD","no.int_Pb_TD"), method=ML)

plot(allEffects(mix.int_Pb_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Pb_TD))

Anova(mix.int_Pb_TD)

with(summary(mix.int_Pb_TD), 1 - deviance/null.deviance)

### Pb Pbl Plots ###

Pb_p=ggplot(data=Pb_DF,aes(F_Pb,Pb))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Pb, y=Pb, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Pb  (TD, mg/l)",    
       y= "Pb Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Pb_p

names(Pb_DF)

Pb_loess <- loess(Pb ~ F_Pb * TURNOVER, data = Pb_DF)

# Create a sequence of incrementPbly increPbing (by 0.3 units) vPbues for both wt and hp

quantile(Pb_DF$F_Pb, probs = seq(0, 1, 1/5))
quantile(Pb_DF$TURNOVER, probs = seq(0, 1, 1/5))

Pb_xgrid <-  seq(min(Pb_DF$F_Pb), max(Pb_DF$F_Pb), (max((Pb_DF$F_Pb-min(Pb_DF$F_Pb))/10)))
Pb_ygrid <-  seq(min(Pb_DF$TURNOVER), max(Pb_DF$TURNOVER), (max((Pb_DF$TURNOVER-min(Pb_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Pb.fit <-  expand.grid(F_Pb = Pb_xgrid, TURNOVER = Pb_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Pb_3d <-  predict(Pb_loess, newdata = Pb.fit)

# Abbreviated display of finPb matrix

Pb_3d[1:4, 1:4]

# Transform data to long form
Pb.melt <- melt(Pb_3d, id.vars = c('F_Pb', 'TURNOVER'), mePbure.vars = 'Pb')

names(Pb.melt) <- c('F_Pb', 'TURNOVER', 'Pb')

#Return data to numeric form

Pb.melt$F_Pb     <-  as.numeric(str_sub(Pb.melt$F_Pb, str_locate(Pb.melt$F_Pb, '=')[1,1] + 1))
Pb.melt$TURNOVER    <-  as.numeric(str_sub(Pb.melt$TURNOVER, str_locate(Pb.melt$TURNOVER, '=')[1,1] + 1))


fig_Pb <- plot_ly(Pb.melt, x = ~F_Pb, y = ~TURNOVER, z = ~Pb, type = "contour",
                  width = 600, height = 500)


#### Se ####

Se_DF <- (COMPARTMENTS_AVG[,c("F_Se82_Oshift","TURNOVER","Se")])

Se_DF <- Se_DF[complete.cases(Se_DF), ]

Se_DF$TURNOVER <- asin(sqrt(Se_DF$TURNOVER ))

quantile(Se_DF$TURNOVER, probs = seq(0, 1, 1/5))

Se_DF$cat <- cut(Se_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Se_DF)<- c("F_Se", "TURNOVER", "Se", "cat")

Se_TD <- glm(Se ~ F_Se, data = Se_DF, family=gaussian)

no.int_Se_TD  <- glm(Se ~ (F_Se + TURNOVER), data = Se_DF, family=gaussian)

mix.int_Se_TD <- glm(Se ~ (F_Se * TURNOVER), data = Se_DF, family=gaussian)

models <- list(Se_TD, mix.int_Se_TD, no.int_Se_TD)

aictab(c(models), modnames=c("Se_TD","mix.int_Se_TD","no.int_Se_TD"), method=ML)

plot(allEffects(mix.int_Se_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Se_TD))

Anova(mix.int_Se_TD)

with(summary(mix.int_Se_TD), 1 - deviance/null.deviance)

### Se Sel Plots ###

Se_p=ggplot(data=Se_DF,aes(F_Se,Se))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Se, y=Se, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Se  (TD, mg/l)",    
       y= "Se Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Se_p

names(Se_DF)

Se_loess <- loess(Se ~ F_Se * TURNOVER, data = Se_DF)

# Create a sequence of incrementSely increSeing (by 0.3 units) vSeues for both wt and hp

quantile(Se_DF$F_Se, probs = seq(0, 1, 1/5))
quantile(Se_DF$TURNOVER, probs = seq(0, 1, 1/5))

Se_xgrid <-  seq(min(Se_DF$F_Se), max(Se_DF$F_Se), (max((Se_DF$F_Se-min(Se_DF$F_Se))/10)))
Se_ygrid <-  seq(min(Se_DF$TURNOVER), max(Se_DF$TURNOVER), (max((Se_DF$TURNOVER-min(Se_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Se.fit <-  expand.grid(F_Se = Se_xgrid, TURNOVER = Se_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Se_3d <-  predict(Se_loess, newdata = Se.fit)

# Abbreviated display of finSe matrix

Se_3d[1:4, 1:4]

# Transform data to long form
Se.melt <- melt(Se_3d, id.vars = c('F_Se', 'TURNOVER'), meSeure.vars = 'Se')

names(Se.melt) <- c('F_Se', 'TURNOVER', 'Se')

#Return data to numeric form

Se.melt$F_Se     <-  as.numeric(str_sub(Se.melt$F_Se, str_locate(Se.melt$F_Se, '=')[1,1] + 1))
Se.melt$TURNOVER    <-  as.numeric(str_sub(Se.melt$TURNOVER, str_locate(Se.melt$TURNOVER, '=')[1,1] + 1))


fig_Se <- plot_ly(Se.melt, x = ~F_Se, y = ~TURNOVER, z = ~Se, type = "contour",
                  width = 600, height = 500)


#### Zn ####

Zn_DF <- (COMPARTMENTS_AVG[,c("F_Zn","TURNOVER","Zn")])

Zn_DF <- Zn_DF[complete.cases(Zn_DF), ]

Zn_DF$TURNOVER <- asin(sqrt(Zn_DF$TURNOVER ))

quantile(Zn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Zn_DF$cat <- cut(Zn_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Zn_DF)<- c("F_Zn", "TURNOVER", "Zn", "cat")

Zn_TD <- glm(Zn ~ F_Zn, data = Zn_DF, family=gaussian)

no.int_Zn_TD  <- glm(Zn ~ (F_Zn + TURNOVER), data = Zn_DF, family=gaussian)

mix.int_Zn_TD <- glm(Zn ~ (F_Zn * TURNOVER), data = Zn_DF, family=gaussian)

models <- list(Zn_TD, mix.int_Zn_TD, no.int_Zn_TD)

aictab(c(models), modnames=c("Zn_TD","mix.int_Zn_TD","no.int_Zn_TD"), method=ML)

plot(allEffects(mix.int_Zn_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Zn_TD))

Anova(mix.int_Zn_TD)

with(summary(mix.int_Zn_TD), 1 - deviance/null.deviance)

### Zn Znl Plots ###

Zn_p=ggplot(data=Zn_DF,aes(F_Zn,Zn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Zn, y=Zn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Zn  (TD, mg/l)",    
       y= "Zn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Zn_p

names(Zn_DF)

Zn_loess <- loess(Zn ~ F_Zn * TURNOVER, data = Zn_DF)

# Create a sequence of incrementZnly increZning (by 0.3 units) vZnues for both wt and hp

quantile(Zn_DF$F_Zn, probs = seq(0, 1, 1/5))
quantile(Zn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Zn_xgrid <-  seq(min(Zn_DF$F_Zn), max(Zn_DF$F_Zn), (max((Zn_DF$F_Zn-min(Zn_DF$F_Zn))/10)))
Zn_ygrid <-  seq(min(Zn_DF$TURNOVER), max(Zn_DF$TURNOVER), (max((Zn_DF$TURNOVER-min(Zn_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Zn.fit <-  expand.grid(F_Zn = Zn_xgrid, TURNOVER = Zn_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Zn_3d <-  predict(Zn_loess, newdata = Zn.fit)

# Abbreviated display of finZn matrix

Zn_3d[1:4, 1:4]

# Transform data to long form
Zn.melt <- melt(Zn_3d, id.vars = c('F_Zn', 'TURNOVER'), meZnure.vars = 'Zn')

names(Zn.melt) <- c('F_Zn', 'TURNOVER', 'Zn')

#Return data to numeric form

Zn.melt$F_Zn     <-  as.numeric(str_sub(Zn.melt$F_Zn, str_locate(Zn.melt$F_Zn, '=')[1,1] + 1))
Zn.melt$TURNOVER    <-  as.numeric(str_sub(Zn.melt$TURNOVER, str_locate(Zn.melt$TURNOVER, '=')[1,1] + 1))


fig_Zn <- plot_ly(Zn.melt, x = ~F_Zn, y = ~TURNOVER, z = ~Zn, type = "contour",
                  width = 600, height = 500)





#### Ca ####

Ca_DF <- (COMPARTMENTS_AVG[,c("F_Ca","TURNOVER","Car")])

Ca_DF <- Ca_DF[complete.cases(Ca_DF), ]

Ca_DF$TURNOVER <- asin(sqrt(Ca_DF$TURNOVER ))

quantile(Ca_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ca_DF$cat <- cut(Ca_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Ca_DF)<- c("F_Ca", "TURNOVER", "Ca", "cat")

Ca_TD <- glm(Ca ~ F_Ca, data = Ca_DF, family=gaussian)

no.int_Ca_TD  <- glm(Ca ~ (F_Ca + TURNOVER), data = Ca_DF, family=gaussian)

mix.int_Ca_TD <- glm(Ca ~ (F_Ca * TURNOVER), data = Ca_DF, family=gaussian)

models <- list(Ca_TD, mix.int_Ca_TD, no.int_Ca_TD)

aictab(c(models), modnames=c("Ca_TD","mix.int_Ca_TD","no.int_Ca_TD"), method=ML)

plot(allEffects(mix.int_Ca_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Ca_TD))

Anova(mix.int_Ca_TD)

with(summary(mix.int_Ca_TD), 1 - deviance/null.deviance)

### Ca Cal Plots ###

Ca_p=ggplot(data=Ca_DF,aes(F_Ca,Ca))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Ca, y=Ca, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ca  (TD, mg/l)",    
       y= "Ca Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ca_p

names(Ca_DF)

Ca_loess <- loess(Ca ~ F_Ca * TURNOVER, data = Ca_DF)

# Create a sequence of incrementCaly increCaing (by 0.3 units) vCaues for both wt and hp

quantile(Ca_DF$F_Ca, probs = seq(0, 1, 1/5))
quantile(Ca_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ca_xgrid <-  seq(min(Ca_DF$F_Ca), max(Ca_DF$F_Ca), (max((Ca_DF$F_Ca-min(Ca_DF$F_Ca))/10)))
Ca_ygrid <-  seq(min(Ca_DF$TURNOVER), max(Ca_DF$TURNOVER), (max((Ca_DF$TURNOVER-min(Ca_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Ca.fit <-  expand.grid(F_Ca = Ca_xgrid, TURNOVER = Ca_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Ca_3d <-  predict(Ca_loess, newdata = Ca.fit)

# Abbreviated display of finCa matrix

Ca_3d[1:4, 1:4]

# Transform data to long form
Ca.melt <- melt(Ca_3d, id.vars = c('F_Ca', 'TURNOVER'), meCaure.vars = 'Ca')

names(Ca.melt) <- c('F_Ca', 'TURNOVER', 'Ca')

#Return data to numeric form

Ca.melt$F_Ca     <-  as.numeric(str_sub(Ca.melt$F_Ca, str_locate(Ca.melt$F_Ca, '=')[1,1] + 1))
Ca.melt$TURNOVER    <-  as.numeric(str_sub(Ca.melt$TURNOVER, str_locate(Ca.melt$TURNOVER, '=')[1,1] + 1))


fig_Ca <- plot_ly(Ca.melt, x = ~F_Ca, y = ~TURNOVER, z = ~Ca, type = "contour",
                  width = 600, height = 500)


#### Co ####

Co_DF <- (COMPARTMENTS_AVG[,c("F_Co","TURNOVER","Co")])

Co_DF <- Co_DF[complete.cases(Co_DF), ]

Co_DF$TURNOVER <- asin(sqrt(Co_DF$TURNOVER ))

quantile(Co_DF$TURNOVER, probs = seq(0, 1, 1/5))

Co_DF$cat <- cut(Co_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Co_DF)<- c("F_Co", "TURNOVER", "Co", "cat")

Co_TD <- glm(Co ~ F_Co, data = Co_DF, family=gaussian)

no.int_Co_TD  <- glm(Co ~ (F_Co + TURNOVER), data = Co_DF, family=gaussian)

mix.int_Co_TD <- glm(Co ~ (F_Co * TURNOVER), data = Co_DF, family=gaussian)

models <- list(Co_TD, mix.int_Co_TD, no.int_Co_TD)

aictab(c(models), modnames=c("Co_TD","mix.int_Co_TD","no.int_Co_TD"), method=ML)

plot(allEffects(mix.int_Co_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Co_TD))

Anova(mix.int_Co_TD)

with(summary(mix.int_Co_TD), 1 - deviance/null.deviance)

### Co Col Plots ###

Co_p=ggplot(data=Co_DF,aes(F_Co,Co))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Co, y=Co, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Co  (TD, mg/l)",    
       y= "Co Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Co_p

names(Co_DF)

Co_loess <- loess(Co ~ F_Co * TURNOVER, data = Co_DF)

# Create a sequence of incrementColy increCoing (by 0.3 units) vCoues for both wt and hp

quantile(Co_DF$F_Co, probs = seq(0, 1, 1/5))
quantile(Co_DF$TURNOVER, probs = seq(0, 1, 1/5))

Co_xgrid <-  seq(min(Co_DF$F_Co), max(Co_DF$F_Co), (max((Co_DF$F_Co-min(Co_DF$F_Co))/10)))
Co_ygrid <-  seq(min(Co_DF$TURNOVER), max(Co_DF$TURNOVER), (max((Co_DF$TURNOVER-min(Co_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Co.fit <-  expand.grid(F_Co = Co_xgrid, TURNOVER = Co_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Co_3d <-  predict(Co_loess, newdata = Co.fit)

# Abbreviated display of finCo matrix

Co_3d[1:4, 1:4]

# Transform data to long form
Co.melt <- melt(Co_3d, id.vars = c('F_Co', 'TURNOVER'), meCoure.vars = 'Co')

names(Co.melt) <- c('F_Co', 'TURNOVER', 'Co')

#Return data to numeric form

Co.melt$F_Co     <-  as.numeric(str_sub(Co.melt$F_Co, str_locate(Co.melt$F_Co, '=')[1,1] + 1))
Co.melt$TURNOVER    <-  as.numeric(str_sub(Co.melt$TURNOVER, str_locate(Co.melt$TURNOVER, '=')[1,1] + 1))


fig_Co <- plot_ly(Co.melt, x = ~F_Co, y = ~TURNOVER, z = ~Co, type = "contour",
                  width = 600, height = 500)


#### Cr ####

Cr_DF <- (COMPARTMENTS_AVG[,c("F_Cr_NH3","TURNOVER","Cr")])

Cr_DF <- Cr_DF[complete.cases(Cr_DF), ]

Cr_DF$TURNOVER <- asin(sqrt(Cr_DF$TURNOVER ))

quantile(Cr_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cr_DF$cat <- cut(Cr_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Cr_DF)<- c("F_Cr", "TURNOVER", "Cr", "cat")

Cr_TD <- glm(Cr ~ F_Cr, data = Cr_DF, family=gaussian)

no.int_Cr_TD  <- glm(Cr ~ (F_Cr + TURNOVER), data = Cr_DF, family=gaussian)

mix.int_Cr_TD <- glm(Cr ~ (F_Cr * TURNOVER), data = Cr_DF, family=gaussian)

models <- list(Cr_TD, mix.int_Cr_TD, no.int_Cr_TD)

aictab(c(models), modnames=c("Cr_TD","mix.int_Cr_TD","no.int_Cr_TD"), method=ML)

plot(allEffects(mix.int_Cr_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Cr_TD))

Anova(mix.int_Cr_TD)

with(summary(mix.int_Cr_TD), 1 - deviance/null.deviance)

### Cr Crl Plots ###

Cr_p=ggplot(data=Cr_DF,aes(F_Cr,Cr))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Cr, y=Cr, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Cr  (TD, mg/l)",    
       y= "Cr Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Cr_p

names(Cr_DF)

Cr_loess <- loess(Cr ~ F_Cr * TURNOVER, data = Cr_DF)

# Create a sequence of incrementCrly increCring (by 0.3 units) vCrues for both wt and hp

quantile(Cr_DF$F_Cr, probs = seq(0, 1, 1/5))
quantile(Cr_DF$TURNOVER, probs = seq(0, 1, 1/5))

Cr_xgrid <-  seq(min(Cr_DF$F_Cr), max(Cr_DF$F_Cr), (max((Cr_DF$F_Cr-min(Cr_DF$F_Cr))/10)))
Cr_ygrid <-  seq(min(Cr_DF$TURNOVER), max(Cr_DF$TURNOVER), (max((Cr_DF$TURNOVER-min(Cr_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Cr.fit <-  expand.grid(F_Cr = Cr_xgrid, TURNOVER = Cr_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Cr_3d <-  predict(Cr_loess, newdata = Cr.fit)

# Abbreviated display of finCr matrix

Cr_3d[1:4, 1:4]

# Transform data to long form
Cr.melt <- melt(Cr_3d, id.vars = c('F_Cr', 'TURNOVER'), meCrure.vars = 'Cr')

names(Cr.melt) <- c('F_Cr', 'TURNOVER', 'Cr')

#Return data to numeric form

Cr.melt$F_Cr     <-  as.numeric(str_sub(Cr.melt$F_Cr, str_locate(Cr.melt$F_Cr, '=')[1,1] + 1))
Cr.melt$TURNOVER    <-  as.numeric(str_sub(Cr.melt$TURNOVER, str_locate(Cr.melt$TURNOVER, '=')[1,1] + 1))


fig_Cr <- plot_ly(Cr.melt, x = ~F_Cr, y = ~TURNOVER, z = ~Cr, type = "contour",
                  width = 600, height = 500)

#### Mn ####

Mn_DF <- (COMPARTMENTS_AVG[,c("F_Mn","TURNOVER","Mnr")])

Mn_DF <- Mn_DF[complete.cases(Mn_DF), ]

Mn_DF$TURNOVER <- asin(sqrt(Mn_DF$TURNOVER ))

quantile(Mn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mn_DF$cat <- cut(Mn_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Mn_DF)<- c("F_Mn", "TURNOVER", "Mn", "cat")

Mn_TD <- glm(Mn ~ F_Mn, data = Mn_DF, family=gaussian)

no.int_Mn_TD  <- glm(Mn ~ (F_Mn + TURNOVER), data = Mn_DF, family=gaussian)

mix.int_Mn_TD <- glm(Mn ~ (F_Mn * TURNOVER), data = Mn_DF, family=gaussian)

models <- list(Mn_TD, mix.int_Mn_TD, no.int_Mn_TD)

aictab(c(models), modnames=c("Mn_TD","mix.int_Mn_TD","no.int_Mn_TD"), method=ML)

plot(allEffects(mix.int_Mn_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Mn_TD))

Anova(mix.int_Mn_TD)

with(summary(mix.int_Mn_TD), 1 - deviance/null.deviance)

### Mn Mnl Plots ###

Mn_p=ggplot(data=Mn_DF,aes(F_Mn,Mn))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Mn, y=Mn, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Mn  (TD, mg/l)",    
       y= "Mn Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Mn_p

names(Mn_DF)

Mn_loess <- loess(Mn ~ F_Mn * TURNOVER, data = Mn_DF)

# Create a sequence of incrementMnly increMning (by 0.3 units) vMnues for both wt and hp

quantile(Mn_DF$F_Mn, probs = seq(0, 1, 1/5))
quantile(Mn_DF$TURNOVER, probs = seq(0, 1, 1/5))

Mn_xgrid <-  seq(min(Mn_DF$F_Mn), max(Mn_DF$F_Mn), (max((Mn_DF$F_Mn-min(Mn_DF$F_Mn))/10)))
Mn_ygrid <-  seq(min(Mn_DF$TURNOVER), max(Mn_DF$TURNOVER), (max((Mn_DF$TURNOVER-min(Mn_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Mn.fit <-  expand.grid(F_Mn = Mn_xgrid, TURNOVER = Mn_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Mn_3d <-  predict(Mn_loess, newdata = Mn.fit)

# Abbreviated display of finMn matrix

Mn_3d[1:4, 1:4]

# Transform data to long form
Mn.melt <- melt(Mn_3d, id.vars = c('F_Mn', 'TURNOVER'), meMnure.vars = 'Mn')

names(Mn.melt) <- c('F_Mn', 'TURNOVER', 'Mn')

#Return data to numeric form

Mn.melt$F_Mn     <-  as.numeric(str_sub(Mn.melt$F_Mn, str_locate(Mn.melt$F_Mn, '=')[1,1] + 1))
Mn.melt$TURNOVER    <-  as.numeric(str_sub(Mn.melt$TURNOVER, str_locate(Mn.melt$TURNOVER, '=')[1,1] + 1))


fig_Mn <- plot_ly(Mn.melt, x = ~F_Mn, y = ~TURNOVER, z = ~Mn, type = "contour",
                  width = 600, height = 500)



#### Ni ####

Ni_DF <- (COMPARTMENTS_AVG[,c("F_Ni","TURNOVER","Ni")])

Ni_DF <- Ni_DF[complete.cases(Ni_DF), ]

Ni_DF$TURNOVER <- asin(sqrt(Ni_DF$TURNOVER ))

quantile(Ni_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ni_DF$cat <- cut(Ni_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Ni_DF)<- c("F_Ni", "TURNOVER", "Ni", "cat")

Ni_TD <- glm(Ni ~ F_Ni, data = Ni_DF, family=gaussian)

no.int_Ni_TD  <- glm(Ni ~ (F_Ni + TURNOVER), data = Ni_DF, family=gaussian)

mix.int_Ni_TD <- glm(Ni ~ (F_Ni * TURNOVER), data = Ni_DF, family=gaussian)

models <- list(Ni_TD, mix.int_Ni_TD, no.int_Ni_TD)

aictab(c(models), modnames=c("Ni_TD","mix.int_Ni_TD","no.int_Ni_TD"), method=ML)

plot(allEffects(mix.int_Ni_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Ni_TD))

Anova(mix.int_Ni_TD)

with(summary(mix.int_Ni_TD), 1 - deviance/null.deviance)

### Ni Nil Plots ###

Ni_p=ggplot(data=Ni_DF,aes(F_Ni,Ni))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_Ni, y=Ni, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "Ni  (TD, mg/l)",    
       y= "Ni Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

Ni_p

names(Ni_DF)

Ni_loess <- loess(Ni ~ F_Ni * TURNOVER, data = Ni_DF)

# Create a sequence of incrementNily increNiing (by 0.3 units) vNiues for both wt and hp

quantile(Ni_DF$F_Ni, probs = seq(0, 1, 1/5))
quantile(Ni_DF$TURNOVER, probs = seq(0, 1, 1/5))

Ni_xgrid <-  seq(min(Ni_DF$F_Ni), max(Ni_DF$F_Ni), (max((Ni_DF$F_Ni-min(Ni_DF$F_Ni))/10)))
Ni_ygrid <-  seq(min(Ni_DF$TURNOVER), max(Ni_DF$TURNOVER), (max((Ni_DF$TURNOVER-min(Ni_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

Ni.fit <-  expand.grid(F_Ni = Ni_xgrid, TURNOVER = Ni_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

Ni_3d <-  predict(Ni_loess, newdata = Ni.fit)

# Abbreviated display of finNi matrix

Ni_3d[1:4, 1:4]

# Transform data to long form
Ni.melt <- melt(Ni_3d, id.vars = c('F_Ni', 'TURNOVER'), meNiure.vars = 'Ni')

names(Ni.melt) <- c('F_Ni', 'TURNOVER', 'Ni')

#Return data to numeric form

Ni.melt$F_Ni     <-  as.numeric(str_sub(Ni.melt$F_Ni, str_locate(Ni.melt$F_Ni, '=')[1,1] + 1))
Ni.melt$TURNOVER    <-  as.numeric(str_sub(Ni.melt$TURNOVER, str_locate(Ni.melt$TURNOVER, '=')[1,1] + 1))


fig_Ni <- plot_ly(Ni.melt, x = ~F_Ni, y = ~TURNOVER, z = ~Ni, type = "contour",
                  width = 600, height = 500)


#### S ####

S_DF <- (COMPARTMENTS_AVG[,c("F_S","TURNOVER","S.2")])

S_DF <- S_DF[complete.cases(S_DF), ]

S_DF$TURNOVER <- asin(sqrt(S_DF$TURNOVER ))

quantile(S_DF$TURNOVER, probs = seq(0, 1, 1/5))

S_DF$cat <- cut(S_DF$TURNOVER, 
                breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(S_DF)<- c("F_S", "TURNOVER", "S", "cat")

S_TD <- glm(S ~ F_S, data = S_DF, family=gaussian)

no.int_S_TD  <- glm(S ~ (F_S + TURNOVER), data = S_DF, family=gaussian)

mix.int_S_TD <- glm(S ~ (F_S * TURNOVER), data = S_DF, family=gaussian)

models <- list(S_TD, mix.int_S_TD, no.int_S_TD)

aictab(c(models), modnames=c("S_TD","mix.int_S_TD","no.int_S_TD"), method=ML)

plot(allEffects(mix.int_S_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_S_TD))

Anova(mix.int_S_TD)

with(summary(mix.int_S_TD), 1 - deviance/null.deviance)

### S Sl Plots ###

S_p=ggplot(data=S_DF,aes(F_S,S))+
  geom_point(aes(color=factor(cat)),size=4)+
  geom_smooth(aes(x=F_S, y=S, group = factor(cat), color=factor(cat), fill=factor(cat)), 
              method="glm", 
              formula = y ~ x, 
              method.args=list(family="gaussian")) +
  labs(title="Title",
       x= "S  (TD, mg/l)",    
       y= "S Content (mg/g)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

S_p

names(S_DF)

S_loess <- loess(S ~ F_S * TURNOVER, data = S_DF)

# Create a sequence of incrementSly increSing (by 0.3 units) vSues for both wt and hp

quantile(S_DF$F_S, probs = seq(0, 1, 1/5))
quantile(S_DF$TURNOVER, probs = seq(0, 1, 1/5))

S_xgrid <-  seq(min(S_DF$F_S), max(S_DF$F_S), (max((S_DF$F_S-min(S_DF$F_S))/10)))
S_ygrid <-  seq(min(S_DF$TURNOVER), max(S_DF$TURNOVER), (max((S_DF$TURNOVER-min(S_DF$TURNOVER))/10)))

# Generate a dataframe with every possible combination of wt and hp

S.fit <-  expand.grid(F_S = S_xgrid, TURNOVER = S_ygrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp

S_3d <-  predict(S_loess, newdata = S.fit)

# Abbreviated display of finS matrix

S_3d[1:4, 1:4]

# Transform data to long form
S.melt <- melt(S_3d, id.vars = c('F_S', 'TURNOVER'), meSure.vars = 'S')

names(S.melt) <- c('F_S', 'TURNOVER', 'S')

#Return data to numeric form

S.melt$F_S     <-  as.numeric(str_sub(S.melt$F_S, str_locate(S.melt$F_S, '=')[1,1] + 1))
S.melt$TURNOVER    <-  as.numeric(str_sub(S.melt$TURNOVER, str_locate(S.melt$TURNOVER, '=')[1,1] + 1))


fig_S <- plot_ly(S.melt, x = ~F_S, y = ~TURNOVER, z = ~S, type = "contour",
                 width = 600, height = 500)
