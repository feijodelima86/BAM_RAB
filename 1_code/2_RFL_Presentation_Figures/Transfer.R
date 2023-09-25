#### Al ####

Al_DF <- (COMPARTMENTS_AVG[,c("F_Al","TURNOVER","Alr")])

Al_DF <- Al_DF[complete.cases(Al_DF), ]

Al_DF$TURNOVER <- asin(sqrt(Al_DF$TURNOVER ))

quantile(Al_DF$TURNOVER, probs = seq(0, 1, 1/5))

Al_DF$cat <- cut(Al_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(Al_DF)<- c("F_Al", "TURNOVER", "Al", "cat")

Al_TD <- glm(Al ~ F_Al, data = Al_DF, family=gaussian)

no.int_Al_TD  <- glm(Al ~ (F_Al + TURNOVER), data = Al_DF, family=gaussian)

mix.int_Al_TD <- glm(Al ~ (F_Al * TURNOVER), data = Al_DF, family=gaussian)

models <- list(Al_TD, mix.int_Al_TD, no.int_Al_TD)

aictab(c(models), modnames=c("Al_TD","mix.int_Al_TD","no.int_Al_TD"), method=ML)

plot(allEffects(mix.int_Al_TD), lines = list(multiline = T), confint = list(style = "auto"))

plot(allEffects(mix.int_Al_TD))

Anova(mix.int_Al_TD)

with(summary(mix.int_Al_TD), 1 - deviance/null.deviance)

### Al all Plots ###

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

As_DF <- (COMPARTMENTS_AVG[,c("F_As_Oshift","TURNOVER","As")])

As_DF <- As_DF[complete.cases(As_DF), ]

As_DF$TURNOVER <- asin(sqrt(As_DF$TURNOVER ))

quantile(As_DF$TURNOVER, probs = seq(0, 1, 1/5))

As_DF$cat <- cut(As_DF$TURNOVER, 
                 breaks=c(-Inf, 0.2186771    , 0.4699721    , 0.5016606   , 0.5487757, 0.5981569   , Inf), 
                 labels=c("NA", "0.21 - 0.46","0.46 - 0.50","0.50 - 0.54","0.54 - 0.59","0.59 - 0.76"))

names(As_DF)<- c("F_As", "TURNOVER", "As", "cat")

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

As_xgrid <-  seq(min(As_DF$F_As), max(As_DF$F_As), (max((As_DF$F_As-min(As_DF$F_As))/10)))
As_ygrid <-  seq(min(As_DF$TURNOVER), max(As_DF$TURNOVER), (max((As_DF$TURNOVER-min(As_DF$TURNOVER))/10)))

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
                  width = 600, height = 500)

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

Cd_xgrid <-  seq(min(Cd_DF$F_Cd), max(Cd_DF$F_Cd), (max((Cd_DF$F_Cd-min(Cd_DF$F_Cd))/10)))
Cd_ygrid <-  seq(min(Cd_DF$TURNOVER), max(Cd_DF$TURNOVER), (max((Cd_DF$TURNOVER-min(Cd_DF$TURNOVER))/10)))

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
                  width = 600, height = 500)
