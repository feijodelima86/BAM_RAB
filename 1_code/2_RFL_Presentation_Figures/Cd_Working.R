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
