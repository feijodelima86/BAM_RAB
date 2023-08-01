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

### To Do: 

### 1: Perform tests on PCA to see if different compartments are truly different.

library("corrplot")
library("vegan")
library("pairwiseAdonis")


alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

nrow(COMPARTMENTS)
names(COMPARTMENTS)

####Big 5+Fe+Se+Mo####

LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(3,4,6,13,19,22,25,31,35,40,49)]))
LABELS<-LABELS[order(LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
PCA.DF<-as.data.frame(LABELS[,c(4:ncol(LABELS))])
names(PCA.DF)<- c("As","Cd","Cu","Fe","Mo","Pb","Se","Zn")

#dev.new()
cor(PCA.DF)
corrplot(cor(PCA.DF), method = "circle", lwd=2) 


###PCA###

PCA.DF<-log(PCA.DF)
my.data <- PCA.DF
my.data <- as.matrix(scale(my.data, center = TRUE, scale = TRUE))

my.prc <- prcomp(na.omit(my.data))

groups <- factor(LABELS$SAMPLE_DESCRIPTOR)

dis <- vegdist(my.data, method="euclidean")

pairwise.adonis2(dis ~ groups, perm=999, method="euclidean")

plot(my.prc)

pca_scores<-scores(my.prc)

sd <- my.prc$sdev
correlations <- t(t(my.prc$rotation)*sd)

summary(my.prc)

#####PLOT####


colors <- c(colors()[89], "gold", "chartreuse")
colors <- colors[as.numeric(factor(na.omit(LABELS$SAMPLE_DESCRIPTOR)))]

#dev.new()

my.prc$x[,1]<-my.prc$x[,1]*-1
correlations[,1]<-correlations[,1]*-1

plot	(my.prc$x[,1], my.prc$x[,2],
      pch=1,
      xlim=c(-max(abs(range(my.prc$x[,1]))),max(abs(range(my.prc$x[,1])))),	
      ylim=c(-max(abs(range(my.prc$x[,2]))),max(abs(range(my.prc$x[,2])))), 
      cex.lab=1.5, cex.axis=1.5, cex.main=1.4, cex.sub=1.4, 
      col="white",
      xlab="PC 1",
      ylab="PC 2"
)
#	col="black")

tab <- matrix(c(my.prc$x[,1], my.prc$x[,2]), ncol=2)

legend(x = "topleft", legend = levels(factor(LABELS$SAMPLE_DESCRIPTOR)), 
       pch=c(23,24,25), 
       cex = 1, 
       box.lwd=3,
       col = "black", 
       pt.bg=c(colors()[89], "gold", "chartreuse"), pt.lwd=3)


panel.first= {
  ordiellipse(tab,LABELS$SAMPLE_DESCRIPTOR,conf=0.95, draw ="polygon", col=c(adjustcolor(colors()[89], alpha.f=0.25),adjustcolor("gold", alpha.f=0.25),adjustcolor("chartreuse", alpha.f=0.25)), border=c(colors()[89],"gold","chartreuse"), lwd=3)
}

points(my.prc$x[,1],my.prc$x[,2], pch=c(23,24,25)[as.numeric(factor(LABELS$SAMPLE_DESCRIPTOR))],col="black", bg=colors,    cex=1.5,     lwd=2)

par(new=TRUE)

plot	(correlations[1,],
      correlations[2,],
      col="NA",
      xlim=c(-1,1),
      ylim=c(-1,1),
      xaxt="n",
      yaxt="n",
      xlab="",
      ylab="")

arrows(0, 0, 
       x1 = correlations[,1], 
       y1 = correlations[,2], 
       length = 0.1, angle = 30, lwd =2.0,
       code = 2, col = par("fg"), 
       lty = par("lty"))

text(	x=correlations[,1]*1.075,
      y=correlations[,2]*1.05,
      labels = names(PCA.DF), 
      cex=1.5,
      font=2)

rownames(my.prc$rotation)

box(lwd=3)

#legend('topright', legend = levels(wa$Stream), col = 1:4, cex = 1, pch = 19)

coordnames<-data.frame(correlations[,1]*1.1,correlations[,2]*1.1)

coordnames

axis(3, cex.axis=1.5)
mtext("y2",side=4,line=2, col="Blue")
axis(4, cex.axis=1.5)
mtext("y2",side=4,line=2, col="Blue")

abline(h=c(0,0))
abline(v=c(0,0))

my.prc


### 2: Perform GLMs with totals instead of subtracted numbers for different size fractions to assess if different filtration methods wyeld different results.

### 3: Add Alice's true numbers to assess metals turnover for different compartments. 

##### Study #####

## Euclidean distances between samples

dis <- vegdist(my.data, method="euclidean")

## First 16 sites grazed, remaining 8 sites ungrazed

adonis2(groups ~ dis, strata=dat$field, perm=999)

## Calculate multivariate dispersions
mod <- betadisper(dis, groups, sqrt.dist=T, bias.adjust=T, type = "centroid")
mod

## Perform test

anova(mod)

names(mod)

mod$centroids

## Permutation test for F

permutest(mod, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences

(mod.HSD <- TukeyHSD(mod))

#dev.new()

plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.95) # 90% data ellipse

## can also specify which axes to plot, ordering respected
plot(mod, axes = c(2,1), seg.col = "forestgreen", seg.lty = "dashed")

## Draw a boxplot of the distances to centroid for each group

boxplot(mod)

## `scores` and `eigenvals` also work
scrs <- scores(mod)
str(scrs)
head(scores(mod, 1:4, display = "sites"))

# group centroids/medians 
scores(mod, 1:4, display = "centroids")

# eigenvalues from the underlying principal coordinates analysis
eigenvals(mod) 

## try out bias correction; compare with mod3
(mod3B <- betadisper(dis, groups, type = "median", bias.adjust=TRUE))
anova(mod3B)
permutest(mod3B, permutations = 99)
mod.HSD3 <- TukeyHSD(mod3B)

plot(mod3B, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod3B, ellipse = TRUE, hull = FALSE, conf = 0.95) # 90% data ellipse



## should always work for a single group
group <- factor(rep("grazed", NROW(varespec)))
(tmp <- betadisper(dis, group, type = "median"))
(tmp <- betadisper(dis, group, type = "centroid"))

## simulate missing values in 'd' and 'group'
## using spatial medians
groups[c(2,20)] <- NA
dis[c(2, 20)] <- NA
mod2 <- betadisper(dis, groups) ## messages
mod2
permutest(mod2, permutations = 99)
anova(mod2)
plot(mod2)
boxplot(mod2)
plot(TukeyHSD(mod2))

## Using group centroids
mod3 <- betadisper(dis, groups, type = "centroid")
mod3
permutest(mod3, permutations = 99)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3))
