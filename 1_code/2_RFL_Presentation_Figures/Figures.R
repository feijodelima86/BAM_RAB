library(plyr)
library(tidyverse)
library(readr)


alldata <- read.csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE)

names(alldata)

# Selecting variables

n<-7
n2<-49 #(As:13 Cd:19,Cu:22 Fe_r234:25, Pb:35, Zn:49)

#Names of axes

Xlabel=expression(bold(paste("Standing Crop"~"(g/m"^"2"~")")))
Ylabel=expression(bold("Zn Burden (mg/g)"))

#Subsetting by compartment

x <- alldata[,n] 
y <- alldata[,n2] 

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

#Graphical paramenters

BW1<-5     #Smoothing for x axis dispersion
BW2<-.0004      #Smoothing for y axis dispersion (As:0.009 Cd:0.0004,Cu:0.035,Fe:1.5, Pb:0.006, Zn:2)
ALPHA<-0.6  #Transparency of dispersion curves

dev.new()

#jpeg(file='3_products/Manuscript/Figures/Cd_Burdens300_2.jpg',height = 5, width = 5, units = 'in', res = 300)

#Matrix layout for plots

layout(matrix(c(1,3,3,1,3,3,0,2,2),ncol=3), c(3.5,1), c(1,3))

layout.show(n=3)
#Plot1

par(mar=c(0.5, 4.5, 0.5, 0.5))

plot(density(dfEPIP[,n], na.rm = TRUE, bw = BW1), 
     xlim = range(alldata[,n],na.rm=TRUE), 
     col=NA,
     lwd=3,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     main=NA,
     ylim=c(0.001,max(density(dfEPIP[,n], na.rm = TRUE, bw = BW1)$y)),
     
)     


aty <- seq(0.0, max(c(density(dfEPIP[,n], na.rm = TRUE, bw = BW1)$y,
                      density(dfEPIL[,n], na.rm = TRUE, bw = BW1)$y,
                      density(dfFILA[,n], na.rm = TRUE, bw = BW1)$y)), length.out=5)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

polygon(density(dfEPIP[,n], na.rm = TRUE, bw = BW1), 
        xlim = range(alldata[,n],na.rm=TRUE), 
        col=adjustcolor("gold",alpha.f=ALPHA),
        #     border="gold",
        lwd=2)  
polygon(density(dfEPIL[,n], na.rm = TRUE, bw = BW1), col=adjustcolor(col=colors()[89],alpha.f=ALPHA),
        #        border=colors()[89],
        lwd=2)          
polygon(density(dfFILA[,n], na.rm = TRUE, bw = BW1), col = adjustcolor("chartreuse3",alpha.f=ALPHA),
        #        border="chartreuse3",
        lwd=2) 

box(lwd=2)

# Plot 2 


dnEPIP<-density(dfEPIP[,n2], na.rm = TRUE, bw = BW2)
x1 <- dnEPIP$y
y1 <- dnEPIP$x

dnEPIL<-density(dfEPIL[,n2], na.rm = TRUE, bw = BW2)
x2 <- dnEPIL$y
y2 <- dnEPIL$x


dnFILA<-density(dfFILA[,n2], na.rm = TRUE, bw = BW2)
x3 <- dnFILA$y
y3 <- dnFILA$x


par(mar=c(4.5, 0.5, 0.5, 0.5))

plot(density(dfEPIP[,n2], na.rm = TRUE, bw = BW1)$y, 
     xlim = c(0.000,max(c(x1,x2,x3))), 
     col=NA,
     lwd="",
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     main=NA,
     ylim=c(min(alldata[,n2], na.rm = TRUE),max(alldata[,n2], na.rm = TRUE)),
     )     



polygon(dnEPIP$y,
        dnEPIP$x, 
#        xlim = range(x1), 
        col=adjustcolor("gold",alpha.f=ALPHA),
#        border="gold",
        lwd=2)  


polygon(dnEPIL$y,
        dnEPIL$x, 
#        xlim = range(x2), 
        col=adjustcolor(colors()[89],alpha.f=ALPHA),
        #     border="gold",
        lwd=2)  


polygon(dnFILA$y,
        dnFILA$x, 
#       xlim = range(x3), 
        col=adjustcolor("chartreuse3",alpha.f=ALPHA),
        #     border="gold",
        lwd=2)  

box(lwd=2)

aty1 <- seq(0.0, max(c(x1,x2,x3)), length.out=5)

axis(side = 1, at = aty1, labels=format(aty1, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

#Plot 3

dev.new()

par(mar=c(4.5, 4.5, 0.5, 0.5))

plot(alldata[,n], alldata[,n2], 
     ylim = range(alldata[,n2],na.rm=TRUE),
     xlim = range(alldata[,n],na.rm=TRUE),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)

aty <- seq(0, max(alldata[,n2], na.rm=TRUE), length.out=5)
axis(side = 1, at = seq(0, 200, length.out=5), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5)
title(ylab=Ylabel, line=2.7, cex.lab=1.5)


points(dfEPIL[,7], dfEPIL[,n2],pch=23, cex=1.5,col="black", bg=colors()[89],lwd=3)

points(dfEPIP[,7], dfEPIP[,n2],pch=23, cex=1.5,col="black", bg="gold",lwd=3)

points(dfFILA[,7], dfFILA[,n2],pch=23, cex=1.5,col="black", bg="chartreuse3",lwd=3)

box(lwd=2)

#dev.off()

