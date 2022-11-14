alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

nrow(COMPARTMENTS)

#Big 5

#PCA.DF<-as.data.frame(na.omit(COMPARTMENTS[,c(13,19,22,25,35,40,49)]))

#PCA.DF<-log(PCA.DF)

#names(PCA.DF)<- c("As","Cd","Cu","Fe","Pb","Se","Zn")

#cor(PCA.DF)

#corrplot(cor(PCA.DF), method = "circle", lwd=2) 

#Least VIF

PCA.DF<-as.data.frame(na.omit(COMPARTMENTS[,c(17,18,20,31,35,36,40,41,44)]))

nrow(PCA.DF)

PCA.DF<-log(PCA.DF)

nrow(PCA.DF)

names(PCA.DF)<- c("Be","Ca","Co","Mo","Pb","S","Se","Si","Sn")

cor(PCA.DF)

corrplot(cor(PCA.DF), method = "circle", lwd=2) 

names(COMPARTMENTS)

colors <- c("blue", "darkgreen", "red")
colors <- colors[as.numeric(factor(na.omit(COMPARTMENTS$SAMPLE_DESCRIPTOR)))]

mydata <- PCA.DF

my.scaled.data <- as.matrix(scale(mydata, center = TRUE, scale = TRUE))

view(my.scaled.data)

ncol(my.scaled.data)

my.prc <- prcomp(na.omit(my.scaled.data), center = TRUE, scale = TRUE)

wa<-data.frame(colors,my.prc$x[,1],my.prc$x[,2])

plot(my.prc)

sd <- my.prc$sdev
correlations <- t(t(my.prc$rotation)*sd)

summary(my.prc)

#dev.new()

plot	(my.prc$x,
      pch=c(0,1,2)[as.numeric(factor(COMPARTMENTS$SAMPLE_DESCRIPTOR))],
      xlim=c(-max(abs(range(my.prc$x[,1]))),max(abs(range(my.prc$x[,1])))),	
      ylim=c(-max(abs(range(my.prc$x[,2]))),max(abs(range(my.prc$x[,2])))), 
      cex.lab=1.2, cex.axis=1.1, cex.main=1.4, cex.sub=1.4, 
      col=colors
        )
#	col="black")


par(new=TRUE)

plot	(correlations[1,],
      correlations[2,],
      col="White",
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

text(	x=correlations[,1]*1.05,
      y=correlations[,2]*1.05,
      labels = names(PCA.DF), 
      cex=0.7)

rownames(my.prc$rotation)

#legend('topright', legend = levels(wa$Stream), col = 1:4, cex = 1, pch = 19)

coordnames<-data.frame(correlations[,1]*1.1,correlations[,2]*1.1)

coordnames

axis(3)
mtext("y2",side=4,line=2, col="Blue")
axis(4)
mtext("y2",side=4,line=2, col="Blue")

abline(h=c(0,0))
abline(v=c(0,0))

