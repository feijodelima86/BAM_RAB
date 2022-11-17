library("readr")
library("corrplot")
library("vegan")

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]
COMPARTMENTS <- COMPARTMENTS[-217,]

nrow(COMPARTMENTS)
names(COMPARTMENTS)

####Least VIF####

LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,6),c(17,18,20,31,35,36,40,41,44)+39)]))
LABELS<-LABELS[order(LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
PCA.DF<-as.data.frame(LABELS[,c(4:ncol(LABELS))])
names(LABELS)
names(PCA.DF)<- c("Be","Ca","Co","Mo","Pb","S","Se","Si","Sn")

cor(PCA.DF)
corrplot(cor(PCA.DF), method = "circle", lwd=2) 


####Big 5+Fe+Se####

LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,6),c(13,19,22,25,35,40,49)+39)]))
LABELS<-LABELS[order(LABELS$SAMPLE_DESCRIPTOR, decreasing = F), ]
PCA.DF<-as.data.frame(LABELS[,c(4:ncol(LABELS))])
names(PCA.DF)<- c("As","Cd","Cu","Fe","Pb","Se","Zn")

cor(PCA.DF)
corrplot(cor(PCA.DF), method = "circle", lwd=2) 


###PCA###

PCA.DF<-log(PCA.DF)
my.data <- PCA.DF
my.data <- as.matrix(scale(my.data, center = TRUE, scale = TRUE))
nrow(my.data)
my.prc <- prcomp(my.data)

plot(my.prc)

sd <- my.prc$sdev
correlations <- t(t(my.prc$rotation)*sd)

summary(my.prc)

#####PLOT####

colors <- c(colors()[89], "gold", "chartreuse")
colors <- colors[as.numeric(factor(na.omit(LABELS$SAMPLE_DESCRIPTOR)))]

dev.new()

plot	(my.prc$x,
      pch=c(23,24,25)[as.numeric(factor(LABELS$SAMPLE_DESCRIPTOR))],
      xlim=c(-max(abs(range(my.prc$x[,1]))),max(abs(range(my.prc$x[,1])))),	
      ylim=c(-max(abs(range(my.prc$x[,2]))),max(abs(range(my.prc$x[,2])))), 
      #      asp = 1,
      cex.lab=1.5, cex.axis=1.5, cex.main=1.4, cex.sub=1.4, 
      col="black",
      bg=colors,
      cex=1.5,
      lwd=2
)
#	col="black")

tab <- matrix(c(my.prc$x[,1], my.prc$x[,2]), ncol=2)

panel.first= {
  ordiellipse(tab,LABELS$SAMPLE_DESCRIPTOR,conf=0.95, draw ="polygon", col=c(adjustcolor(colors()[89], alpha.f=0.25),adjustcolor("gold", alpha.f=0.25),adjustcolor("chartreuse", alpha.f=0.25)), border=c(colors()[89],"gold","chartreuse"))
}


legend(x = "topright", legend = levels(factor(LABELS$SAMPLE_DESCRIPTOR)), 
       pch=c(23,24,25), 
       cex = 1.5, 
       box.lwd=3,
       col = "black", 
       pt.bg=c(colors()[89], "gold", "chartreuse"), pt.lwd=3)

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

text(	x=correlations[,1]*1.05,
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

