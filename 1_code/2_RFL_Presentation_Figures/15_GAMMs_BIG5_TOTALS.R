library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw())
library("scales")
library("AICcmodavg")
library("gam")
library("visreg")
library("mgcv")

#alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP.csv"))
alldata <- read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

#COMPARTMENTS <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"| alldata$SAMPLE_DESCRIPTOR == "EPIP" | alldata$SAMPLE_DESCRIPTOR == "FILA"),]

names(COMPARTMENTS)<-c("ROWNUM",names(COMPARTMENTS[c(2:ncol(COMPARTMENTS))]))

COMPARTMENTS$ROWNUM==1

COMPARTMENTS <- COMPARTMENTS[COMPARTMENTS$ROWNUM != 1, ]

COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]


#COMPARTMENTS <- COMPARTMENTS[which(COMPARTMENTS$SITE != "WS"),]

#COMPARTMENTS$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[COMPARTMENTS$SITE]

#COMPARTMENTS<-COMPARTMENTS[-1,]#<- removed outlier*

####Big 5+Fe+Se####

#+39

BIG5.LABELS<-as.data.frame(na.omit(COMPARTMENTS[,c(c(3,4,5,6,7,89),c(13,19,22,25,35,40,49)+39)]))
BIG5.LABELS[,c(5:12)] %>% mutate_if(is.character, as.numeric)
BIG5.LABELS[,3]<-as.character(BIG5.LABELS[,3])
BIG5.LABELS$SITE<-factor(BIG5.LABELS$SITE)
BIG5.LABELS<-distinct(BIG5.LABELS, SAMPLING_DATE, SITE,FIELD.REP, .keep_all= TRUE)
names(BIG5.LABELS)[c(7:ncol(BIG5.LABELS))]<- c("As","Cd","Cu","Fe","Pb","Se","Zn")
BIG5.GLM <- BIG5.LABELS[,-4] %>% group_by(SAMPLING_DATE, SITE, FIELD.REP) %>% summarise_all(sum)

names(BIG5.GLM)

K=7

mydata<-na.omit(BIG5.GLM[,c(1,5,6)])




#mydata[,3]<-log10(mydata[,3])

names(mydata)

xlab="Distance (Km)"
ylab="Julian Day"

#mydata[,1] %>% mutate_if(is.character, as.numeric)

x <- mydata[ ,1]
y <- mydata[ ,2]
z <- mydata[ ,3]

#z <- log10(mydata[ ,3])

mydata <- transform(mydata, ndate = as.numeric(x),
                    nyear  = as.numeric(format(x, '%Y')),
                    nmonth = as.numeric(format(x, '%m')),
                    doy    = as.numeric(format(x, '%j')))

names(mydata)<-c("x","y","z")

summary(lm4)

anova(lm4)


K=7

?gam

mod_gam <- gam(z ~ s(ndate, k=8)+s(y, k=6), data=mydata)


dat <- data.frame(mydata[,3])
rbPal <- colorRampPalette(c('blue','red'))
dat$Col <- rbPal(10)[as.numeric(cut(mydata[,3],breaks = 10))]

p <- quote({
  axis(1, at = mx, labels = lx)
  axis(2, at = my, labels = ly)
  with(mydata, points(x, y, pch=19, col="black", cex=0.8))
})

#dev.new(width=12, height=10)



plot<-visreg2d(mod_gam, xvar='ndate', yvar='y', scale='response',
               #			plot.type="persp",
               xlab=xlab, ylab=ylab, 
               main=NULL, 
               plot.axes=p, 
               #               zlim=c(min(z),max(z)),
               #               zlim=c(0,30000),
               color.palette=colorRampPalette(c("white","red", "darkred")),
               #               color.palette=colorRampPalette(c("white","gold", "peru")),
               #               color.palette=colorRampPalette(c("white","blue", "navyblue")),
               #               color.palette=colorRampPalette(c("lightblue","blue", "navyblue")),
               font.lab = 2,
               font.axis = 2,
               #			cex.axis = 0.8,
)

summary(mod_gam)








