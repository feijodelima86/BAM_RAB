library(tidyverse)
library(readr)
library(plyr)
library(dplyr)

alldata <- data.frame(read_csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_Dont_Use.csv"))

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

se <- function(x, ...) sqrt(var(x, ...)/length(x))

#alldata[,c(11-87)]<-alldata[,c(13-87)]

#alldata<-replace(alldata[c(13:88)], alldata[c(13:87)] < 0, 0)

# Variable date to factor

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

#Selected site 

ssite<-"GC"

# Selected element

names(alldata)

n1<-87

mult=1

EPIL.SUM<-data.frame(aggregate(as.numeric(dfEPIL[,n1]) ~ SAMPLING_DATE+SITE, dfEPIL, mean))
EPIL.SUM[,4]<-aggregate(as.numeric(dfEPIL[,n1]) ~ SAMPLING_DATE+SITE, dfEPIL, se)[,3]
EPIL.SUM <- subset(EPIL.SUM, SITE == ssite)
names(EPIL.SUM)<-c("DATE","SITE","MEAN.EPIL","STDER.EPIL")
EPIL.SUM[,5]<--2
EPIL.SUM[,c(3,4)]<-EPIL.SUM[,c(3,4)]*mult 

FILA.SUM<-data.frame(aggregate(as.numeric(dfFILA[,n1]) ~ SAMPLING_DATE+SITE, dfFILA, mean))
FILA.SUM[,4]<-aggregate(as.numeric(dfFILA[,n1]) ~ SAMPLING_DATE+SITE, dfFILA, se)[,3]
FILA.SUM <- subset(FILA.SUM, SITE == ssite)
names(FILA.SUM)<-c("DATE","SITE","MEAN.FILA","STDER.FILA")
FILA.SUM[,5]<--2
FILA.SUM[,c(3,4)]<-FILA.SUM[,c(3,4)]*mult 

EPIP.SUM<-data.frame(aggregate(as.numeric(dfEPIP[,n1]) ~ SAMPLING_DATE+SITE, dfEPIP, mean))
EPIP.SUM[,4]<-aggregate(as.numeric(dfEPIP[,n1]) ~ SAMPLING_DATE+SITE, dfEPIP, se)[,3]
EPIP.SUM <- subset(EPIP.SUM, SITE == ssite)
names(EPIP.SUM)<-c("DATE","SITE","MEAN.EPIP","STDER.EPIP")
EPIP.SUM[,5]<--2
EPIP.SUM[,c(3,4)]<-EPIP.SUM[,c(3,4)]*mult 


ALL.SUM<-list(EPIP.SUM,FILA.SUM,EPIL.SUM)

ALL.SUM<-ALL.SUM %>% reduce(full_join)

ALL.SUM[is.na(ALL.SUM)] <- 0

L=nrow(ALL.SUM)

#Column Width

j=2

#each 'polygon' is inside a list with xx and yy coordinates

dat1 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.SUM[x,1]+j, ALL.SUM[x,1]-j, ALL.SUM[x,1]-j, ALL.SUM[x,1]+j),
              yy=c(ALL.SUM[x,5], ALL.SUM[x,5], ALL.SUM[x,8], ALL.SUM[x,8]))
  return(res)
})

dat2 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.SUM[x,1]+j, ALL.SUM[x,1]-j, ALL.SUM[x,1]-j, ALL.SUM[x,1]+j),
              yy=c(ALL.SUM[x,5], ALL.SUM[x,5], ALL.SUM[x,6]+ALL.SUM[x,8], ALL.SUM[x,6]+ALL.SUM[x,8]))
  return(res)
})

dat3 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.SUM[x,1]+j, ALL.SUM[x,1]-j, ALL.SUM[x,1]-j, ALL.SUM[x,1]+j),
              yy=c(ALL.SUM[x,5], ALL.SUM[x,5], ALL.SUM[x,3]+ALL.SUM[x,6]+ALL.SUM[x,8], ALL.SUM[x,3]+ALL.SUM[x,6]+ALL.SUM[x,8]))
  return(res)
})

ardat1 <- lapply(1:L,function(x){
  res <- list(x0=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y0=c(ALL.SUM[x,8]),
              x1=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y1=c(ALL.SUM[x,8]+ALL.SUM[x,9]))
  return(res)
})

ardat2 <- lapply(1:L,function(x){
  res <- list(x0=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y0=c(ALL.SUM[x,6]+ ALL.SUM[x,8]),
              x1=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y1=c(ALL.SUM[x,6]+ALL.SUM[x,8]+ALL.SUM[x,7]))
  return(res)
})
ardat3 <- lapply(1:L,function(x){
  res <- list(x0=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y0=c(ALL.SUM[x,6]+ALL.SUM[x,8]+ALL.SUM[x,3]),
              x1=c(ALL.SUM[x,1], ALL.SUM[x,1]),
              y1=c(ALL.SUM[x,6]+ALL.SUM[x,8]+ALL.SUM[x,3]+ALL.SUM[x,4]))
  return(res)
})


#create empty plot

#dev.new(width=10, height=7, noRStudioGD = TRUE)
plot(ALL.SUM[,1], ALL.SUM[,5],
     type='n',
     ylim=c(20, plyr::round_any(max(ALL.SUM[3]+ALL.SUM[6]+ALL.SUM[8]+ALL.SUM[4]), 100, f = ceiling)),
     xlim= c(as.Date("2021-06-20"),as.Date("2021-10-30")),
     ylab=NA, 
     yaxt = "n",
     cex=1.5,
     xaxt = "n",
     xlab = NA,
     cex.lab=1, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
)


axis(side = 2, las=1, font.axis=1, cex.axis=1.5)


plyr::round_any(range(ALL.SUM[3]+ALL.SUM[6]+ALL.SUM[8]), 10, f = ceiling)

at = seq(0, plyr::round_any(max(ALL.SUM[3]+ALL.SUM[6]+ALL.SUM[8]), 175, f = ceiling), length.out=5)

axis.Date(1, at=seq(as.Date("2021-06-20"), as.Date("2021-10-30"), "1 week"), las=1, font.axis=2, cex.axis=1.5, format="%d-%b")

#add polygons

for (i in 1:L) polygon(unlist(dat3[[i]]$xx),unlist(dat3[[i]]$yy),col="gold", border="black")
for (i in 1:L) polygon(unlist(dat2[[i]]$xx),unlist(dat2[[i]]$yy),col="chartreuse3", border="black")
for (i in 1:L) polygon(unlist(dat1[[i]]$xx),unlist(dat1[[i]]$yy),col=colors()[89],border="black")

warnings()

for (i in 1:L) arrows(unlist(ardat3[[i]]$x0),unlist(ardat3[[i]]$y0),unlist(ardat3[[i]]$x1),unlist(ardat3[[i]]$y1), length=0.075, angle=90, code=2, lwd=2,col="black")
for (i in 1:L) arrows(unlist(ardat2[[i]]$x0),unlist(ardat2[[i]]$y0),unlist(ardat2[[i]]$x1),unlist(ardat2[[i]]$y1), length=0.075, angle=90, code=2, lwd=2,col="black")
for (i in 1:L) arrows(unlist(ardat1[[i]]$x0),unlist(ardat1[[i]]$y0),unlist(ardat1[[i]]$x1),unlist(ardat1[[i]]$y1), length=0.075, angle=90, code=2, lwd=2,col="black")


box(lwd=4)
options("warn"=1)


