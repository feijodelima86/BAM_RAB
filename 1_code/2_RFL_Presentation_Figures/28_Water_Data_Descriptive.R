library("readr")
library("dplyr")
library("ggplot2"); theme_set(theme_bw() +
                                theme(axis.line = element_line(color='black'),
                                      plot.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank()))
library("scales")
library("splines")
library("effects")
library("mgcv")

COMPARTMENTS <- read.csv("2_incremental/water5.csv")

colnames(COMPARTMENTS)[c(2:3)]<-c("SAMPLING_DATE","SITE")

COMPARTMENTS$SITE <- factor(COMPARTMENTS$SITE, labels = c("WS","DL","GR","GC","BG","BN"))


#Selected site 

ssite<-"WS"

# Selected element

names(COMPARTMENTS)

# P= 4 Cd = 12 Fe= 6 As = 9

n1<-12

mult=1

#Column Width

j=2

# Summary frames for each size fraction

TD.SUM<-data.frame(aggregate(as.numeric(COMPARTMENTS[,n1+20]) ~ SAMPLING_DATE+SITE, COMPARTMENTS, mean))
CL.SUM<-data.frame(aggregate(as.numeric(COMPARTMENTS[,n1+10]) ~ SAMPLING_DATE+SITE, COMPARTMENTS, mean))
SS.SUM<-data.frame(aggregate(as.numeric(COMPARTMENTS[,n1]) ~ SAMPLING_DATE+SITE, COMPARTMENTS, mean))

ALL.wSUM<-list(TD.SUM,CL.SUM,SS.SUM)

ALL.wSUM<-ALL.wSUM %>% purrr::reduce(full_join)

ALL.wSUM <- subset(ALL.wSUM, SITE == ssite)
names(ALL.wSUM)<-c("DATE","SITE","TD.SUM","CL.SUM","SS.SUM")


#each 'polygon' is inside a list with xx and yy coordinates

ALL.wSUM[is.na(ALL.wSUM)] <- 0

ALL.wSUM[,6]<--2

L=nrow(ALL.wSUM)

ALL.wSUM$DATE<-as.Date(ALL.wSUM$DATE, format = "%m/%d/%Y")

dat1 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.wSUM[x,1]+j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]+j),
              yy=c(ALL.wSUM[x,6], ALL.wSUM[x,6], ALL.wSUM[x,3], ALL.wSUM[x,3]))
  return(res)
})

dat2 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.wSUM[x,1]+j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]+j),
              yy=c(ALL.wSUM[x,6], ALL.wSUM[x,6], ALL.wSUM[x,3]+ALL.wSUM[x,4], ALL.wSUM[x,3]+ALL.wSUM[x,4]))
  return(res)
})

dat3 <- lapply(1:L,function(x){
  res <- list(xx=c(ALL.wSUM[x,1]+j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]-j, ALL.wSUM[x,1]+j),
              yy=c(ALL.wSUM[x,6], ALL.wSUM[x,6], ALL.wSUM[x,3]+ALL.wSUM[x,4]+ALL.wSUM[x,5], ALL.wSUM[x,3]+ALL.wSUM[x,4]+ALL.wSUM[x,5]))
  return(res)
})

#create empty plot

#dev.new(width=12, height=6, noRStudioGD = TRUE)
par(mar=c(4,6,4,4))
plot(ALL.wSUM[,1], ALL.wSUM[,5],
     type='n',
     ylim=c(0, plyr::round_any(max(ALL.wSUM[,3]+ALL.wSUM[,4]+ALL.wSUM[,5]), 0.06, f = ceiling)),
     xlim= c(as.Date("2021-06-20"),as.Date("2021-10-30")),
     ylab=NA, 
     yaxt = "n",
     cex=1.5,
     xaxt = "n",
     xlab = NA,
     cex.lab=1, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
)


axis(side = 2, las=1, font.axis=2, cex.axis=2)

at = seq(0, plyr::round_any(max(ALL.wSUM[,3]+ALL.wSUM[,4]+ALL.wSUM[,5]), 175, f = ceiling), length.out=5)

axis.Date(1, at=seq(as.Date("2021-06-20"), as.Date("2021-10-30"), "2 weeks"), las=1, font.axis=2, cex.axis=2, format="%d-%b")

for (i in 1:L) polygon(unlist(dat3[[i]]$xx),unlist(dat3[[i]]$yy),col="tomato4", border="black")
for (i in 1:L) polygon(unlist(dat2[[i]]$xx),unlist(dat2[[i]]$yy),col="plum", border="black")
for (i in 1:L) polygon(unlist(dat1[[i]]$xx),unlist(dat1[[i]]$yy),col="turquoise3",border="black")

box(lwd=4)

