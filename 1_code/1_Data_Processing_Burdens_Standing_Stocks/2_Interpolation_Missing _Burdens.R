library("dplyr")
library("readr")
library("tidyr")
library("mgcv")
library("mgcViz")
require(gam)

#if na on column x then interpolate burdens.

alldata <- read.csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-factor(alldata$SAMPLING_DATE,levels = c("6/22/2021", "7/7/2021", "7/20/2021", "8/3/2021", "8/17/2021","9/9/2021", "9/22/2021", "10/13/2021"), ordered=TRUE)

names(alldata)

dfEPIL <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIL"),]

dfEPIP <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "EPIP"),]

dfFILA <- alldata[which(alldata$SAMPLE_DESCRIPTOR == "FILA"),]

dfALL <- alldata %>% filter(
  SAMPLE_DESCRIPTOR== "EPIL" | SAMPLE_DESCRIPTOR == "EPIP"| SAMPLE_DESCRIPTOR == "FILA",
)

n<-c(13,19,22,25,35,49) ## (As:13 Cd:19,Cu:22 Fe_r234:25, Pb:35, Zn:49)
#c(13,19,22,25,35,49)+39


se <- function(x, ...) sqrt(var(x, ...)/length(x))

ALL.SUM<-data.frame(aggregate(as.numeric(dfALL[,n[1]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean))
ALL.SUM[,5]<- aggregate(as.numeric(dfALL[,n[1]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]
ALL.SUM[,6]<- aggregate(as.numeric(dfALL[,n[2]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean)[,4]
ALL.SUM[,7]<- aggregate(as.numeric(dfALL[,n[2]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]
ALL.SUM[,8]<- aggregate(as.numeric(dfALL[,n[3]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean)[,4]
ALL.SUM[,9]<- aggregate(as.numeric(dfALL[,n[3]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]
ALL.SUM[,10]<-aggregate(as.numeric(dfALL[,n[4]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean)[,4]
ALL.SUM[,11]<-aggregate(as.numeric(dfALL[,n[4]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]
ALL.SUM[,12]<-aggregate(as.numeric(dfALL[,n[5]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean)[,4]
ALL.SUM[,13]<-aggregate(as.numeric(dfALL[,n[5]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]
ALL.SUM[,14]<-aggregate(as.numeric(dfALL[,n[6]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, mean)[,4]
ALL.SUM[,15]<-aggregate(as.numeric(dfALL[,n[6]]) ~ SAMPLING_DATE+SITE+SAMPLE_DESCRIPTOR, dfALL, se)[,4]

names(ALL.SUM)<-c("DATE","SITE","SAMPLE_DESCRIPTOR","MEAN.AS","STDER.AS","MEAN.CD","STDER.CD","MEAN.CU","STDER.CU","MEAN.FE","STDER.FE","MEAN.PB","STDER.PB","MEAN.ZN","STDER.ZN")

head(ALL.SUM)

ALL.SUM<-ALL.SUM %>%
  mutate(SITENUM = case_when(
    endsWith(SITE, "WS") ~ "1",
    endsWith(SITE, "DL") ~ "6",
    endsWith(SITE, "GR") ~ "8",
    endsWith(SITE, "GC") ~ "9",
    endsWith(SITE, "BG") ~ "10",
    endsWith(SITE, "BN") ~ "11"
  ))


ALL.SUM<-ALL.SUM %>% arrange(SAMPLE_DESCRIPTOR,SITENUM,DATE)
ALL.SUM<-ALL.SUM[,-16]


ncol(ALL.SUM)
ALL.SUM %>% subset(SITE == "BN")

EPIL.SUM <- ALL.SUM[which(ALL.SUM$SAMPLE_DESCRIPTOR == "EPIL"),]

EPIP.SUM <- ALL.SUM[which(ALL.SUM$SAMPLE_DESCRIPTOR == "EPIP"),]

FILA.SUM <- ALL.SUM[which(ALL.SUM$SAMPLE_DESCRIPTOR == "FILA"),]

#write.csv(ALL.SUM, paste0("2_incremental/",gsub("-", "", Sys.Date()),"_BURDENS_SS_SUMMARY.csv"))

names(alldata)

dfEPIL$SITE <- as.factor(dfEPIL$SITE) 

color_easy <- c("red", "blue", "yellow", "green", "orange", "purple")[dfEPIL$SITE]

as.numeric(dfEPIL$SAMPLING_DATE)

plot(dfEPIL$SAMPLING_DATE,dfEPIL$As, type="p", col=color_easy)

GAMDATA=data.frame(as.numeric(dfEPIL$SAMPLING_DATE),dfEPIL$SITE,dfEPIL$As)

names(GAMDATA)<-c("DATE","SITE","As")

GAM <- mgcv::gam(As ~ SITE + s(DATE , by = SITE , k=7),
                                       data = GAMDATA)
plot(GAM, pages=1)

summary(GAM)

plot(GAM,residuals=TRUE,pch=19,cex=.5)

pd<-data.frame(SITE = c("WS","DL","GC","BG","BN","GR"),DATE=as.numeric(as.Date("2021/07/20")))
pv <- predict(GAM,pd)

pd[,3]<-pv

pd[,2]<-as.Date(pd[,2], origin='1970-01-01')

plot(pv)

b <- getViz(GAM)

plot(b)


##### FORLOOP TEST #####

## make up some data
library(mgcv)
set.seed(0) 
dat <- gamSim(1,n=200,scale=2)
set.seed(1)
dat2 <- gamSim(1,n=200,scale=2)
names(dat2)[1:5] <- c("y1", paste0("x", 4:7))
d <- cbind(dat[, 1:5], dat2[, 1:5])
d_resp <- d[ c("y", "y1")]
d_pred <- d[, !(colnames(d) %in% c("y", "y1"))]

## create a "matrix" list of dimensions i x j
results_m <- vector("list", length=ncol(d_resp)*ncol(d_pred))
dim(results_m) <- c(ncol(d_resp), ncol(d_pred))

for(i in 1:ncol(d_resp)){
  for(j in 1:ncol(d_pred)){
    results_m[i, j][[1]] <- gamm(d_resp[, i] ~ s(d_pred[, j]))
  }
}

# flatten the "matrix" list
results_l <- do.call("list", results_m)

