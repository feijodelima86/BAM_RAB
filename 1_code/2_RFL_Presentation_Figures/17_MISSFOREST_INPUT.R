library("readr")
library("missMDA")
library("missForest")

alldata <- read_csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

alldata <- transform(alldata, SAMPLING_DATE = as.numeric(alldata$SAMPLING_DATE))

alldata$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[alldata$SITE]

names(alldata)

mydata <- data.frame(alldata$SAMPLING_DATE,alldata$SITE_DISTANCE,alldata$OM.AREA.g.m2,alldata$Cd,alldata$Cu,alldata$Fe,alldata$Pb,alldata$Zn)

mydata$SAMPLE_DESCRIPTOR <- as.factor(alldata$SAMPLE_DESCRIPTOR)

mydata <- mydata[complete.cases(mydata), ]

z_scores <- as.data.frame(sapply(mydata[,c(1:ncol(mydata)-1)], function(mydata) (abs(mydata-mean(mydata))/sd(mydata))))

mydata <- mydata[!rowSums(z_scores>2.0), ]

names(mydata) <- c("SAMPLING_DATE","SITE_DISTANCE","OM.AREA.g.m2","Cd","Cu","Fe","Pb","Zn","SAMPLE_DESCRIPTOR")

mydata_miss<-mydata

mydata_miss$Fe <- as.numeric((lapply(mydata_miss$Fe, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.70, 0.30), size = length(cc), replace = TRUE)])))

#mydata_miss <- prodNA(mydata, noNA = 0.2)

mydata_miss_imp <- missForest(mydata_miss, xtrue = mydata, verbose = TRUE)

summary(lm(mydata_miss_imp$ximp[,"Fe"] ~ mydata[,"Fe"]))

write.csv(mydata_miss_imp$ximp, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"completedata.csv"))



