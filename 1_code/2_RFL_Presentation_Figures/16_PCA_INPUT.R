library("readr")
library("missMDA")
library("missForest")

alldata <- read_csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

alldata <- transform(alldata, ndate = as.numeric(alldata$SAMPLING_DATE))

alldata$SITE_DISTANCE<-c(WS=0,DL=44.9,GR=64.82,GC=89.22,BG=144.32,BN=167.82)[alldata$SITE]

names(alldata)

EPIP <- ifelse(alldata$SAMPLE_DESCRIPTOR == 'EPIP', 1, 0)

EPIL <- ifelse(alldata$SAMPLE_DESCRIPTOR == 'EPIL', 1, 0)

FILA <- ifelse(alldata$SAMPLE_DESCRIPTOR == 'FILA', 1, 0)

mydata <- data.frame(alldata$ndate,alldata$SITE_DISTANCE,EPIP,EPIL,FILA,alldata$OM.AREA.g.m2,alldata$As)

names(mydata) <- c("ndate","SITE_DISTANCE","EPIP","EPIL","FILA","OM.AREA.g.m2","As")

mydata <- mydata[complete.cases(mydata), ]

z_scores <- as.data.frame(sapply(mydata, function(mydata) (abs(mydata-mean(mydata))/sd(mydata))))

mydata <- mydata[!rowSums(z_scores>2.5), ]

write.csv(mydata, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"completedata_AS.csv"))

mydata_05<-mydata

mydata_05$As <- as.numeric((lapply(mydata_05$As, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.80, 0.20), size = length(cc), replace = TRUE)])))

mydata_05_input <- imputePCA(mydata_05, ncp = 5)

summary(lm(mydata_05_input$completeObs[,"As"] ~ mydata_05_input$fittedX[,7]))

plot(mydata_05_input$completeObs[,"As"],mydata_05_input$fittedX[,7])

range(mydata_05_input$completeObs[,"As"])
range(mydata_05_input$fittedX[,7])

write.csv(mydata_05, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"mydata_05_As.csv"))

write.csv(mydata_05_input, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"mydata_05_As_input.csv"))

