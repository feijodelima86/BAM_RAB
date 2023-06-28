library("readr")
library("missMDA")
library("missForest")

#Shotgun

alldata <- read_csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

mydata <- data.frame(alldata[ ,c(3:49)])

mydata_05$As <- as.numeric((lapply(mydata_05$As, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.80, 0.20), size = length(cc), replace = TRUE)])))

mydata$SITE <- as.factor(mydata$SITE)

mydata <- mydata[which(mydata$SAMPLE_DESCRIPTOR == "EPIL"| mydata$SAMPLE_DESCRIPTOR == "EPIP" | mydata$SAMPLE_DESCRIPTOR == "FILA"),]

mydata$SAMPLE_DESCRIPTOR <- as.factor(mydata$SAMPLE_DESCRIPTOR)

mydata <- filter(mydata, FIELD.REP %in% c(1:5)) 

mydata$FIELD.REP<- as.factor(mydata$FIELD.REP)

mydata$SAMPLING_DATE <- as.factor(mydata$SAMPLING_DATE)

names(mydata)

lapply(mydata, class)

sum(is.na(mydata))/prod(dim(mydata))

mydata_miss_imp <- missForest(mydata, xtrue = mydata, verbose = TRUE)

write.csv(mydata_miss_imp$ximp, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"completedata_SHOTGUN_ALL.csv"))

