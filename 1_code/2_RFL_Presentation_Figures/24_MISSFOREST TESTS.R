library("readr")
library("missMDA")
library("missForest")

#Shotgun

alldata <- read_csv("2_incremental/20220420_STANDING_CROP.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

names(alldata)

mydata <- data.frame(alldata[ ,c(3:49)])

mydata <- mydata[complete.cases(mydata), ]

mydata <- filter(mydata, FIELD.REP %in% c(1:5)) 

mydata$SITE <- as.factor(mydata$SITE)

mydata$SAMPLING_DATE <- as.numeric(mydata$SAMPLING_DATE)

mydata <- mydata[which(mydata$SAMPLE_DESCRIPTOR == "EPIL"| mydata$SAMPLE_DESCRIPTOR == "EPIP" | mydata$SAMPLE_DESCRIPTOR == "FILA"),]

mydata$SITE_DISTANCE <- as.numeric(str_replace_all(mydata$SITE, c("WS"="0","DL"="44.9","GR"="64.82","GC"="89.22","BG"="144.32","BN"="167.82")))

mydata$SAMPLE_DESCRIPTOR <- as.factor(mydata$SAMPLE_DESCRIPTOR)

mydata$FIELD.REP<- as.factor(mydata$FIELD.REP)

sum(is.na(mydata))/prod(dim(mydata))

names(mydata)

z_scores <- as.data.frame(sapply(mydata[,c(5:47)], function(mydata) (abs(mydata-mean(mydata))/sd(mydata))))

mydata <- mydata[!rowSums(z_scores>2.5), ] 

mydata_05<-data.frame(mydata)

mydata_05$As <- replace(mydata_05$As , sample(seq_along(mydata_05$As), 35), NA)

#mydata_05 <- data.frame(sapply(mydata_05, \(x) replace(x, sample(seq_along(x), 3), NA)))

sum(is.na(mydata_05$As))/length(mydata_05$As)

mydata_miss_imp_05 <- missForest(mydata_05, xtrue = mydata, verbose = TRUE, replace = T)

mydata_miss_imp_05$ximp$As

setdiff(mydata_miss_imp_05$ximp$As, mydata$As)

setdiff(mydata$As, mydata_miss_imp_05$ximp$As)

LM_05<-lm(setdiff(mydata_miss_imp_05$ximp$As, mydata$As) ~ setdiff(mydata$As, mydata_miss_imp_05$ximp$As))

summary(LM_05)

plot(setdiff(mydata_miss_imp_05$ximp$As, mydata$As), setdiff(mydata$As, mydata_miss_imp_05$ximp$As))

sum(is.na(mydata_05$As))/length(mydata_05$As)

#write.csv(mydata_miss_imp$ximp, paste0("2_incremental/PCA_INPUT_TEST/",gsub("-", "", Sys.Date()),"completedata_SHOTGUN_ALL.csv"))

