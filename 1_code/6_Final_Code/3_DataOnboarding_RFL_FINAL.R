library(readr)
library(tidyverse)


algaedata <- data.frame(read.csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_2.csv"))

waterdata <- data.frame(read.csv("2_incremental/230801_Water_Metals_Size_diff.csv"))

turnoverdata <- read_csv("2_incremental/TURNOVER_REDUX.csv")

### prepping algae data ###

names(algaedata)

algaedata$SAMPLING_DATE<-as.Date(algaedata$SAMPLING_DATE, format = "%m/%d/%Y")

algaedata<-algaedata[,-c(1,2)]

### prepping water data ###

names(waterdata)

colnames(waterdata)[c(2:3)]<-c("SITE","SAMPLING_DATE")

water_calc<-waterdata

factor(waterdata$SAMPLING_DATE)

factor(waterdata$SITE)

water_calc<-data.frame(waterdata)

water_calc$SITE <- factor(water_calc$SITE, labels = c("WS","DL","GR","GC","BG","BN"))

waterdata<-data.frame(water_calc)

colnames(water_calc)

factor(waterdata$SAMPLING_DATE)

waterdata$SAMPLING_DATE<-as.Date(factor(waterdata$SAMPLING_DATE, labels = c("2021-06-22", "2021-07-07", "2021-07-20", "2021-08-03", "2021-08-17", "2021-08-30", "2021-09-09", "2021-09-22", "2021-10-13")))

### prepping turnover data ###

names(turnoverdata)

turnoverdata$SAMPLING_DATE<-as.Date(turnoverdata$SAMPLING_DATE, format = "%m/%d/%Y")

### full join reduce ###

add.comp<-list(algaedata,waterdata,turnoverdata)

alldata<-add.comp %>% reduce(full_join)

write.csv(alldata, "2_incremental/TURNOVER_Full_Dataset.csv")

### aggredating dataset by means###

COMPARTMENTS_AVG <-aggregate(x = alldata[,colnames(alldata) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
                             by = list(alldata$SAMPLING_DATE,alldata$SITE,alldata$SAMPLE_DESCRIPTOR),
                             FUN = mean,
                             na.rm = T)
COMPARTMENTS_AVG

names(COMPARTMENTS_AVG)

colnames(COMPARTMENTS_AVG)[c(1:3)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

write.csv(COMPARTMENTS_AVG, "2_incremental/TURNOVER_Full_Dataset_AVG.csv")

