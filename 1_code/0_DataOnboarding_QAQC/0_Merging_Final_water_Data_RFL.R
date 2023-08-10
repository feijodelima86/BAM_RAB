library("readr")
library("missMDA")
library("missForest")
library("stringr")
library("tidyverse")

alldata   <- read.csv("2_incremental/20220420_STANDING_CROP.csv")
waterdata <- read.csv("2_incremental/Final_water_data_Redux.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

alldata$SITE_DISTANCE <- as.numeric(str_replace_all(alldata$SITE, c("WS"="0","DL"="44.9","GR"="64.82","GC"="89.22","BG"="144.32","BN"="167.82")))

waterdata$date <- as.Date(waterdata$date, format = "%m/%d/%Y")

waterdata$SITE <- factor(waterdata$site, labels = c("WS","DL","GR","GC","BG","BN"))

colnames(waterdata)[c(2:3)]<-c("SITENUM","SAMPLING_DATE")

names(waterdata)

WTR.DATA<-waterdata[,c(3,28,4:27)]

names(alldata)

ALG.DATA<-alldata[,c(3,4,89,5,6,7:88)]


matchColClasses <- function(df1, df2) {
  
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  
  return(df2)
}

matchColClasses(ALG.DATA, WTR.DATA)

# Creating list

WATER.ALGAE<-list(ALG.DATA, WTR.DATA)

#Joining data frames

WATER.ALGAE<-WATER.ALGAE %>% reduce(full_join)

names(WATER.ALGAE)

write.csv(WATER.ALGAE, "2_incremental/wateralgae_Final_1.csv")


