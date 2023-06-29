library("readr")
library("missMDA")
library("missForest")
library("stringr")
library(tidyverse)

alldata   <- read.csv("2_incremental/20220420_STANDING_CROP.csv")
waterdata <- read.csv("2_incremental/water4.csv")

alldata$SAMPLING_DATE<-as.Date(alldata$SAMPLING_DATE, format = "%m/%d/%Y")

alldata$SITE_DISTANCE <- as.numeric(str_replace_all(alldata$SITE, c("WS"="0","DL"="44.9","GR"="64.82","GC"="89.22","BG"="144.32","BN"="167.82")))

waterdata$Group.1 <- as.Date(waterdata$Group.1, format = "%m/%d/%Y")

waterdata$SITE <- factor(waterdata$Group.2, labels = c("WS","DL","GR","GC","BG","BN"))

colnames(waterdata)[c(2:3)]<-c("SAMPLING_DATE","SITENUM")

names(waterdata)

names(alldata)

ALG.DATA<-alldata[,c(3,4,89,5,6,7:88)]

WTR.DATA<-waterdata[,c(2,34,4:33)]

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

write.csv(WATER.ALGAE, "2_incremental/wateralgae13.csv")


