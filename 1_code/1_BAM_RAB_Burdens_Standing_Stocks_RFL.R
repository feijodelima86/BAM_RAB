library(tidyverse)
library(readr)
library(lubridate) 

BB.RAW <- read.csv("2_incremental/20220312_MDL_SAMPLES.csv")

#Creating dataset with relevant columns and correcting for dilution (50ml) to obtain mg/sample

BB.PRUNE <- cbind(BB.RAW[,c(42,43,44,45,53)], BB.RAW[,c(3:41)]*50/1000)

#Dividing all metal masses by sample mass in g. Burdens are in mg/g.

BB.FINAL <- cbind(BB.PRUNE[,c(1,2,3,4)],(BB.PRUNE[,c(6:44)]/BB.PRUNE[,5]))

write.csv(BB.FINAL, paste0("2_incremental/",gsub("-", "", Sys.Date()),"_BB_CURRENT.csv"))

#Bringing in LTREB biomass spreadsheets. 

INTEGRATED_BIOMASS_2021_03_10 <- read.csv("0_data/external/FLIPPED_BIOMASS_2021_03_10.csv", fileEncoding="UTF-8-BOM")

#Renaming columns so reduce works

names(INTEGRATED_BIOMASS_2021_03_10)[names(INTEGRATED_BIOMASS_2021_03_10) == 'DATE'] <- 'SAMPLING_DATE'
names(INTEGRATED_BIOMASS_2021_03_10)[names(INTEGRATED_BIOMASS_2021_03_10) == 'SAMPLE'] <- 'FIELD.REP'

#Adjusting date/time formats

INTEGRATED_BIOMASS_2021_03_10$SAMPLING_DATE <-mdy(INTEGRATED_BIOMASS_2021_03_10$SAMPLING_DATE) 
BB.FINAL$SAMPLING_DATE <- dmy(BB.FINAL$SAMPLING_DATE)

#Function to make column classes the same

matchColClasses <- function(df1, df2) {
  
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  
  return(df2)
}


matchColClasses(INTEGRATED_BIOMASS_2021_03_10, BB.FINAL)

# Creating list

STANDING.CROP<-list(INTEGRATED_BIOMASS_2021_03_10, BB.FINAL)

#Joining data frames

STANDING.CROP<-STANDING.CROP %>% reduce(full_join)

#Adding columns with calculated standing crops. Final values are in mg/m2 (mg/g *g/m2)

names(STANDING.CROP)

STANDING.CROP[,c(49:87)]<-STANDING.CROP[,c(10:48)]*STANDING.CROP[,6]

names(STANDING.CROP[,c(49:87)])<-paste(names(STANDING.CROP[,c(49:87)]), "SC")

#CSV output.

write.csv(STANDING.CROP, paste0("2_incremental/",gsub("-", "", Sys.Date()),"_STANDING_CROP.csv"))

