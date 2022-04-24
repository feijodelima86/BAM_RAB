### This code produces a spreadsheet with the calculated burdens for every sample. Did not handle QAQC data.

library(openxlsx)
library(tidyverse)
library(readxl)
library(lubridate) 


### Function to match data frame classes

matchColClasses <- function(df1, df2) {
  
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  
  return(df2)
}

### Set a working file path for OES spreadsheets. MANI(pulated) because I changed file names to have consistent dates and 
### the OES term so the list function could find them.

path1<-"0_Data/Filemergetest/OES_MANI/"

### This command finds all spreadsheets with OES in their name within the folder specified above.

OES.LIST<-c(list.files(path= path1, full.names=TRUE,pattern=c("OES")))

### This command creates a list of objects with all those files.

OES<-lapply(OES.LIST,function(filename){
  print(paste("Merging",filename,sep.=""))
  data.frame(read.csv(filename, fileEncoding="UTF-8-BOM"))
})

### This removes the first line that contains only "mg/l" from all data frames within the list.

OES <- lapply(OES,function(x){
  x[-c(1,0),]
})

### This is a trick to remove numbers from all column headings. The make.unique function then adds unique identifiers to elements that are repeated.
### E.G. Iron has 3 types of detection (1 axial 2 radial), Those become Fe, Fer and Fer1. 

OES <- lapply(OES, function(x) 
{names(x) <- make.unique(gsub("[^[:alpha:]]", "", names(x)));x})

### This is where all spreadsheets are merged into one. 

OES.CALC<-OES%>%reduce(full_join)

### First output with samples and all QAQC samples associated.

write.csv(OES.CALC, "2_Incremental/Calc/OES_CALC.csv")

### Removing all QAQC samples from the dataset. 

OES.SAMPLES<-data.frame(OES.CALC)

OES.SAMPLES<-OES.SAMPLES[!is.na(as.numeric(OES.SAMPLES$SampleID)), ]

### Second output with samples without QAQC info. 

write.csv(OES.SAMPLES, "2_Incremental/Calc/OES_CALC_SAMPLES.csv")


###Bringing in sample specific biomass data

path2<-"0_Data/Filemergetest/BIO_MANI/"

### This command finds all spreadsheets with BIOMASS in their name within the folder specified above.

BIO.LIST<-c(list.files(path= path2, full.names=TRUE,pattern=c("BIOMASS")))

### This command creates a list of objects with all those files.

BIO<-lapply(BIO.LIST,function(filename){
  print(paste("Merging",filename,sep.=""))
  data.frame(read.xlsx(filename))
})

### This command reduces the list to the columns of interest 

BIO.CALC<-lapply(BIO, "[", c("DATE","SITE","SAMPLE","COMPARTMENT","TOTAL.DRY.WEIGHT..g."))

### Making classes of variables be tha same for merger

class(BIO.CALC[[1]][,3]) = "character"
class(BIO.CALC[[2]][,3]) = "character"
class(BIO.CALC[[3]][,3]) = "character"

### Merging all datasets

BIO.CALC<-BIO.CALC%>%reduce(full_join)

### First output with all samples containing date, site, sample, compartment and Dry weights stacked.

write.csv(BIO.CALC, "2_Incremental/Calc/BIO_CALC_SAMPLES.csv")

###Bringing in sample log file to make translations between sample specific biomass and OES files

SLOG<-read.xlsx("0_Data/Filemergetest/SLOG_MANI/20220216_SAMPLE_LOG_RFL_DTW_2022_3_3_MASS.xlsx")

SLOG<-SLOG[!is.na(as.numeric(SLOG$FIELD.REP)),]

names(SLOG)

#Renaming columns so reduce works

SLOG.SEL<-SLOG[,c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR","FIELD.REP","OES_DATE","TUBE_NUMBER","ICP_NUMBER")]

names(BIO.CALC)<-c("SAMPLING_DATE","SITE","FIELD.REP","SAMPLE_DESCRIPTOR","TOTAL.DRY.WEIGHT..g.")

SLOG.BIO <- list(SLOG.SEL, BIO.CALC)

#merge all data frames in list

SLOG.BIO<-SLOG.BIO %>% reduce(full_join)

write.csv(SLOG.BIO, "2_Incremental/Calc/SLOG_BIO_SAMPLES.csv")

#merge all data frames


matchColClasses <- function(df1, df2) {
  
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  
  return(df2)
}



# note added ugly formats below

SLOG.BIO$OES_DATE<-as.Date(as.numeric(SLOG.BIO$OES_DATE), origin = "1899-12-30")

SLOG.BIO$OES_DATE <- ymd(SLOG.BIO$OES_DATE) 
OES.SAMPLES$Date <- mdy(OES.SAMPLES$Date) 

names(OES.SAMPLES)[names(OES.SAMPLES) == 'Date'] <- 'OES_DATE'
names(OES.SAMPLES)[names(OES.SAMPLES) == 'SampleID'] <- 'ICP_NUMBER'


matchColClasses(SLOG.BIO, OES.SAMPLES)

print(lapply(OES.SAMPLES, class))
print(lapply(SLOG.BIO, class))

BB.SAMPLES <- list(OES.SAMPLES, SLOG.BIO)

BB.SAMPLES<-BB.SAMPLES %>% reduce(full_join)

write.csv(BB.SAMPLES, "2_Incremental/Calc/BB_SAMPLES.csv")
