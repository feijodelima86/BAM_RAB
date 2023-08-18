library(readr)

algaedata <- data.frame(read.csv("2_incremental/20220420_STANDING_CROP_Rafa_Interpolation_2.csv"))

waterdata <- data.frame(read.csv("2_incremental/230731_Water_Metals_long_diff.csv"))

turnoverdata <- read_csv("2_incremental/TURNOVER_REDUX.csv")

### prepping algae data ###

names(algaedata)

algaedata$SAMPLING_DATE<-as.Date(algaedata$SAMPLING_DATE, format = "%m/%d/%Y")

algaedata<-algaedata[,-c(1,2)]

### prepping water data ###

names(waterdata)

water_seperate <- waterdata %>%
  # Separate the date, size, and site using str_extract 
  mutate(
    date = str_extract(Sample.Id, "\\d+_\\d+_\\d+"),             # Extract date in the format: 6_21_21
    filtersize = str_extract(Sample.Id, "(W|UF|F)"),                   # Extract size (W, UF, or F)
    site = str_extract(Sample.Id, "(?<=_(W|UF|F)_)(\\d+)")       # Extract the site using lookbehind regex
  ) %>%
  # Convert date to the actual date type using lubridate package
  mutate(date = mdy(date))

water_calc <- water_seperate %>%
  pivot_wider(id_cols = c(site, date),values_from = Concentration, names_from = c(filtersize, Element), values_fn = mean)


colnames(water_calc)[c(1:2)]<-c("SITE","SAMPLING_DATE")

water_calc$SITE <- factor(water_calc$SITE, labels = c("WS","BG","BN","DL","GR","GC"))

water_calc$SAMPLING_DATE<-as.Date(water_calc$SAMPLING_DATE, format = "%m/%d/%Y")

water_calc<-data.frame(water_calc)

colnames(water_calc)

waterdata<-water_calc[,c(1,2,6,13,19,20,23,29,31,33,35,39,46,52,53,56,62,64,66,68,72,79,85,86,89,95,97,99,101)]

factor(waterdata$SAMPLING_DATE)

waterdata$SAMPLING_DATE<-as.Date(factor(waterdata$SAMPLING_DATE, labels = c("2021-06-22", "2021-07-07", "2021-07-20", "2021-08-03", "2021-08-17", "2021-08-30", "2021-09-09", "2021-09-22", "2021-10-13")))

### prepping turnover data ###

names(turnoverdata)

turnoverdata$SAMPLING_DATE<-as.Date(turnoverdata$SAMPLING_DATE, format = "%m/%d/%Y")

### full join reduce ###

add.comp<-list(algaedata,waterdata,turnoverdata)

alldata<-add.comp %>% reduce(full_join)

### aggredating dataset ###

COMPARTMENTS_AVG <-aggregate(x = alldata[,colnames(alldata) != c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")],          
                             by = list(alldata$SAMPLING_DATE,alldata$SITE,alldata$SAMPLE_DESCRIPTOR),
                             FUN = mean,
                             na.rm = T)
COMPARTMENTS_AVG

names(COMPARTMENTS_AVG)

colnames(COMPARTMENTS_AVG)[c(1:3)]<-c("SAMPLING_DATE","SITE","SAMPLE_DESCRIPTOR")

write.csv(COMPARTMENTS_AVG, "2_incremental/TURNOVER_Full_Dataset.csv")

