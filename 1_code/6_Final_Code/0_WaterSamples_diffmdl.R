## Script for processing water samples during the BAMRAB time period 

#install and load packages 
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(lubridate)

options(scipen=100)

# check sesion info for reporting
sessionInfo()

#writing micron
x <- "\u00b5"
x
#load data and clean ----

water1 <- read.csv(file.path("./0_data/icp_ms_internal","230516_BAMRAB1REPRO1.csv"), header=T, na.strings = c(""))
water1.2 <- read.csv(file.path("./0_data/icp_ms_internal","230516_BAMRAB1REPROend.csv"), header=T, na.strings = c(""))
water2 <- read.csv(file.path("./0_data/icp_ms_internal","230608_BAMRABrun2.1.csv"), header=T, na.strings = c(""))
water2.2 <- read.csv(file.path("./0_data/icp_ms_internal","230609_BAMRAB2.2.csv"), header=T, na.strings = c(""))
water3 <- read.csv(file.path("./0_data/icp_ms_internal","230616_BAMRAB3.1.csv"), header=T, na.strings = c(""))



# load sample info
#sampleinfo <- read.csv(file.path("./0_data","SampleInfo.csv"), header=T, na.strings = c(""))

#remove rows that we do not want to use (rerun needed or a split reprocess)

water1 <- water1[c(1:177),] 
water1.2 <- water1.2[c(119:151),] 
water2 <- water2[c(15:132),] 
water2.2 <- water2.2[c(152:293),] 

#remove rows that failed anlaysis
#water1= water1 %>% filter(QC.Status == "Passed")
#this dataset has failed values that we want

#select sample id column and columns that are elements

water1 <- water1 %>% select(c(Sample.Id, contains("ppb.")))

water1.2 <- water1.2 %>% select(c(Sample.Id, contains("ppb.")))

water2 <- water2 %>% select(c(Sample.Id, contains("ppb.")))

water2.2 <- water2.2 %>% select(c(Sample.Id, contains("ppb.")))

water3 <- water3 %>% select(c(Sample.Id, contains("ppb.")))

#remove any columns(elements) that are not in all the datasets

water2names <- colnames(water2)

#remove old int.stds.

water1 <- water1 %>% select(any_of(water2names))
water1.2 <- water1.2 %>% select(any_of(water2names))


#bind data frames together
water1 <- rbind(water1, water1.2 )
water2 <- rbind(water2, water2.2)

water2 <-water2[,c(1:3,5:36)] #remove Al std mode

water <- rbind(water1, water2)
water <- water[,c(1:21,23:35)] #remove Rh oshift

water3 <- water3[c(1:105, 108:235),] # T-239s that were labeled incorectly as Dblanks

names <- colnames(water)
print(names)
#rename columns for water1
water_names <- water %>% rename( c(Ca = starts_with("Ca"),Mg = starts_with("Mg"),Al = starts_with("Al"), S = starts_with("S."),P = starts_with("P."),K = starts_with("K"),Na = starts_with("Na"), V = starts_with("V"),Cr_KED =starts_with("Cr..52"),Cr_NH3 = starts_with("Cr..53"),Fe_NH3 = starts_with("Fe..54"),Mn = starts_with("Mn"),Fe_KED = starts_with("Fe..56"),Co = starts_with("Co"),Ni = starts_with("Ni"),Cu_KED = starts_with("Cu..63"),Cu_NH3 = starts_with("Cu..65"),Zn = starts_with("Zn"),As_STD = starts_with("As.75"),As_KED = starts_with("As.KED"),As_Oshift = starts_with("As.Oshift"),Se77_Oshift = starts_with("Se..77..Oshift"),Se78_Oshift = starts_with("Se..78..Oshift"),Se80_Oshift = starts_with("Se..80..Oshift"),Se82_Oshift = starts_with("Se..82..Oshift"),Se77_Oshift = starts_with("Se..77..Oshift"),Mo_NH3 = starts_with("Mo..95..NH3"),Mo_Oshift = starts_with("Mo..95..Oshift"),Cd111 = starts_with("Cd..111"),Cd114 = starts_with("Cd..114"),Pb_STD = starts_with("Pb.208"),Pb_NH3 = starts_with("Pb.NH3")))

#rename water3

water3_names <- water3 %>% rename( c(Ca = starts_with("Ca"),Mg = starts_with("Mg"),Al = starts_with("Al"),S = starts_with("S."),P = starts_with("P."),K = starts_with("K"),Na = starts_with("Na"), V = starts_with("V"),Cr_KED =starts_with("Cr..52"),Cr_NH3 = starts_with("Cr..53"),Fe_NH3 = starts_with("Fe..54"),Mn = starts_with("Mn"),Fe_KED = starts_with("Fe..56"),Co = starts_with("Co"),Ni = starts_with("Ni"),Cu_KED = starts_with("Cu..63"),Cu_NH3 = starts_with("Cu..65"),Zn = starts_with("Zn"),As_STD = starts_with("As.75"),As_KED = starts_with("As.KED"),As_Oshift = starts_with("As.Oshift"),Se77_Oshift = starts_with("Se..77..Oshift"),Se78_Oshift = starts_with("Se..78..Oshift"),Se80_Oshift = starts_with("Se..80..Oshift"),Se82_Oshift = starts_with("Se..82..Oshift"),Se77_Oshift = starts_with("Se..77..Oshift"),Mo_NH3 = starts_with("Mo..95..NH3"),Mo_Oshift = starts_with("Mo..95..Oshift"),Cd111 = starts_with("Cd..111"),Cd114 = starts_with("Cd..114"),Pb_STD = starts_with("Pb.208"),Pb_NH3 = starts_with("Pb.NH3")))

#merge water 1 and 3

water_all <- rbind(water_names,water3_names)

# seperate samples from QC and samples into treatments for water 1 ----
dblanks <- water_all %>% filter(grepl("DBLANK", Sample.Id) | grepl("DB$", Sample.Id))

USGS_T239 <- water_all %>% filter(grepl("T-239", Sample.Id) | grepl("239$", Sample.Id))

USGS_M224 <- water_all %>% filter(grepl("M-224", Sample.Id) | grepl("224$", Sample.Id))

bottleblank <- water_all %>% filter( grepl("BB$", Sample.Id))

filterblank <- water_all %>% filter( grepl("FB$", Sample.Id))

labblank <- water_all %>% filter( grepl("LAB", Sample.Id))

MFB <- water_all %>% filter(grepl("MFB$",Sample.Id) | grepl("MFB",Sample.Id))

#merge QC data

water_QAQC <- rbind(dblanks, USGS_T239, USGS_M224,bottleblank,filterblank,MFB, labblank)


# filter out ICPMS QAQC

water_data <- water_all %>% filter(!grepl("^CCV", Sample.Id)) %>%
  filter(!grepl("CCV$", Sample.Id)) %>%
  filter(!grepl("^STD", Sample.Id)) %>%
  filter(!grepl("^std", Sample.Id)) %>%
  filter(!grepl("SP$", Sample.Id)) %>%
  filter(!grepl("SPIKE$", Sample.Id)) %>%
  filter(!grepl("LD$", Sample.Id)) %>%
  filter(!grepl("LSP$", Sample.Id)) %>%
  filter(!grepl("CAL BLANK", Sample.Id)) %>%
  filter(!grepl("^Cal", Sample.Id)) %>%
  filter(!grepl("cal", Sample.Id)) %>%
  filter(!grepl("STD$", Sample.Id)) %>%
  filter(!grepl("^ICV", Sample.Id)) %>%
  filter(!grepl("ICV$", Sample.Id)) %>%
  filter( !grepl("LFB", Sample.Id)) %>%
  filter(!grepl("LBLANK", Sample.Id)) %>%
  filter(!grepl("^LOW", Sample.Id)) %>%
  filter(!grepl("^sac", Sample.Id)) %>%
  filter(!grepl("^LLOQ", Sample.Id))

#filter out all QAQC
qaqcnames<- unique(water_QAQC$Sample.Id)

data <-  water_data %>% filter(!grepl("DBLANK", Sample.Id)) %>%
  filter(!grepl("DB$", Sample.Id))%>%
  filter(!grepl("T-239", Sample.Id)) %>%
  filter(!grepl("239$", Sample.Id)) %>% 
  filter(!grepl("M-224", Sample.Id)) %>%
  filter(!grepl("224$", Sample.Id)) %>% 
  filter( !grepl("BB$", Sample.Id)) %>% 
  filter( !grepl("FB$", Sample.Id)) %>% 
  filter( !grepl("LAB", Sample.Id)) %>% 
  filter(!grepl("MFB$",Sample.Id)) %>%
  filter(!grepl("MFB",Sample.Id))

# some data was duplicated from int. std. reruns so remove duplicates

water_data <- data %>%
  distinct(Sample.Id, .keep_all = TRUE)


###calculate MDLs using digestion blanks ----

summary(dblanks)


# Function to replace outliers with NAs based on z-score

replace_outliers_with_na <- function(x, threshold = 2) {
  z_scores <- abs(scale(x))
  x[z_scores > threshold] <- NA
  return(x)
}

# Apply the function to each column in the data frame

for (col in colnames(dblanks[,c(2:ncol(dblanks))])) {
  dblanks[[col]] <- replace_outliers_with_na(dblanks[[col]])
}



#DF is number of rows -1

#automate with R's quantile function and the number or rows in dataframe

tvalue <- qt(.99,df=(colSums(!is.na(dblanks[,2:34])))-1)

#standard deviation blanks across of all elements 

blanks.stdev <- dblanks[,2:34] %>%
  summarise_all(sd, na.rm = TRUE)

#mean of all elements 

blanks.mean <- dblanks[,2:34] %>%
  summarise_all(sd, na.rm = TRUE)

#rotate data frame to calc mdl

long.stdev <- pivot_longer(blanks.stdev,cols = everything(), names_to = "element", values_to = "stdev")
long.mean <- pivot_longer(blanks.mean,cols = everything(), names_to = "element", values_to = "mean")

#join

long.mdlcalc <- left_join(long.stdev,long.mean, by="element")

#calc mdl

calc.mdl <-long.mdlcalc %>%
  mutate(mdl = stdev*tvalue)

#remove other mean and stdev

mdl <- calc.mdl[,c(1,4)]

#rotate dataframe back 

mdl <- mdl %>%
  pivot_wider(names_from = "element", values_from = "mdl")

### Censor data based on calculated MDL values ----

#chat gpt

process_samples <- function(SAMPLES, mdl) {
  mutate_columns <- colnames(mdl)
  
  for (col in mutate_columns) {
    SAMPLES <- SAMPLES %>%
      mutate(!!col := ifelse(.data[[col]] < mdl[[col]], 0.5 * mdl[[col]], .data[[col]]))
  }
  
  return(SAMPLES)
}

# Assuming you have SAMPLES and mdl data frames already defined
SAMPLES <- water_data
water.mdl <- process_samples(SAMPLES, mdl)



#### calculate true concentrations using sample mass and dilution volume----

#make data long form
water_long <- water.mdl %>% pivot_longer(cols= 2:34, names_to = "Element", values_to = "conc")

# 4.95 mL into 10mL, W samples diluted an extra 5x

water_correct <- water_long %>%
  mutate(Concentration = ifelse(str_detect(Sample.Id, "W"), conc * ((10/4.95)*5), conc * (10/4.95)))


#filter out field dups
water_correct <- water_correct %>%
  filter(!grepl("DUP",Sample.Id))

summary(water_correct)

## Save cleaned data---- 


write.csv(water_correct, file.path("./2_incremental","230731_Water_Metals_long_diff.csv"))

write.csv(water_correct, file.path("./2_incremental","230731_Water_Metals_long.csv"))

## Extract date, size, and site from the Sample.Id----
water_seperate <- water_correct %>%
  # Separate the date, size, and site using str_extract 
  mutate(
    date = str_extract(Sample.Id, "\\d+_\\d+_\\d+"),             # Extract date in the format: 6_21_21
    filtersize = str_extract(Sample.Id, "(W|UF|F)"),                   # Extract size (W, UF, or F)
    site = str_extract(Sample.Id, "(?<=_(W|UF|F)_)(\\d+)")       # Extract the site using lookbehind regex
  ) %>%
  # Convert date to the actual date type using lubridate package
  mutate(date = mdy(date))  # Assuming the format is month-day-year (e.g., 6_21_21)


## calculate size fractions----

water_calc <- water_seperate %>%
  pivot_wider(id_cols = c(site, date,Element),values_from = Concentration, names_from = filtersize)


water_size <- water_calc %>%
  mutate(spm = W - F) %>%
  mutate(colloidal = F - UF) %>%
  mutate(trulydissolved = UF)


##save finished data ----


write.csv(water_size, file.path("./2_incremental","230801_Water_Metals_Size_diff.csv"))


