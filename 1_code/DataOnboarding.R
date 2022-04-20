### This is the script for post processing ICP data for the BAM-RAB LTREB Biomass project

# Data consists of minimally processed ICP-OES data and excell files containing sample mass
# There is an ICP data file with unique detection limits for each run.
# There are multiple runs for each collection date

# this metal concentration data will be combined with biomass and chrolorphyll data to relate metal concentrations to primary productivity

### I. read in ICP files
### II. bring in sample information file
### III. QAQC tests
### IV. select beamlines and calculate detection limits (MDL)
### V. censor data bellow MDL values 
### VI. Bring in sample mass data and calculate per gram concentrations


library(tidyverse)
library(naniar)
library(ggplot2)
library(openxlsx)
library(readxl)
library(lubridate) 
options(scipen=100)

### Main Functions ----


### Currently not needed
### Function to match data frame classes

matchColClasses <- function(df1, df2) {
  
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  
  return(df2)
}



### Read in ICP files as a list and merge runs together ----




### Set a working file path for OES spreadsheets. 

path1<-"0_Data/updated_parsed_icp_data"

### This command finds all spreadsheets with Biomass in their name within the folder specified above.

ICP.LIST<-c(list.files(path= path1, full.names=TRUE,pattern=c("Biomass")))

### This command creates a list of objects with all those files.

ICP<-lapply(ICP.LIST,function(filename){
  print(paste("Merging",filename,sep.=""))
  data.frame(read.csv(filename, fileEncoding="UTF-8-BOM"))
})

### This removes the first line that contains only "mg/l" from all data frames within the list.

ICP <- lapply(ICP,function(x){
  x[-c(1,0),]
})

### This is a trick to remove numbers from all column headings. The make.unique function then adds unique identifiers to elements that are repeated.
### E.G. Iron has 3 types of detection (1 axial 2 radial), Those become Fe, Fer and Fer1. 

ICP <- lapply(ICP, function(x) 
{names(x) <- make.unique(gsub("[^[:alpha:]]", "", names(x)));x})

### This is where all spreadsheets are merged into one. 

ICP.ALL <- ICP %>% reduce(full_join)

#replace blanks with NAs 
ind <- which(ICP.ALL == "", arr.ind = TRUE)
ICP.ALL[ind] <- NA

#change concentration columns to numeric ----

summary(ICP.ALL) #currently character becuase of calibrations

#as.numeric doesn't like pipe, so use lapply
ICP.ALL[,6:44] <- lapply(ICP.ALL[,6:44], as.numeric) #these are the columns with concentrations 
#nans introduced by coersion from the calibration standards

#check transformation
summary(ICP.ALL)


### Remove calibration standards ----

#remove row denoting units and calibration standards 
ICP_samples <- ICP.ALL %>% 
  filter(SampleID !="NA")%>%
  filter(SampleID != "Calib Blank 1")%>%
  filter(SampleID != "SedSTD_1")%>%
  filter(SampleID != "SedSTD_2")%>%
  filter(SampleID != "SedSTD_3")%>%
  filter(SampleID != "SedSTD_4")%>%
  filter(SampleID != "SedSTD_5")%>%
  filter(SampleID != "SedSTD_S")%>%
  filter(SampleID != "SED Std 1")%>%
  filter(SampleID != "SED Std 2")%>%
  filter(SampleID != "SED Std 3")%>%
  filter(SampleID != "SED Std 4")%>%
  filter(SampleID != "SED Std 5")%>%
  filter(SampleID != "SED Std S")


### Merge concentration data with sample information and masses----

#load in sample information, masses are already attached
info_mass <- as.data.frame(read.csv(file.path("./0_data/sample_mass","20220307_SAMPLE_LOG_MASS.csv"), header=T, na.strings = c("")))

### Variable Classes (column types) need to be identical to join
#what are the variable classes
summary(info_mass)
summary(ICP_samples)

#currently all the same so proceed


#join sample information to ICP concentrations
SAMPLE_MASTER <- left_join(ICP_samples, info_mass, by = c( "Date" = "OES_DATE", "SampleID" = "ICP_NUMBER"), keep = TRUE )

#Inspect new file to make sure there are no columns with .y and .x which indicates incomplete merge
#is oldest ICP date corresponing to a sampling date of 6/22


#identify column numbers
summary(SAMPLE_MASTER)
#Remove unneeded columns 
# RFL: Final version of mass file had column "NOTES" removed and was making line 174 fail. Removed it.

prunedSAMPLE_MASTER <- select(SAMPLE_MASTER,c(SampleID, Date, Al:Znr, SAMPLING_DATE:OES_DATE, ICP_NUMBER, LAB.DUPLICATE:ICP.SPIKE, DRY_MASS))

### filter QAQC samples----

QAQC_samples <- prunedSAMPLE_MASTER %>%
  filter(SampleID == "CCV" | SampleID == "CV" | SampleID == "S IPC" | SampleID == "Lblank" | SampleID == "LFB" | SampleID == "cal std 1" | SampleID == "cal std 5" | SampleID == "cal blank" | SAMPLE_DESCRIPTOR == "STSD2" | SAMPLE_DESCRIPTOR == "METHOD FORTIFIED BLANK" | SAMPLE_DESCRIPTOR == "LABORATORY FORTIFIED SAMPLE" | SAMPLE_DESCRIPTOR == "DIGEST BLANK" | SAMPLE_DESCRIPTOR == "FILTER BLANK" | SampleID == "MFB 1.1" | SampleID == "MFB 1.8" | SampleID == "MFB 2.1" | SampleID == "MFB 2.8")

#save off QAQC file to incremental
write.csv(QAQC_samples,"2_incremental\\20220307_QAQC_SAMPLES.csv", row.names = FALSE)



###sample file

SAMPLES <- prunedSAMPLE_MASTER %>%
  filter(SAMPLE_DESCRIPTOR != "STSD2" ) %>%
  filter(DRY_MASS != "NA" ) %>%
  filter(SAMPLE_DESCRIPTOR != "STSD2" ) %>%
  filter(SAMPLE_DESCRIPTOR != "METHOD FORTIFIED BLANK" ) %>%
  filter(SAMPLE_DESCRIPTOR != "LABORATORY FORTIFIED SAMPLE" ) %>%
  filter(SAMPLE_DESCRIPTOR != "DIGEST BLANK" ) %>%
  filter(SAMPLE_DESCRIPTOR != "FILTER BLANK") 


#save off sample file

write.csv(SAMPLES,"2_incremental\\20220307_FIELD_SAMPLES.csv", row.names = FALSE)



### Calculate MDL values and apply them to censor concentrations bellow the detection limit ----


###calculate MDLs

# load in QAQC data
qaqc.samples <- as.data.frame(read.csv(file.path("./2_incremental","20220307_QAQC_SAMPLES.csv"), header=T, na.strings = c("")))

summary(qaqc.samples)

#select Digestion blank
blanks <- qaqc.samples %>% filter(SAMPLE_DESCRIPTOR == "DIGEST BLANK")

#17 blanks

#calculate MDL off of digestion blanks MDL = students t-value (99% CI for n-1 DF)*standard deviation. t value is 3.14 for 7 replicates

#DF is number of rows -1

#automate with R's quantile function and the number or rows in dataframe
tvalue <- qt(.99,df=(length(blanks[ ,1])-1))
#hard calculation t value is 2.583 for n=17 (df=16) double check automate worked


#standard deviation blanks across of all elements and beam lines
blanks.stdev <- blanks[,3:41] %>%
  summarise_all(sd)

#rotate data frame to calc mdl
long.stdev <- pivot_longer(blanks.stdev,cols = everything(), names_to = "element", values_to = "stdev")

#calc mdl
calc.mdl <-long.stdev %>%
  mutate(mdl = stdev*tvalue)

#rotate dataframe back 

mdl <- calc.mdl %>%
  pivot_wider(id_cols=c(1,3),names_from = "element", values_from = "mdl")

#rotate long again

long.mdl <- mdl %>% pivot_longer(cols = everything(), names_to = "element", values_to = "MDL")

#select beamlines based on lowest



### Censor data based on calculated MDL values ----


mdl.samples <- SAMPLES %>%
  mutate(As = ifelse(As < mdl$As,
                     0.5*mdl$As,
                     As)) %>%
  mutate(Cd = ifelse(Cd < mdl$Cd,
                     0.5*mdl$Cd,
                     Cd)) %>%
  mutate(Cu = ifelse(Cu < mdl$Cu,
                     0.5*mdl$Cu,
                     Cu)) %>%
  mutate(Zn = ifelse(Zn < mdl$Zn,
                     0.5*mdl$Zn,
                     Zn)) %>%
  mutate(Pb = ifelse(Pb < mdl$Pb,
                     0.5*mdl$Pb,
                     Pb)) %>%
  mutate(Al = ifelse(Al < mdl$Al,
                     0.5*mdl$Al,
                     Al)) %>%
  mutate(Alr = ifelse(Alr < mdl$Alr,
                      0.5*mdl$Alr,
                      Alr)) %>%
  mutate(B = ifelse(B < mdl$B,
                    0.5*mdl$B,
                    B)) %>%
  mutate(Br = ifelse(Br < mdl$Br,
                     0.5*mdl$Br,
                     Br)) %>%
  mutate(Ba = ifelse(Ba < mdl$Ba,
                     0.5*mdl$Ba,
                     Ba)) %>%
  mutate(Be = ifelse(Be < mdl$Be,
                     0.5*mdl$Be,
                     Be)) %>%
  mutate(Car = ifelse(Car < mdl$Car,
                      0.5*mdl$Car,
                      Car)) %>%
  mutate(Co = ifelse(Co < mdl$Co,
                     0.5*mdl$Co,
                     Co)) %>%
  mutate(Cr = ifelse(Cr < mdl$Cr,
                     0.5*mdl$Cr,
                     Cr)) %>%
  mutate(Fe = ifelse(Fe < mdl$Fe,
                     0.5*mdl$Fe,
                     Fe)) %>%
  mutate(Fer = ifelse(Fer < mdl$Fer,
                      0.5*mdl$Fer,
                      Fer)) %>%
  mutate(Fer.1 = ifelse(Fer.1 < mdl$Fer.1,
                        0.5*mdl$Fer.1,
                        Fer.1)) %>%
  mutate(Kr = ifelse(Kr < mdl$Kr,
                     0.5*mdl$Kr,
                     Kr)) %>%
  mutate(Lir = ifelse(Lir < mdl$Lir,
                      0.5*mdl$Lir,
                      Lir)) %>%
  mutate(Mgr = ifelse(Mgr < mdl$Mgr,
                      0.5*mdl$Mgr,
                      Mgr)) %>%
  mutate(Mn = ifelse(Mn < mdl$Mn,
                     0.5*mdl$Mn,
                     Mn)) %>%
  mutate(Mnr = ifelse(Mnr < mdl$Mnr,
                      0.5*mdl$Mnr,
                      Mnr)) %>%
  mutate(Mo = ifelse(Mo < mdl$Mo,
                     0.5*mdl$Mo,
                     Mo)) %>%
  mutate(Nar = ifelse(Nar < mdl$Nar,
                      0.5*mdl$Nar,
                      Nar)) %>%
  mutate(Ni = ifelse(Ni < mdl$Ni,
                     0.5*mdl$Ni,
                     Ni)) %>%
  mutate(P = ifelse(P < mdl$P,
                    0.5*mdl$P,
                    P)) %>%
  mutate(S = ifelse(S < mdl$S,
                    0.5*mdl$S,
                    S)) %>%
  mutate(S.1 = ifelse(S.1 < mdl$S.1,
                      0.5*mdl$S.1,
                      S.1)) %>%
  mutate(Sr = ifelse(Sr < mdl$Sr,
                     0.5*mdl$Sr,
                     Sr)) %>%
  mutate(Sb = ifelse(Sb < mdl$Sb,
                     0.5*mdl$Sb,
                     Sb)) %>%
  mutate(Se = ifelse(Se < mdl$Se,
                     0.5*mdl$Se,
                     Se)) %>%
  mutate(Si = ifelse(Si < mdl$Si,
                     0.5*mdl$Si,
                     Si)) %>%
  mutate(Sir = ifelse(Sir < mdl$Sir,
                      0.5*mdl$Sir,
                      Sir)) %>%
  mutate(Sir.1 = ifelse(Sir.1 < mdl$Sir.1,
                        0.5*mdl$Sir.1,
                        Sir.1)) %>%
  mutate(Sn = ifelse(Sn < mdl$Sn,
                     0.5*mdl$Sn,
                     Sn)) %>%
  mutate(Srr = ifelse(Srr < mdl$Srr,
                      0.5*mdl$Srr,
                      Srr)) %>%
  mutate(Ti = ifelse(Ti < mdl$Ti,
                     0.5*mdl$Ti,
                     Ti)) %>%
  mutate(Tl = ifelse(Tl < mdl$Tl,
                     0.5*mdl$Tl,
                     Tl)) %>%
  mutate(Znr = ifelse(Znr < mdl$Znr,
                      0.5*mdl$Znr,
                      Znr)) 

summary(mdl.samples)

### Save file of censored mdl samples

write.csv(mdl.samples,"2_incremental\\20220312_MDL_SAMPLES.csv", row.names = FALSE)


#create BDL dataframe to verify code is working and visually inspect detection limit issues


bdl.samples <- SAMPLES %>%
  mutate(As = ifelse(As < mdl$As,
                     "BDL",
                     As)) %>%
  mutate(Cd = ifelse(Cd < mdl$Cd,
                     "BDL",
                     Cd)) %>%
  mutate(Cu = ifelse(Cu < mdl$Cu,
                     "BDL",
                     Cu)) %>%
  mutate(Zn = ifelse(Zn < mdl$Zn,
                     "BDL",
                     Zn)) %>%
  mutate(Pb = ifelse(Pb < mdl$Pb,
                     "BDL",
                     Pb)) %>%
  mutate(Al = ifelse(Al < mdl$Al,
                     "BDL",
                     Al)) %>%
  mutate(Alr = ifelse(Alr < mdl$Alr,
                      "BDL",
                      Alr)) %>%
  mutate(B = ifelse(B < mdl$B,
                    "BDL",
                    B)) %>%
  mutate(Br = ifelse(Br < mdl$Br,
                     "BDL",
                     Br)) %>%
  mutate(Ba = ifelse(Ba < mdl$Ba,
                     "BDL",
                     Ba)) %>%
  mutate(Be = ifelse(Be < mdl$Be,
                     "BDL",
                     Be)) %>%
  mutate(Car = ifelse(Car < mdl$Car,
                      "BDL",
                      Car)) %>%
  mutate(Co = ifelse(Co < mdl$Co,
                     "BDL",
                     Co)) %>%
  mutate(Cr = ifelse(Cr < mdl$Cr,
                     "BDL",
                     Cr)) %>%
  mutate(Fe = ifelse(Fe < mdl$Fe,
                     "BDL",
                     Fe)) %>%
  mutate(Fer = ifelse(Fer < mdl$Fer,
                      "BDL",
                      Fer)) %>%
  mutate(Fer.1 = ifelse(Fer.1 < mdl$Fer.1,
                        "BDL",
                        Fer.1)) %>%
  mutate(Kr = ifelse(Kr < mdl$Kr,
                     "BDL",
                     Kr)) %>%
  mutate(Lir = ifelse(Lir < mdl$Lir,
                      0.5*mdl$Lir,
                      Lir)) %>%
  mutate(Mgr = ifelse(Mgr < mdl$Mgr,
                      0.5*mdl$Mgr,
                      Mgr)) %>%
  mutate(Mn = ifelse(Mn < mdl$Mn,
                     0.5*mdl$Mn,
                     Mn)) %>%
  mutate(Mnr = ifelse(Mnr < mdl$Mnr,
                      0.5*mdl$Mnr,
                      Mnr)) %>%
  mutate(Mo = ifelse(Mo < mdl$Mo,
                     0.5*mdl$Mo,
                     Mo)) %>%
  mutate(Nar = ifelse(Nar < mdl$Nar,
                      0.5*mdl$Nar,
                      Nar)) %>%
  mutate(Ni = ifelse(Ni < mdl$Ni,
                     0.5*mdl$Ni,
                     Ni)) %>%
  mutate(P = ifelse(P < mdl$P,
                    0.5*mdl$P,
                    P)) %>%
  mutate(S = ifelse(S < mdl$S,
                    0.5*mdl$S,
                    S)) %>%
  mutate(S.1 = ifelse(S.1 < mdl$S.1,
                      0.5*mdl$S.1,
                      S.1)) %>%
  mutate(Sr = ifelse(Sr < mdl$Sr,
                     0.5*mdl$Sr,
                     Sr)) %>%
  mutate(Sb = ifelse(Sb < mdl$Sb,
                     0.5*mdl$Sb,
                     Sb)) %>%
  mutate(Se = ifelse(Se < mdl$Se,
                     0.5*mdl$Se,
                     Se)) %>%
  mutate(Si = ifelse(Si < mdl$Si,
                     0.5*mdl$Si,
                     Si)) %>%
  mutate(Sir = ifelse(Sir < mdl$Sir,
                      0.5*mdl$Sir,
                      Sir)) %>%
  mutate(Sir.1 = ifelse(Sir.1 < mdl$Sir.1,
                        0.5*mdl$Sir.1,
                        Sir.1)) %>%
  mutate(Sn = ifelse(Sn < mdl$Sn,
                     0.5*mdl$Sn,
                     Sn)) %>%
  mutate(Srr = ifelse(Srr < mdl$Srr,
                      0.5*mdl$Srr,
                      Srr)) %>%
  mutate(Ti = ifelse(Ti < mdl$Ti,
                     0.5*mdl$Ti,
                     Ti)) %>%
  mutate(Tl = ifelse(Tl < mdl$Tl,
                     0.5*mdl$Tl,
                     Tl)) %>%
  mutate(Znr = ifelse(Znr < mdl$Znr,
                      0.5*mdl$Znr,
                      Znr)) 





### Inital QAQC tests ----

duplicates <- prunedSAMPLE_MASTER %>% 
  filter(grepl("dup$",SampleID ) | grepl("DUP$",SampleID))


#grepl searches for string matches in the vector $ goes after for ends with or before for begins with. This can also be accomplished using strdetect or strsub
#grepl returns locigal vector
#apropos returns character vector

#mutate case_when could work

full.dup <- SAMPLES %>% filter(grepl()
                               
                               #
                               perdiff2.1 <- run2 %>%
                                 mutate(pct_diff = (Profit/lead(Profit) - 1) * 100)
                               
                               
                               #SRM
                               #MFB
                               #vaiance across digests
                               
                               #how
                               
                               
                               
                               
                               
                               ### ask stack exchange about looping
                               x=c(1,2,3,4,5)
                               y=c(3,4,1,7,8)
                               z=c("ab","cd","bc","cd","de")
                               samples <- 
                                 
                                 for (i in seq_along(SAMPLES)){
                                   loop.mdl.samples <- SAMPLES %>%
                                     mutate("MDL.{{i}}" := ifelse(SAMPLES[i] < mdl[[i]], 0.5*mdl[[i]], SAMPLES[i]))}
                               
                               
                               
                               mdl.samples <- SAMPLES[,3:41] %>%
                                 lapply( function(x){
                                   mutate("{x}.mdl" := ifelse(only.samples < long.mdl, 0.5*mdl[[i]], only.samples[i]))long.mdl[as.character(Year)])
                                 })

#samples.mdl <- SAMPLES %>%
# mutate(across(3:41, ifelse(.[[i]] < mdl$[[i]], 0.5*mdl$[[i]], .[[i]]), .names = "MDL_{.col}")) 




# Messing around with your code. MWAHAHAHAHAHAHA!
go plot yourself!!!!
  