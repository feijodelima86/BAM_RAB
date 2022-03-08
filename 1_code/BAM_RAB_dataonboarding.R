
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
options(scipen=100)

### I. Read in ICP files----

run1 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","210702_LTREBBiomass_update.csv"), header=T, na.strings = c("")))

run2 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","210714_LTREBBiomass2_update.csv"), header=T, na.strings = c("")))
run3 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","210715_LTREBBiomass3_update.csv"), header=T, na.strings = c("")))
run4 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","210719_LTREBBiomass4_update.csv"), header=T, na.strings = c("")))

run5 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass1_update.csv"), header=T, na.strings = c("")))
run6 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass2_update.csv"), header=T, na.strings = c("")))
run7 <- as.data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass3_update.csv"), header=T, na.strings = c("")))

### Change data types, names, and remove calibration standards ----

#rename columns and elements 
colnames(run1) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER", "Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,	"S181" ,	"S_r",	"S180",	"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
        
colnames(run2) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER","Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,	"S181" ,	"S_r",	"S180",	"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
colnames(run3) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER", "Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,	"S181" ,	"S_r",	"S180",	"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
colnames(run4) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER", "Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,	"S181" ,	"S_r",	"S180",	"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
colnames(run5) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER","Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,"S180",	"S181" ,	"S_r",		"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
colnames(run6) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks" ,	"RUN_NUMBER", "Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,"S180",	"S181" ,	"S_r",		"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)
colnames(run7) <- c("SAMPLE_ID",	"ICP_DATE",	"TIME",	"Remarks",	"RUN_NUMBER","Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,"S180",	"S181" ,	"S_r",		"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" 
)

#now that column names are the same join dataframes

ICP_MASTER <- rbind.data.frame(run1,run2)
ICP_MASTER <- rbind.data.frame(ICP_MASTER,run3)
ICP_MASTER <- rbind.data.frame(ICP_MASTER,run4)
ICP_MASTER <- rbind.data.frame(ICP_MASTER,run5)
ICP_MASTER <- rbind.data.frame(ICP_MASTER,run6)
ICP_MASTER <- rbind.data.frame(ICP_MASTER,run7)



#remove row denoting units as it forces character
ICP_MASTER <- filter(ICP_MASTER,SAMPLE_ID !="NA")

#filter out calibration standards for same reason
ICP_full_samples <- ICP_MASTER %>% 
  filter(SAMPLE_ID != "Calib Blank 1")%>%
  filter(SAMPLE_ID != "SedSTD_1")%>%
  filter(SAMPLE_ID != "SedSTD_2")%>%
  filter(SAMPLE_ID != "SedSTD_3")%>%
  filter(SAMPLE_ID != "SedSTD_4")%>%
  filter(SAMPLE_ID != "SedSTD_5")%>%
  filter(SAMPLE_ID != "SedSTD_S")%>%
  filter(SAMPLE_ID != "SED Std 1")%>%
  filter(SAMPLE_ID != "SED Std 2")%>%
  filter(SAMPLE_ID != "SED Std 3")%>%
  filter(SAMPLE_ID != "SED Std 4")%>%
  filter(SAMPLE_ID != "SED Std 5")%>%
  filter(SAMPLE_ID != "SED Std S")


#change concentrations to numeric

ICP_samples <- as.data.frame(ICP_full_samples,stringsAsFactors = FALSE)

class(ICP_samples)

# is data frame but as.numeric doesn't like me and thinks its a list when I pipeline ----

#ICP_samples <- ICP_samples %>% 
 # as.numeric("Al") %>%
  #as.numeric("Al_r")%>%
  #as.numeric("As")%>%
  #as.numeric("B")%>%
  #as.numeric("B_r")%>%
  #as.numeric("Ba")%>%
  #as.numeric("Be")%>%
  #as.numeric("Ca_r")%>%
  #as.numeric("Cd")%>%
  #as.numeric("Co")%>%
  #as.numeric("Cr")%>%
  #as.numeric("Cu")%>%
  #as.numeric("Fe")%>%
  #as.numeric("Fe_r238")%>%
  #as.numeric("Fe_r234")%>%
  #as.numeric("K_r")%>%
  #as.numeric("Li_r")%>%
  #as.numeric("Mg_r")%>%
  #as.numeric("Mn")%>%
  #as.numeric("Mn_r")%>%
  #as.numeric("Mo")%>%
  #as.numeric("Na_r")%>%
  #as.numeric("Ni")%>%
  #as.numeric("P")%>%
  #as.numeric("Pb")%>%
  #as.numeric("S180")%>%
  #as.numeric("S181")%>%
  #as.numeric("S_r")%>%
  #as.numeric("Sb")%>%
  #as.numeric("Se")%>%
  #as.numeric("Si")%>%
  #as.numeric("Si_r251")%>%
  #as.numeric("Si_r212")%>%
  #as.numeric("Sn")%>%
  #as.numeric("Sr_r")%>%
  #as.numeric("Ti")%>%
  #as.numeric("Tl")%>%
  #as.numeric("Zn")%>%
  #as.numeric("Zn_r")


#instead select columns----
ICP_samples[,6:44] <- lapply(ICP_samples[,6:44], as.numeric) #these are the columns with concentrations 
  

#check transformation
summary(ICP_samples)

### Merge concentration data with sample information and masses----

#load in sample information, masses are already attached
info_mass <- as.data.frame(read.csv(file.path("./0_data/sample_mass","20220307_SAMPLE_LOG_MASS.csv"), header=T, na.strings = c("")))

#column types need to be identical to join
summary(info_mass)

summary(ICP_samples)

#change Run number to character as its descriptive
ICP_samples$RUN_NUMBER <- as.character(ICP_samples$RUN_NUMBER) 

#change ICP_Number to character as itsdescriptive
info_mass$ICP_NUMBER <- as.character(info_mass$ICP_NUMBER)


#Check transformation
summary(ICP_samples)
summary(info_mass)

#join sample information to ICP concentrations
SAMPLE_MASTER <- left_join(ICP_samples, info_mass, by = c( "ICP_DATE" = "OES_DATE", "SAMPLE_ID" = "ICP_NUMBER"), keep = TRUE )

#Inspect new file to make sure there are no columns with .y and .x which indicates incomplete merge
#is oldest ICP date corresponing to a sampling date of 6/22


#identify column numbers
summary(SAMPLE_MASTER)
#Remove unneeded columns
prunedSAMPLE_MASTER <- select(SAMPLE_MASTER,c(SAMPLE_ID, ICP_DATE, Al:Zn_r, SAMPLING_DATE:OES_DATE, ICP_NUMBER, LAB.DUPLICATE:ICP.SPIKE, DRY_MASS, NOTES))


### filter QAQC samples----

QAQC_samples <- prunedSAMPLE_MASTER %>%
  filter(SAMPLE_ID == "CCV" | SAMPLE_ID == "CV" | SAMPLE_ID == "S IPC" | SAMPLE_ID == "Lblank" | SAMPLE_ID == "LFB" | SAMPLE_ID == "cal std 1" | SAMPLE_ID == "cal std 5" | SAMPLE_ID == "cal blank" | SAMPLE_DESCRIPTOR == "STSD2" | SAMPLE_DESCRIPTOR == "METHOD FORTIFIED BLANK" | SAMPLE_DESCRIPTOR == "LABORATORY FORTIFIED SAMPLE" | SAMPLE_DESCRIPTOR == "DIGEST BLANK" | SAMPLE_DESCRIPTOR == "FILTER BLANK" | SAMPLE_ID == "MFB 1.1" | SAMPLE_ID == "MFB 1.8" | SAMPLE_ID == "MFB 2.1" | SAMPLE_ID == "MFB 2.8")

#save off QAQC file to incremental
write.csv(QAQC_samples,"2_incremental\\20220307_QAQC_SAMPLES.csv", row.names = FALSE)


duplicates <- prunedSAMPLE_MASTER %>% 
  filter(grepl("dup$",SAMPLE_ID )) #grepl searches for string matches in the vector $ goes after for ends with or before for begins with. This can also be accomplished using strdetect or strsub

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




### Inital QAQC tests ----


#using the original imported files compare analysis stage duplicates.
run2perdiff33 <- run2[33,6:45]
run2perdiff34 <- run2[34,6:45]

#
perdiff2.1 <- run2 %>%
  mutate(pct_diff = (Profit/lead(Profit) - 1) * 100)
    


#SRM
#MFB
#vaiance across digests

#how
### Calculate detection limits and select beam lines to use ----

# beamlines are selected based on accuracy to CCV and how low their detection limit is 

#elements that have multiple beamlines Al(2), B(2), Fe(3), Mn(2), S(3), Si(3), Tl(2), Zn(2)


###calculate MDLs

# load in QAQC data
qaqc.samples <- as.data.frame(read.csv(file.path("./2_incremental","20220307_QAQC_SAMPLES.csv"), header=T, na.strings = c("")))


#select Digestion blank
blanks <- qaqc.samples %>% filter(SAMPLE_DESCRIPTOR == "DIGEST BLANK")

#17 blanks

#calculate MDL off of digestion blanks MDL= students t-value (99% CI for n-1 DF)*standard deviation. t value is 3.14 for 7 replicates

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
long.mdl <-long.stdev %>%
  mutate(mdl = stdev*tvalue)

#rotate dataframe back 

mdl <- long.mdl %>%
  pivot_wider(id_cols=c(1,3),names_from = "element", values_from = "mdl")

#select beamlines based on lowest




### Censor data based on calculated MDL values ----

#make a list of all elements
#a <- c(Al,	Al_r,	As,	B,	B_r, 	Ba,	Be,	Ca_r,	Cd,	Co,	Cr, Cu, Fe,	Fe_r238 ,	Fe_r234 ,	K_r,	Li_r,	Mg_r,	Mn ,	Mn_r ,	Mo,	Na_r ,	Ni ,	P ,	Pb ,	S181 ,	S_r,	S180,	Sb ,	Se ,	Si, Si_r251, Si_r212, Sn, Sr_r, Ti, Tl, Zn, Zn_r )

a <-c("Al",	"Al_r",	"As",	"B",	"B_r", 	"Ba",	"Be",	"Ca_r",	"Cd",	"Co",	"Cr", "Cu", "Fe",	"Fe_r238" ,	"Fe_r234" ,	"K_r",	"Li_r",	"Mg_r",	"Mn" ,	"Mn_r" ,	"Mo",	"Na_r" ,	"Ni" ,	"P" ,	"Pb" ,	"S181" ,	"S_r",	"S180",	"Sb" ,	"Se" ,	"Si","Si_r251","Si_r212","Sn","Sr_r","Ti","Tl","Zn","Zn_r" )


for (i in a) {
  mdl.censored.samples <- SAMPLES[,3:41] %>%
    mutate_("i_cen = ifelse(i < mdl$i,
                           0.5*mdl$i,
                           i)") 
  
}




# Create example data frame
df <- data.frame(c(1, 2, 3, 4, 5, 6, 7))
colnames(df) <- c("var")

# Your sequence of exponents
alpha <- seq(0.85, 0.95, by= .01)

# Call to the function, assign the return value to df2
df2 <- dblExponential(alpha,df$var)

# print df2
df2

SAMPLES_CENSORED <- ()


samples.limit <- SAMPLES %>%
  mutate(As_cen = ifelse(as_raw < mdl$As,
                         0.5*mdl$as,
                         As)) %>%
  mutate(cd_cen = ifelse(cd_raw < mdl$cd,
                         0.5*mdl$Cd,
                         Cd)) %>%
  mutate(cu_cen = ifelse(cu_raw < mdl_$cu,
                         0.5*mdl$Cu,
                         Cu_raw)) %>%
  mutate(zn_cen = ifelse(zn_raw < mdl$zn,
                         0.5*mdl$zn,
                         zn_raw)) %>%
  mutate(pb_cen = ifelse(pb_raw < mdl$pb,
                         0.5*mdl$pb,
                         pb_raw))






# Messing around with your code. MWAHAHAHAHAHAHA!
go plot yourself!!!!