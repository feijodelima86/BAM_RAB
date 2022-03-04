
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



run1 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","210702_LTREBBiomass_update.csv"), header=T, na.strings = c("")))


run2 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","210714_LTREBBiomass2_update.csv"), header=T, na.strings = c("")))
run3 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","210715_LTREBBiomass3_update.csv"), header=T, na.strings = c("")))
run4 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","210719_LTREBBiomass4_update.csv"), header=T, na.strings = c("")))

run5 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass1_update.csv"), header=T, na.strings = c("")))
run6 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass2_update.csv"), header=T, na.strings = c("")))
run7 <- data.frame(read.csv(file.path("./0_data/parsed_icp_data","BAAMBiomass3_update.csv"), header=T, na.strings = c("")))


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

ICP_MASTER <- full_join(run1,run2)
ICP_MASTER <- full_join(ICP_MASTER,run3)
ICP_MASTER <- full_join(ICP_MASTER,run4)
ICP_MASTER <- full_join(ICP_MASTER,run5)
ICP_MASTER <- full_join(ICP_MASTER,run6)
ICP_MASTER <- full_join(ICP_MASTER,run7)

#force it to be a data frame not a list
ICP_MASTER = data.frame(ICP_MASTER)

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

ICP_full_samples <- ICP_full_samples %>% 
  as.numeric("Al")) %>%
  as.numeric("Al_r"))%>%
  as.numeric("As")%>%
  as.numeric("B")%>%
  as.numeric("B_r")%>%
  as.numeric("Ba")%>%
  as.numeric("Be")%>%
  as.numeric("Ca_r")%>%
  as.numeric("Cd")%>%
  as.numeric("Co")%>%
  as.numeric("Cr")%>%
  as.numeric("Cu")%>%
  as.numeric("Fe")%>%
  as.numeric("Fe_r238")%>%
  as.numeric("Fe_r234")%>%
  as.numeric("K_r")%>%
  as.numeric("Li_r")%>%
  as.numeric("Mg_r")%>%
  as.numeric("Mn")%>%
  as.numeric("Mn_r")%>%
  as.numeric("Mo")%>%
  as.numeric("Na_r")%>%
  as.numeric("Ni")%>%
  as.numeric("P")%>%
  as.numeric("Pb")%>%
  as.numeric("S180")%>%
  as.numeric("S181")%>%
  as.numeric("S_r")%>%
  as.numeric("Sb")%>%
  as.numeric("Se")%>%
  as.numeric("Si")%>%
  as.numeric("Si_r251")%>%
  as.numeric("Si_r212")%>%
  as.numeric("Sn")%>%
  as.numeric("Sr_r")%>%
  as.numeric("Ti")%>%
  as.numeric("Tl")%>%
  as.numeric("Zn")%>%
  as.numeric("Zn_r")
  
  
  
#check transformation
summary(ICP_full_samples)


#seperate unknown samples from calibration and qaqc samples
# last two unknown samples from each run (execept the first run) are unlabeled qaqc samples. however, keep with samples until sample info is brought in

#unknonwn samples + dups and sample spikes (last 2 for every selection)
sample1 <- run1[14:25, ]
sample2 <- run2[c(24:35,39:50,54:65,69:80), ]
sample3 <- run3[c(24:35,39:50,54:65,69:71,73:81,85:96), ]
sample4 <- run4[c(14:25,29:40,44:49), ]

#rename run 1 samples so they are the same format as newer samples
#rownames(sample1) <- c("26","27","28","29","30","50", "49", "54", "57","60","26spike","26dup")




#qaqc samples CV,CCV, analysis spikes
qaqc1 <- run1[c(11:13,26:28), ]
qaqc2 <- run2[c(19:23,36:38,51:53,66:68,81:83), ]
qaqc3 <- run3[c(19:23,24:28,51:53,66:68,82:84,97:99), ]
qaqc4 <- run4[c(11:13,26:28,41:43,50:52), ]


#IDL Blanks
blanks1 <- run1[c(12,28),]
blanks2 <- run2[c(22,38,53,68,83),]
blanks3 <- run3[c(22,38,53,68,84,99),]
blanks4 <- run4[c(12,28,43,52),]

#MDL Blanks (digestion blanks) are not run everytime, search for them and add them to master MDL
MDL.run1 <- run1[23, ]
MDL.run3 <- run3[c(24,25,40:42), ]


### read in sample information file
sample.meta <- read.csv(file.path("./0_data/sample_info","080421_SAMPLE_LOG.csv"),header=T, na.strings = c(""))

#join sample metadata and ICP data together
summary(sample.meta)
summary(samples.full)
sample.meta$ICP_NUMBER <-as.character(sample.meta$ICP_NUMBER)

raw_data <- inner_join(sample.meta,samples.full, by = c("ICP_NUMBER"="Sample.ID","ICP_DATE"="Date"))

#check to make sure join was successful


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

#combine all MDLs from different runs

mdl.master <- rbind(MDL.run1,MDL.run3)

#get rid of labeling columns leaving only element columns

mdl.master <- mdl.master[, c(6:45)]

#calculate MDL off of digestion blanks MDL= students t-value (99% CI for n-1 DF)*standard deviation. t value is 3.14 for 7 replicates

#students t-test value for number of blanks used 
#DF is number of rows -1
tvalue <- qt(.99,df=(length(mdl.master[ ,1])-1))

#standard deviation blanks across of all elements and beam lines
mdl.master.calc<-mdl.master %>%
  summarise_all(sd)

#rotate data frame to calc mdl

long.mdl <- pivot_longer(mdl.master.calc,cols = everything(), names_to = "element", values_to = "stdev")

#calc mdl
long.mdl.calc <-long.mdl %>%
  mutate(mdl = stdev*tvalue)

#rotate dataframe back 

mdl.full <- long.mdl.calc %>%
  pivot_wider(id_cols=c(1,3),names_from = "element", values_from = "mdl")

#select beamlines based on lowest




### Censor data based on calculated MDL values ----

#combine sample data
samples.full <- rbind(sample1,sample2,sample3,sample4)

# censor data based off calculated MDL 

samples.cen <- samples.full %>%
  mutate(as_cen = ifelse(as_raw < mdl$as,
                         0.5*mdl$as,
                         as_raw)) %>%
  mutate(cd_cen = ifelse(cd_raw < mdl$cd,
                         0.5*mdl$cd,
                         cd_raw)) %>%
  mutate(cu_cen = ifelse(cu_raw < mdl_$cu,
                         0.5*mdl$cu,
                         cu_raw)) %>%
  mutate(zn_cen = ifelse(zn_raw < mdl$zn,
                         0.5*mdl$zn,
                         zn_raw)) %>%
  mutate(pb_cen = ifelse(pb_raw < mdl$pb,
                         0.5*mdl$pb,
                         pb_raw))



