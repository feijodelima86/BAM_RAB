# This script if for analyzing QAQC data from ICP instrumnetations


# load packages ----
library(tidyverse)
library(naniar)
library(ggplot2)
library(openxlsx)
library(readxl)
library(lubridate) 
options(scipen=100)

# load in data ----

bam_data <- as.data.frame(read.csv(file.path("./2_incremental","20220307_QAQC_SAMPLES.csv"), header=T, na.strings = c("")))

SRM_values <- as.data.frame(read.csv(file.path("./0_data","STSD2_reported_values.csv"), header=T, na.strings = c("")))


#adjust column types as needed
summary(bam_data)
bam_data <- bam_data %>% mutate(DRY_MASS = as.numeric(DRY_MASS))%>% 
  mutate(Al = as.numeric(Al))
summary(bam_data)

# filter types of QAQC
srms <- bam_data %>% filter(SAMPLE_DESCRIPTOR == "STSD2")


# SRMs ----
#calculate conc per mass of srm
#remove srms with missing dry mass
srms <- srms %>% drop_na(DRY_MASS)

#pivot to calculate concertration per gram
srms_long <- as.data.frame(pivot_longer(srms, c(3:41), names_to = "element", values_to = "concentration"))

srms_long <- srms_long %>% mutate(DRY_CONC = (concentration * 0.05 )/DRY_MASS)
summary(srms_long)
srms_long <- srms_long %>% drop_na(DRY_CONC)

# adjust srm values to match

SRM_oes <- SRM_values %>% mutate(conc_mg_g = partial_extraction_concentration.ug.g. / 1000)

#pivot back and specify IDs to stop NA filling
srms_wide <- srms_long %>% pivot_wider(id_cols = "SampleID", names_from = element, values_from = DRY_CONC)

write.csv(srms_wide,file.path("./2_incremental","srm_data.csv"))



