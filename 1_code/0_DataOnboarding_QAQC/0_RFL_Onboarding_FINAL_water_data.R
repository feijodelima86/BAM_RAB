library(readr)
library(tidyr)
library(tidyverse)
library(lubridate)

water <- read_csv("2_incremental/230731_Water_Metals_long_diff.csv")

water_seperate <- water %>%
  # Separate the date, size, and site using str_extract 
  mutate(
    date = str_extract(Sample.Id, "\\d+_\\d+_\\d+"),             # Extract date in the format: 6_21_21
    filtersize = str_extract(Sample.Id, "(W|UF|F)"),                   # Extract size (W, UF, or F)
    site = str_extract(Sample.Id, "(?<=_(W|UF|F)_)(\\d+)")       # Extract the site using lookbehind regex
  ) %>%
  # Convert date to the actual date type using lubridate package
  mutate(date = mdy(date))

water_calc <- water_seperate %>%
  pivot_wider(id_cols = c(site, date,Element),values_from = Concentration, names_from = filtersize)

water_calc$site <- factor(water_calc$site , labels = c("WS","DL","GR","GC","BG","BN"))

test <- water_calc %>% pivot_wider(names_from = c("Element"), values_from = c("UF", "F", "W"))

write.csv(test, "2_incremental/Final_water_data.csv")
