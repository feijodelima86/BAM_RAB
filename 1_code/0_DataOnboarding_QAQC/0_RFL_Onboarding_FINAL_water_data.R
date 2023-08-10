library(readr)
library(tidyr)

water <- read_csv("2_incremental/230801_Water_Metals_Size.csv")

names(water)

test <- water %>% pivot_wider(names_from = Element, values_from = c("UF","F","W","spm","colloidal","trulydissolved"))

write.csv(test, "2_incremental/Final_water_data.csv")
