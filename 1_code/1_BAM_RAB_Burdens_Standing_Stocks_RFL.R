library(readr)

BB.RAW <- read.csv("2_incremental/20220307_FIELD_SAMPLES.csv")

#Creating dataset with relevant columns and correcting for dilution (50ml) to obtain g/sample

names(BB.RAW)

BB.PRUNE <- cbind(BB.RAW[,c(42,43,44,45,53)], BB.RAW[,c(3:41)]*50/1000)

names(BB.PRUNE)

BB.FINAL <- cbind(BB.PRUNE[,c(1,2,3,4)],(BB.PRUNE[,c(6:44)]/BB.PRUNE[,5]))

write.csv(BB.FINAL, paste0("2_incremental/",gsub("-", "", Sys.Date()),"_BB_CURRENT.csv"))

