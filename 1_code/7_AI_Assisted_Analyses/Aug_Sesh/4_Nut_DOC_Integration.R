# Read your main dataset
main_data <- speciation_spatial

# Read nutrient data - try different approaches
# Option 1: If it's a space/tab separated file
nutrient_data <- read.csv("2_incremental/AI_Subsets/AnnualAggregateFile_NH4_N_SRP_NO3_N_WY2021.csv")

# Option 2: If the above fails, create from your sample data
# (Replace this with actual file reading once format is confirmed)
if (nrow(nutrient_data) == 0) {
  # Sample data you provided - adjust this to read your actual file
  cat("Nutrient data file not found or empty. Check file name and format.\n")
  cat("Expected columns: dateSample, calendarYearSample, calendarMonthSample, calendarDaySample, siteID, repID, NH4N, SRP, NO3N\n")
  cat("Try: read.table('your_file.txt', header=TRUE) or read.csv('your_file.csv')\n")
}

# Check what siteIDs exist in nutrient data
unique_sites <- unique(nutrient_data$siteID)
print("Unique siteIDs in nutrient data:")
print(unique_sites)

# Create site mapping using character format with leading zeros
# From your description: 1=WS, 6=DL, 8=GR, 9=GC, 10=BG, 11=BN
site_mapping <- data.frame(
  siteID = c("001", "006", "008", "009", "010", "011"),
  Site = c("WS", "DL", "GR", "GC", "BG", "BN"),
  stringsAsFactors = FALSE
)

# Only keep siteIDs that have mapping and exclude NAs
nutrient_data <- nutrient_data[!is.na(nutrient_data$siteID) & 
                                 nutrient_data$siteID %in% site_mapping$siteID, ]

print(paste("Filtered to", nrow(nutrient_data), "nutrient records"))

# Convert dates in main data to Date format
main_data$Date <- as.Date(main_data$Date, format = "%m/%d/%Y")

# Convert dates in nutrient data to Date format  
nutrient_data$Date <- as.Date(as.character(nutrient_data$dateSample), format = "%Y%m%d")

# Add site abbreviations to nutrient data
nutrient_data <- merge(nutrient_data, site_mapping, by = "siteID")

# Aggregate nutrient data by site and date (mean of replicates)
nutrient_summary <- aggregate(cbind(NH4N, SRP, NO3N) ~ Site + Date, 
                              data = nutrient_data, 
                              FUN = mean, na.rm = TRUE)

# Create function to find closest date for any dataset
find_closest_match <- function(algae_site, algae_date, reference_data) {
  site_data <- reference_data[reference_data$Site == algae_site, ]
  if (nrow(site_data) == 0) return(NA)
  
  time_diffs <- abs(as.numeric(site_data$Date - algae_date))
  closest_idx <- which.min(time_diffs)
  
  return(site_data[closest_idx, ])
}

# Apply to each row of main data for nutrients
integrated_data <- main_data
nutrient_matches <- do.call(rbind, lapply(1:nrow(main_data), function(i) {
  find_closest_match(main_data$Site[i], main_data$Date[i], nutrient_summary)
}))

# Add nutrient columns to main data
integrated_data$NH4N <- nutrient_matches$NH4N
integrated_data$SRP <- nutrient_matches$SRP  
integrated_data$NO3N <- nutrient_matches$NO3N
integrated_data$nutrient_date <- nutrient_matches$Date
integrated_data$nutrient_days_diff <- as.numeric(abs(integrated_data$Date - integrated_data$nutrient_date))

# Read DOC data 
doc_data <- read.csv("2_incremental/AI_Subsets/table_doc_wy2021.csv")

# Convert DOC dates to Date format
doc_data$Date <- as.Date(as.character(doc_data$dateString), format = "%Y%m%d")

# Create DOC site mapping (numeric siteID format)
doc_site_mapping <- data.frame(
  siteID = c(1, 6, 8, 9, 10, 11),
  Site = c("WS", "DL", "GR", "GC", "BG", "BN")
)

# Add site abbreviations to DOC data
doc_data <- merge(doc_data, doc_site_mapping, by = "siteID")

# Aggregate DOC data by site and date (mean of replicates)
doc_summary <- aggregate(doc ~ Site + Date, 
                         data = doc_data, 
                         FUN = mean, na.rm = TRUE)

# Initialize integrated data
integrated_data <- main_data

# Add empty columns for nutrients
integrated_data$NH4N <- NA
integrated_data$SRP <- NA
integrated_data$NO3N <- NA
integrated_data$nutrient_date <- as.Date(NA)
integrated_data$nutrient_days_diff <- NA

# Match nutrients by temporal proximity
for (i in 1:nrow(integrated_data)) {
  site_nutrients <- nutrient_summary[nutrient_summary$Site == integrated_data$Site[i], ]
  
  if (nrow(site_nutrients) > 0) {
    time_diffs <- abs(as.numeric(site_nutrients$Date - integrated_data$Date[i]))
    closest_idx <- which.min(time_diffs)
    
    integrated_data$NH4N[i] <- site_nutrients$NH4N[closest_idx]
    integrated_data$SRP[i] <- site_nutrients$SRP[closest_idx]
    integrated_data$NO3N[i] <- site_nutrients$NO3N[closest_idx]
    integrated_data$nutrient_date[i] <- site_nutrients$Date[closest_idx]
    integrated_data$nutrient_days_diff[i] <- time_diffs[closest_idx]
  }
}

# Add empty columns for DOC
integrated_data$DOC <- NA
integrated_data$doc_date <- as.Date(NA)
integrated_data$doc_days_diff <- NA

# Match DOC by temporal proximity
for (i in 1:nrow(integrated_data)) {
  site_doc <- doc_summary[doc_summary$Site == integrated_data$Site[i], ]
  
  if (nrow(site_doc) > 0) {
    time_diffs <- abs(as.numeric(site_doc$Date - integrated_data$Date[i]))
    closest_idx <- which.min(time_diffs)
    
    integrated_data$DOC[i] <- site_doc$doc[closest_idx]
    integrated_data$doc_date[i] <- site_doc$Date[closest_idx]
    integrated_data$doc_days_diff[i] <- time_diffs[closest_idx]
  }
}

# View result
head(integrated_data)