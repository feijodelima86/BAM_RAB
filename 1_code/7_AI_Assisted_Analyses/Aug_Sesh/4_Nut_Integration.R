# Read your main dataset
main_data <- speciation_spatial

# Read nutrient data (assuming it's in a file called nutrients.csv)
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

# Create function to find closest nutrient date for each algae sample
find_closest_nutrient <- function(algae_site, algae_date, nutrient_summary) {
  site_nutrients <- nutrient_summary[nutrient_summary$Site == algae_site, ]
  if (nrow(site_nutrients) == 0) return(NA)
  
  time_diffs <- abs(as.numeric(site_nutrients$Date - algae_date))
  closest_idx <- which.min(time_diffs)
  
  return(site_nutrients[closest_idx, ])
}

# Apply to each row of main data
integrated_data <- main_data
nutrient_matches <- do.call(rbind, lapply(1:nrow(main_data), function(i) {
  find_closest_nutrient(main_data$Site[i], main_data$Date[i], nutrient_summary)
}))

# Add nutrient columns to main data
integrated_data$NH4N <- nutrient_matches$NH4N
integrated_data$SRP <- nutrient_matches$SRP  
integrated_data$NO3N <- nutrient_matches$NO3N
integrated_data$nutrient_date <- nutrient_matches$Date
integrated_data$days_diff <- as.numeric(abs(integrated_data$Date - integrated_data$nutrient_date))

# View result
head(integrated_data)
