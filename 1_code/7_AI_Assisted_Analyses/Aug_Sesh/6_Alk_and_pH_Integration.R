# Load required libraries
library(dplyr)
library(lubridate)

# Read main dataset
main_data <- main_data
main_data$Date <- as.Date(main_data$Date)

# Create alkalinity/pH dataset
alkalinity_data <- data.frame(
  Site_Code = c(1, 6, 8, 9, 10, 11, 1, 6, 8, 9, 10, 11, 1, 6, 8, 9, 10, 11, 
                1, 6, 8, 9, 10, 11, 1, 6, 8, 9, 10, 11),
  Date = as.Date(c("2021-06-21", "2021-06-21", "2021-06-21", "2021-06-21", "2021-06-21", "2021-06-21",
                   "2021-07-19", "2021-07-19", "2021-07-19", "2021-07-19", "2021-07-19", "2021-07-19",
                   "2021-08-30", "2021-08-30", "2021-08-30", "2021-08-30", "2021-08-30", "2021-08-30",
                   "2021-09-15", "2021-09-15", "2021-09-15", "2021-09-15", "2021-09-15", "2021-09-15",
                   "2021-10-15", "2021-10-15", "2021-10-15", "2021-10-15", "2021-10-15", "2021-10-15")),
  Alkalinity = c(1860.8, 2481, 2941.6, 2797.5, 3176.5, 3155.7, 2272, 3181.2, 3412.7, 2914.2,
                 3294.3, 3377.5, 2043.8, 2765.1, 2937.8, 2879.2, 3483.9, 3265, 2140.2, 3103.3,
                 3182.7, 3063.4, 3433.2, 3417, 2244.6, 2960.4, 3050, 3029.6, 3349, 3354.5),
  pH = c(8.1979, 8.5148, 8.5262, 8.4311, 8.5709, 8.6512, 8.4468, 8.2501, 8.3596, 8.5139,
         8.2504, 8.4343, 8.6386, 8.3873, 8.5688, 8.6163, 8.5614, 8.6815, 8.3098, 8.272,
         8.4487, 8.5782, 8.5373, 8.5477, 8.3531, 8.331, 8.787, 8.576, 8.4428, 8.5168)
)

# Site code mapping
site_mapping <- data.frame(
  Site_Code = c(1, 6, 8, 9, 10, 11),
  Site = c("WS", "DL", "GR", "GC", "BG", "BN")
)

# Apply mapping to alkalinity data
alkalinity_data <- alkalinity_data %>%
  left_join(site_mapping, by = "Site_Code") %>%
  select(Site, Date, Alkalinity, pH)

# Initialize result columns in main_data
main_data$Alkalinity <- NA
main_data$pH <- NA
main_data$alkalinity_date <- as.Date(NA)
main_data$alkalinity_days_diff <- NA

# Perform temporal matching
for (i in 1:nrow(main_data)) {
  current_site <- main_data$Site[i]
  current_date <- main_data$Date[i]
  
  # Filter alkalinity data for same site
  site_alkalinity <- alkalinity_data[alkalinity_data$Site == current_site, ]
  
  if (nrow(site_alkalinity) > 0) {
    # Calculate days difference
    days_diff <- abs(as.numeric(site_alkalinity$Date - current_date))
    
    # Find closest match
    closest_idx <- which.min(days_diff)
    
    # Assign values
    main_data$Alkalinity[i] <- site_alkalinity$Alkalinity[closest_idx]
    main_data$pH[i] <- site_alkalinity$pH[closest_idx]
    main_data$alkalinity_date[i] <- site_alkalinity$Date[closest_idx]
    main_data$alkalinity_days_diff[i] <- days_diff[closest_idx]
  }
}

# Final integrated dataset
integrated_data <- main_data

write.csv(main_data, "2_incremental/AI_Subsets/integrated_parsed.csv")
