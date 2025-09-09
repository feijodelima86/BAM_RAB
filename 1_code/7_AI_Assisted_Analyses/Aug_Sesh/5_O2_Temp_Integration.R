# Load required libraries
library(dplyr)
library(lubridate)

# Read the main dataset
main_data <- integrated_data
main_data$Date <- as.Date(main_data$Date)

# Get unique sites
sites <- unique(main_data$Site)

# Initialize empty dataframe to store daily averages
daily_averages <- data.frame()

# Process each site
for(site in sites) {
  
  # Option 1: All sensor files in a subdirectory
  filename <- file.path("2_incremental/AI_Subsets/Minidot", paste0(site, ".csv"))
  

  if(file.exists(filename)) {
    
    # Read the sensor data
    sensor_data <- read.csv(filename)
    
    # Convert datetime to proper format and extract date
    sensor_data$datetime <- as.POSIXct(sensor_data$UTCdateTimeSensor, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    sensor_data$date <- as.Date(sensor_data$datetime)
    
    # Calculate daily averages
    site_daily <- sensor_data %>%
      group_by(date) %>%
      summarise(
        daily_temp = mean(waterTempSensor, na.rm = TRUE),
        daily_DO = mean(DOobs, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Add site column
    site_daily$Site <- site
    
    # Append to daily averages
    daily_averages <- rbind(daily_averages, site_daily)
  }
}

# Calculate accumulated degree days for each site (after all sites processed)
base_temp <- 0  # Base temperature for degree day calculation (adjust as needed)

# Check if daily_averages has data
print(paste("Number of rows in daily_averages:", nrow(daily_averages)))
print("Column names in daily_averages:")
print(colnames(daily_averages))
print("First few rows:")
print(head(daily_averages))

# Only proceed if we have data
if(nrow(daily_averages) > 0) {
  daily_averages <- daily_averages %>%
    arrange(Site, date) %>%
    group_by(Site) %>%
    mutate(
      degree_days = pmax(daily_temp - base_temp, 0),  # Daily degree days above base temp
      accumulated_degree_days = cumsum(degree_days)   # Cumulative sum
    ) %>%
    ungroup()
} else {
  print("No sensor data found! Check file paths and filenames.")
}

# Merge with main dataset by finding closest date match
main_data$daily_temp <- NA
main_data$daily_DO <- NA
main_data$accumulated_degree_days <- NA

for(i in 1:nrow(main_data)) {
  
  site_name <- main_data$Site[i]
  sample_date <- main_data$Date[i]
  
  # Get daily averages for this site
  site_averages <- daily_averages[daily_averages$Site == site_name, ]
  
  if(nrow(site_averages) > 0) {
    
    # Calculate date differences
    date_diffs <- abs(as.numeric(site_averages$date - sample_date))
    
    # Find closest date
    closest_idx <- which.min(date_diffs)
    
    # Assign values
    main_data$daily_temp[i] <- site_averages$daily_temp[closest_idx]
    main_data$daily_DO[i] <- site_averages$daily_DO[closest_idx]
    main_data$accumulated_degree_days[i] <- site_averages$accumulated_degree_days[closest_idx]
  }
}

# View the merged dataset
head(main_data)